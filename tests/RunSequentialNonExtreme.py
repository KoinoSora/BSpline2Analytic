"""Serial cold-context iso-then-fallback regression. IO/audits are outside total_ms."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess


def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write(p, rows):
    if not rows:
        p.write_text('', encoding='utf-8')
        return
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('exe', 'sdk', 'data', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--runs', type=int, default=3)
    p.add_argument('--cases', nargs='*')
    p.add_argument('--resume', action='store_true')
    a = p.parse_args()
    exe, sdk, data, output = (getattr(a, k).resolve() for k in ('exe', 'sdk', 'data', 'output'))
    if not 1 <= a.runs <= 100:
        p.error('runs must be 1..100')
    from RationalDatasetSupport import require_v2_data
    require_v2_data(data/'general-rational-v2')
    records = []
    for suite, count in [('iso', 500), ('main', 600), ('extra', 60)]:
        folder = data/'dataset500' if suite == 'iso' else data/'general-rational-v2'/suite
        manifest = read(folder/'manifest.csv')
        assert [int(r['id']) for r in manifest] == list(range(count))
        for r in manifest:
            i = int(r['id'])
            path = folder/(f'case-{i:03}.b2a' if suite == 'iso' else r['file'])
            expectation = r.get('expectation', 'accept')
            if suite == 'iso' and int(r['collapsed_mask']):
                expectation = 'collapsed_probe'
            category = 'positive' if expectation == 'accept' else 'negative' if expectation == 'reject_infeasible' else 'diagnostic'
            records.append(dict(suite=suite, id=i, expectation=expectation, category=category,
                                sha256=sha(path), group=r.get('group', r.get('type', ''))))
    if a.cases:
        selected = set(a.cases)
        records = [r for r in records if f"{r['suite']}:{r['id']}" in selected]
        assert len(records) == len(selected), 'Unknown/duplicate case key'
    output.mkdir(parents=True, exist_ok=True)
    config = dict(exe_sha256=sha(exe), runs=a.runs, max_case_ms=50, cases=records,
                  options='default; no supplied truth/chart or per-fixture option injection',
                  timing='steady clock around entire API; IO, source certificates and output audits excluded',
                  repetitions='fresh context, no per-case warmup', parallel_processes=1)
    stamp = output/'run.json'
    if stamp.exists():
        assert a.resume and json.loads(stamp.read_text()) == config, 'Use a new output directory'
    stamp.write_text(json.dumps(config, indent=2), encoding='utf-8')
    parts = output/'cases'
    parts.mkdir(exist_ok=True)
    env = dict(os.environ)
    key = next(k for k in env if k.lower() == 'path')
    env[key] = str(sdk/'x64-win/bin')+os.pathsep+env[key]
    rows, errors = [], []
    for n, r in enumerate(records):
        target = parts/f"{r['suite']}-{r['id']:03}.csv"
        try:
            if not (a.resume and target.exists() and len(read(target)) == a.runs):
                with target.with_suffix('.log').open('wb') as log:
                    result = subprocess.run([str(exe), str(data), r['suite'], str(r['id']), str(a.runs), str(target), 'run'],
                                            cwd=sdk, env=env, stdout=log, stderr=subprocess.STDOUT,
                                            timeout=120, creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
                assert result.returncode == 0, f'exit {result.returncode}'
            case = read(target)
            assert len(case) == a.runs and {int(x['run']) for x in case} == set(range(a.runs))
            for x in case:
                assert x['suite'] == r['suite'] and int(x['id']) == r['id']
                assert x['source_certified'] == x['iso_attempted'] == '1'
                ms = float(x['total_ms'])
                assert math.isfinite(ms) and ms > 0 and ms >= float(x['api_ms'])
                assert (ms < 10) == (x['under_10ms'] == '1')
                assert (ms <= 50) == (x['within_50ms'] == '1')
                if ms > 50 or x['time_budget_exceeded'] == '1':
                    assert x['accepted'] == '0'
                if x['accepted'] == '0':
                    assert x['fallback_attempted'] == '1' or x['time_budget_exceeded'] == '1'
                elif x['path'] == 'general_rational':
                    assert x['fallback_attempted'] == '1'
                else:
                    assert x['fallback_attempted'] == '0'
                x.update(expectation=r['expectation'], category=r['category'], group=r['group'])
            rows.extend(case)
        except (AssertionError, OSError, subprocess.TimeoutExpired) as e:
            errors.append(dict(suite=r['suite'], id=r['id'], error=str(e)))
        if (n+1) % 50 == 0 or n+1 == len(records):
            print(f'{n+1}/{len(records)} cases; {len(errors)} execution errors', flush=True)
    write(output/'runs.csv', rows)
    write(output/'execution-errors.csv', errors)
    cases = []
    for r in records:
        rr = [x for x in rows if x['suite'] == r['suite'] and int(x['id']) == r['id']]
        if len(rr) != a.runs:
            continue
        times = [float(x['total_ms']) for x in rr]
        accepted = all(x['accepted'] == '1' for x in rr)
        audit = accepted and all(x['audit'] == 'passed' for x in rr)
        timeout = any(x['time_budget_exceeded'] == '1' for x in rr)
        false_accept = r['category'] == 'negative' and any(x['accepted'] == '1' for x in rr)
        audit_fail = any(x['audit'] == 'failed' for x in rr)
        correct = not timeout and not audit_fail and (audit if r['category'] == 'positive' else
                  not false_accept if r['category'] == 'negative' else False)
        cases.append(dict(suite=r['suite'], id=r['id'], category=r['category'], expectation=r['expectation'],
                          group=r['group'], accepted=accepted, audited_success=audit, correct=correct,
                          false_accept=false_accept, audit_failed=audit_fail, timeout=timeout,
                          fallback=any(x['fallback_attempted'] == '1' for x in rr),
                          path='|'.join(sorted({x['path'] for x in rr})),
                          min_ms=min(times), mean_ms=statistics.mean(times), max_ms=max(times),
                          under_10ms=max(times) < 10, within_50ms=max(times) <= 50))
    write(output/'cases.csv', cases)
    write(output/'failures.csv', [x for x in cases if (x['category'] != 'diagnostic' and not x['correct']) or x['audit_failed'] or x['timeout']])
    write(output/'over-10ms.csv', [x for x in cases if not x['under_10ms']])
    write(output/'over-50ms.csv', [x for x in rows if x['within_50ms'] == '0'])
    write(output/'timeouts.csv', [x for x in rows if x['time_budget_exceeded'] == '1'])
    summary = []
    for suite in ('iso', 'main', 'extra', 'all'):
        cc = [x for x in cases if suite == 'all' or x['suite'] == suite]
        rr = [x for x in rows if suite == 'all' or x['suite'] == suite]
        if not cc:
            continue
        summary.append(dict(suite=suite, cases=len(cc), runs=len(rr),
            positive=sum(x['category'] == 'positive' for x in cc),
            positive_correct=sum(x['category'] == 'positive' and x['correct'] for x in cc),
            negative=sum(x['category'] == 'negative' for x in cc),
            negative_correct=sum(x['category'] == 'negative' and x['correct'] for x in cc),
            diagnostic=sum(x['category'] == 'diagnostic' for x in cc),
            audited_success=sum(x['audited_success'] for x in cc), fallback=sum(x['fallback'] for x in cc),
            timeouts=sum(x['timeout'] for x in cc), false_accepts=sum(x['false_accept'] for x in cc),
            audit_failures=sum(x['audit_failed'] for x in cc),
            under_10ms=sum(x['under_10ms'] for x in cc), within_50ms=sum(x['within_50ms'] for x in cc),
            correct_under_10ms=sum(x['correct'] and x['under_10ms'] for x in cc),
            total_algorithm_ms=sum(float(x['total_ms']) for x in rr),
            mean_ms=statistics.mean(float(x['total_ms']) for x in rr),
            max_ms=max(x['max_ms'] for x in cc)))
    write(output/'summary.csv', summary)
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)
    return 2 if errors else int(any((x['category'] != 'diagnostic' and not x['correct']) or x['audit_failed'] or x['timeout'] for x in cases))


if __name__ == '__main__':
    raise SystemExit(main())
