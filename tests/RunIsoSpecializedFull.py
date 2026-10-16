"""Regress all three fixed specialized suites; never call a retired validator.

Three fresh-context runs per input, serial subprocesses, no per-case warmup.
The process watchdog includes SDK initialization/source certification and is
not an algorithm deadline. Geometric failures return 1; execution errors 2.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
from datetime import datetime, timezone

STAGES = ['preparation_ms', 'position_ms', 'subdivision_ms', 'boundary_ms', 'reporting_ms']
CLASS = ['class_sampling_ms', 'class_algebraic_ms', 'class_geometric_ms',
         'class_scoring_ms', 'class_recovery_ms', 'class_other_ms']
METRICS = ['classification_compute_ms', 'classification_validation_ms', 'algorithm_ms'] + STAGES + CLASS
ATTEMPT_FIELDS = ['suite', 'id', 'run', 'candidates', 'output_type', 'status', 'accepted',
                  'classification_stages', 'classification_compute_ms', 'classification_validation_ms'] + STAGES + [
                  'source_cache_hits', 'mother_bound', 'edge_bound', 'joint_tangent_bound_deg',
                  'joint_tangent_failure', 'message']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_csv(p):
    if not p.exists():
        return []
    with p.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def save_csv(p, rows):
    if not rows:
        p.write_text('', encoding='utf-8')
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def means(rows):
    out = {}
    for k in METRICS:
        if k in STAGES + CLASS:
            # Full-case costs: a polynomial case that never enters validation
            # contributes zero work. Rational/mixed populations are N/A.
            out[k] = statistics.mean(float(r.get(k) or 0) for r in rows) if rows and all(
                r['backend'] == 'mother' for r in rows) else ''
        else:
            out[k] = statistics.mean(float(r[k]) for r in rows if r.get(k)) if any(r.get(k) for r in rows) else ''
    return out


def check(rows, attempts, case, runs):
    assert len(rows) == runs and {int(r['run']) for r in rows} == set(range(runs)), 'missing/duplicate runs'
    for r in rows:
        assert int(r['id']) == case
        for k in METRICS:
            if r.get(k):
                assert math.isfinite(float(r[k])) and float(r[k]) >= -1e-7, ('invalid timing', k, r[k])
        assert abs(float(r['classification_compute_ms']) + float(r['classification_validation_ms']) - float(r['algorithm_ms'])) < 1e-7
        assert r['optimization_ran'] == '0'
        if r['backend'] == 'mother':
            selected = [a for a in attempts if a['run'] == r['run']]
            assert len(selected) == int(r['classification_validation_calls'])
            assert abs(sum(float(a['classification_validation_ms']) for a in selected) - float(r['classification_validation_ms'])) < 1e-7
            for k in STAGES:
                assert abs(sum(float(a[k]) for a in selected) - float(r[k] or 0)) < 1e-7
            assert abs(sum(float(r[k] or 0) for k in STAGES) - float(r['classification_validation_ms'])) < 1e-7
            assert abs(sum(float(r[k]) for k in CLASS) - float(r['classification_compute_ms'])) < 1e-7
        else:
            assert all(not r[k] for k in STAGES), 'rational five stages must be N/A'
        if r['accepted'] == '1':
            assert r['audit'] in ('Passed', 'Failed')
            if r['backend'] == 'mother':
                assert float(r['mother_bound']) < 1e-5 and float(r['edge_bound']) < 1e-5
                assert float(r['joint_tangent_bound_deg']) < .1
                assert float(r['joint_sine_squared_bound']) < math.sin(math.radians(.1)) ** 2
                assert not r['joint_tangent_failure']


def summarize(rows, records, output):
    cases = []
    for record in records:
        group = [r for r in rows if r['suite'] == record['suite'] and int(r['id']) == record['id']]
        times = [float(r['algorithm_ms']) for r in group]
        first = group[0] if group else {}
        cases.append(dict(suite=record['suite'], id=record['id'], group=record['group'],
            type=first.get('type', ''), backend=first.get('backend', ''), collapsed=first.get('collapsed', ''),
            degree_u=first.get('degree_u', ''), degree_v=first.get('degree_v', ''), controls=first.get('controls', ''),
            runs=len(group), accepted=sum(r['accepted'] == '1' for r in group),
            audited=sum(r['audit'] == 'Passed' for r in group),
            expected=sum(r['expected_outcome'] == '1' for r in group),
            qualified=first.get('qualified', ''), statuses=' | '.join(sorted({r['status'] for r in group})),
            output_types=' | '.join(sorted({r['output_type'] for r in group if r['output_type']})),
            same_type=all(r['same_type'] == '1' for r in group) if group else '',
            process_error=record.get('error', ''), **means(group),
            total_min_ms=min(times) if times else '', total_max_ms=max(times) if times else '',
            runs_under_10ms=sum(t < 10 for t in times),
            validation_attempts_max=max((int(r['classification_validation_calls']) for r in group), default=0),
            messages=' | '.join(sorted({r['message'] for r in group}))))
    groups = []
    for suite, backend, kind in sorted({(r['suite'], r['backend'], r['type']) for r in rows}):
        for population in ('all', 'accepted_audited'):
            g = [r for r in rows if (r['suite'], r['backend'], r['type']) == (suite, backend, kind)
                 and (population == 'all' or r['accepted'] == '1' and r['audit'] == 'Passed')]
            if g:
                groups.append(dict(suite=suite, backend=backend, type=kind, population=population,
                    cases=len({r['id'] for r in g}), measurements=len(g),
                    validation_measurements=sum(r['validation_ran'] == '1' for r in g), **means(g),
                    max_ms=max(float(r['algorithm_ms']) for r in g),
                    runs_under_10ms=sum(float(r['algorithm_ms']) < 10 for r in g)))
    save_csv(output/'per-case.csv', cases)
    save_csv(output/'modules-by-type.csv', groups)
    save_csv(output/'failures.csv', [c for c in cases if c['process_error'] or c['expected'] != c['runs']])
    save_csv(output/'over-10ms.csv', [r for r in rows if float(r['algorithm_ms']) >= 10])
    summary = {}
    for suite in sorted({r['suite'] for r in records}):
        cs = [c for c in cases if c['suite'] == suite]
        rs = [r for r in rows if r['suite'] == suite]
        summary[suite] = dict(cases=len(cs), measurements=len(rs),
            all_runs_accepted_audited=sum(c['runs'] > 0 and c['audited'] == c['runs'] and not c['process_error'] for c in cs),
            all_runs_expected=sum(c['runs'] > 0 and c['expected'] == c['runs'] and not c['process_error'] for c in cs),
            process_errors=sum(bool(c['process_error']) for c in cs),
            audit_failures=sum(r['audit'] == 'Failed' for r in rs),
            runs_under_10ms=sum(float(r['algorithm_ms']) < 10 for r in rs),
            max_ms=max((float(r['algorithm_ms']) for r in rs), default=None), **means(rs))
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    return cases, summary


def main():
    parser = argparse.ArgumentParser(__doc__)
    for key in ('exe', 'sdk', 'work', 'output'):
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--timeout', type=float, default=120)
    parser.add_argument('--suites', nargs='+', choices=['dataset500', 'business-extreme200', 'extreme200'],
                        default=['dataset500', 'business-extreme200', 'extreme200'])
    parser.add_argument('--cases', nargs='+', type=int, help='Optional smoke-test subset; omitted for full regression')
    args = parser.parse_args()
    assert 1 <= args.runs <= 100 and args.timeout > 0
    root = Path(__file__).resolve().parent.parent
    for k in ('exe', 'sdk', 'work', 'output'):
        setattr(args, k, getattr(args, k).resolve())
    args.work.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    sources = sorted(list((root/'src').glob('*.[ch]*')) + list((root/'include').glob('*.h')))
    sources += [root/'CMakeLists.txt', Path(__file__).resolve()] + [root/'tests'/f for f in (
        'MotherDistanceDatasetTests.cpp', 'Dataset500.cpp', 'ExtremeDataset.cpp', 'BusinessExtremeDataset.cpp', 'PairedSurfaceBounds.cpp')]
    manifest = dict(started_utc=datetime.now(timezone.utc).isoformat(), runs=args.runs, timeout_seconds=args.timeout,
        source_commit=subprocess.check_output(['git', '-c', 'safe.directory='+root.as_posix(), 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
        source_sha256={p.relative_to(root).as_posix(): sha(p) for p in sources}, executable_sha256=sha(args.exe),
        contract='incremental isoparametric classification + mother distance / exact Recognize; no optimizer; no retired validators',
        records=[])
    env = os.environ.copy()
    env['PATH'] = str(args.sdk/'x64-win/bin') + os.pathsep + env['PATH']
    rows, attempts = [], []
    (args.output/'runs.log').write_text('', encoding='utf-8')
    for suite in args.suites:
        data = root/'tests/data'/suite
        count = 500 if suite == 'dataset500' else 200
        prefix = 'case' if suite == 'dataset500' else 'extreme'
        assert len(list(data.glob(prefix+'-*.b2a'))) == count
        metadata = {int(r['id']): r for r in read_csv(data/'manifest.csv')}
        for case in args.cases if args.cases is not None else range(count):
            assert 0 <= case < count
            dest = args.work/f'{suite}-{case:03}.csv'
            input_file = data/f'{prefix}-{case:03}.b2a'
            command = [str(args.exe), '--case', str(data), str(case), str(dest), str(args.runs)]
            record = dict(suite=suite, id=case, group=metadata.get(case, {}).get('group', 'standard'),
                          input_sha256=sha(input_file), command=command)
            log = ''
            try:
                result = subprocess.run(command, cwd=args.sdk, env=env, capture_output=True, text=True,
                                        encoding='utf-8', errors='replace', timeout=args.timeout)
                log = result.stdout + result.stderr
                record['returncode'] = result.returncode
                assert result.returncode in (0, 2), f'process exit {result.returncode}'
                current, parts = read_csv(dest), read_csv(Path(str(dest)+'.attempts.csv'))
                check(current, parts, case, args.runs)
            except (subprocess.TimeoutExpired, AssertionError, ValueError, KeyError) as e:
                record['error'] = str(e)
                if isinstance(e, subprocess.TimeoutExpired):
                    log = (e.stdout or b'').decode('utf-8', errors='replace') + (e.stderr or b'').decode('utf-8', errors='replace')
                print(f'PROCESS_ERROR {suite}/{case}: {e}', flush=True)
            finally:
                current, parts = read_csv(dest), read_csv(Path(str(dest)+'.attempts.csv'))
                rows.extend(dict(suite=suite, group=record['group'], **r) for r in current)
                attempts.extend({k: (suite if k == 'suite' else r.get(k, '')) for k in ATTEMPT_FIELDS} for r in parts)
                with (args.output/'runs.log').open('a', encoding='utf-8') as stream:
                    stream.write(f'[{suite}/{case:03}]\n'+log.rstrip()+'\n')
                manifest['records'].append(record)
            if (case+1) % 25 == 0 or args.cases:
                print(f'{suite}: {case+1}/{count}; measurements={len(rows)}; errors={sum(bool(r.get("error")) for r in manifest["records"])}', flush=True)
                save_csv(args.output/'measurements.csv', rows)
                save_csv(args.output/'attempts.csv', attempts)
                (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    assert all(sha(root/k) == v for k, v in manifest['source_sha256'].items()), 'sources changed during regression'
    save_csv(args.output/'measurements.csv', rows)
    save_csv(args.output/'attempts.csv', attempts)
    cases, summary = summarize(rows, manifest['records'], args.output)
    manifest['completed_utc'] = datetime.now(timezone.utc).isoformat()
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)
    if any(r.get('error') for r in manifest['records']):
        return 2
    return int(any(c['expected'] != c['runs'] for c in cases))


if __name__ == '__main__':
    raise SystemExit(main())
