"""Serial paired measurements of two CURRENT general rational verifiers.

Alternate backend order per case; no retired validator is invoked. Both
executables run the same existing full-case driver, with certificates and SDK
audits outside the algorithm timer. Preserve every raw measurement.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess

CASES = [('main', i) for i in (0, 98, 100, 198, 218, 238, 298, 318, 398, 438, 498,
                               104, 244, 254, 304, 444, 500)]
CASES += [('synthetic', i) for i in range(5)]
CASES += [('extreme', i) for i in (40, 152, 153, 156, 157, 246, 254, 255, 314)]
CASES += [('extra', i) for i in (40, 50, 55)]


def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write(p, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('baseline', 'updated', 'sdk', 'data', 'extreme-data', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--runs', type=int, default=5)
    args = p.parse_args()
    assert args.runs > 0
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ)
    key = next(k for k in env if k.lower() == 'path')
    sdk = args.sdk.resolve()
    env[key] = str(sdk/'x64-win/bin')+os.pathsep+env[key]
    rows, attempts, summaries = [], [], []
    manifest = dict(runs=args.runs, cases=CASES, baseline_sha256=sha(args.baseline),
                    updated_sha256=sha(args.updated), inputs={}, retired_validators_run=False,
                    parallel_processes=1, cpu_affinity_fixed=False, exclusive_machine=False)
    for index, (suite, case) in enumerate(CASES):
        data = (args.extreme_data if suite == 'extreme' else args.data).resolve()
        if suite != 'synthetic':
            name = {'main': 'manifest-general600.csv', 'extra': 'manifest-general-extra60.csv',
                    'extreme': 'manifest.csv'}[suite]
            meta = next(x for x in read(data/name) if int(x['id']) == case)
            manifest['inputs'][f'{suite}:{case}'] = sha(data/meta['file'])
        pair = {}
        for backend in (('baseline', 'updated') if index % 2 == 0 else ('updated', 'baseline')):
            target = output/f'{backend}-{suite}-{case}.csv'
            with target.with_suffix('.log').open('wb') as log:
                subprocess.run([str(getattr(args, backend).resolve()), '--rational-full-case',
                                str(data), str(target), suite, str(case), str(args.runs), '--rational-only'],
                               cwd=sdk, env=env, stdout=log, stderr=subprocess.STDOUT,
                               timeout=120, check=True, creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
            part = read(target)
            assert len(part) == args.runs and [int(x['run']) for x in part] == list(range(args.runs))
            for r in part:
                assert r['suite'] == suite and int(r['id']) == case
                assert r['accepted'] != '1' or r['audit_status'] == 'passed', (suite, case, backend, 'audit')
                if r['accepted'] == '1':
                    assert float(r['boundary_bound']) < 1e-5 and float(r['tangent_deg']) < .1
                r['backend'] = backend
            rows.extend(part)
            for r in read(target.with_name(target.stem+'-attempts.csv')):
                r['backend'] = backend
                attempts.append(r)
            pair[backend] = part
        for a, b in zip(pair['baseline'], pair['updated']):
            for field in ('status', 'accepted', 'output_type', 'attempts', 'chart_attempts'):
                assert a[field] == b[field], (suite, case, field, a[field], b[field])
        summary = dict(suite=suite, id=case, accepted=pair['updated'][0]['accepted'],
                       degree_u=part[0]['degree_u'], degree_v=part[0]['degree_v'])
        for backend, part in pair.items():
            for field in ('total_ms', 'boundary_ms', 'prepare_ms', 'nodes', 'refinements'):
                summary[backend+'_'+field] = statistics.median(float(r[field]) for r in part)
            summary[backend+'_max_ms'] = max(float(r['total_ms']) for r in part)
        for field in ('envelope_checks', 'envelope_accepted', 'envelope_fallbacks', 'envelope_position_witnesses', 'envelope_tangent_witnesses'):
            summary[field] = pair['updated'][0][field]
        summaries.append(summary)
        write(output/'measurements.csv', rows)
        write(output/'attempts.csv', attempts)
        write(output/'summary.csv', summaries)
        print(f'{index+1}/{len(CASES)} {suite}:{case} status unchanged; '
              f'total {summary["baseline_total_ms"]:.3f} -> {summary["updated_total_ms"]:.3f} ms; '
              f'boundary {summary["baseline_boundary_ms"]:.3f} -> {summary["updated_boundary_ms"]:.3f} ms', flush=True)
    manifest['evidence_sha256'] = {f.name: sha(f) for f in output.glob('*.csv')}
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
