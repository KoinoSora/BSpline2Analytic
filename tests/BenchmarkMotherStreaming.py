"""Run only selected mother-distance fixtures, then audit extraction envelopes."""
import argparse
import csv
from decimal import localcontext
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
from AuditMotherExtraction import audit

CASES = [1, 101, 201, 301, 351, 381, 383, 391, 395, 401, 451, 481, 485, 489, 491, 497]
METRICS = ['classification_ms', 'validation_ms', 'algorithm_ms', 'preparation_ms',
           'position_ms', 'subdivision_ms', 'boundary_ms', 'reporting_ms',
           'source_extraction_ms', 'source_boundary_ms', 'residual_build_ms']

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    for name in ['baseline', 'streaming', 'sdk', 'dataset', 'work', 'output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--runs', type=int, default=5)
    args = parser.parse_args()
    for name in ['baseline', 'streaming', 'sdk', 'dataset', 'work', 'output']:
        setattr(args, name, getattr(args, name).resolve())
    args.work.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['PATH'] = str(args.sdk/'x64-win/bin')+os.pathsep+env['PATH']
    rows = {'baseline': [], 'streaming': []}
    audit_log, logs, count, ratio = [], [], 0, 0
    manifest = {'baseline_commit': 'ed4d63b9de7d8cb0f705ff04eefa0a40bb68e6f5',
                'runs': args.runs, 'cases': CASES,
                'executable_sha256': {k: sha(getattr(args, k)) for k in rows},
                'fixture_sha256': {str(i): sha(args.dataset/f'case-{i:03}.b2a') for i in CASES},
                'commands': []}
    for index, case in enumerate(CASES):
        order = ['baseline', 'streaming'] if index % 2 == 0 else ['streaming', 'baseline']
        for mode in order:
            dest = (args.work/f'{mode}-{case:03}.csv').resolve()
            command = [str(getattr(args, mode).resolve()), str(args.dataset.resolve()),
                       str(case), str(dest), str(args.runs)]
            manifest['commands'].append(command)
            run = subprocess.run(command, cwd=args.sdk, env=env, capture_output=True,
                                 text=True, errors='replace', timeout=120)
            logs.append(f'{mode} {case}\n{run.stdout}\n{run.stderr}')
            if run.returncode:
                raise RuntimeError(logs[-1])
            with dest.open(newline='') as stream:
                current = list(csv.DictReader(stream))
            if len(current) != 2*args.runs or any(r['status'] != 'Passed' for r in current):
                raise AssertionError(f'{mode} {case}: unexpected result count or status')
            rows[mode].extend(current)
            if mode == 'streaming':
                with localcontext() as ctx:
                    ctx.prec = 80
                    n, r = audit(Path(str(dest)+'.extraction.json'))
                count += n
                ratio = max(ratio, r)
                audit_log.append(f'case {case}: coordinates={n} max_error/envelope={r:.9g} PASSED')
            print(f'{mode} {case} PASSED', flush=True)
    for mode, records in rows.items():
        with (args.output/f'{mode}.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, records[0].keys())
            writer.writeheader()
            writer.writerows(records)
    summary = []
    for mode, records in rows.items():
        for case in CASES:
            for cache in ['cold', 'warm']:
                group = [r for r in records if int(r['id']) == case and r['cache'] == cache]
                row = dict(backend=mode, id=case, type=group[0]['type'], cache=cache, runs=len(group))
                for metric in METRICS:
                    values = [float(r[metric]) for r in group if r[metric]]
                    row[metric] = statistics.mean(values) if values else ''
                summary.append(row)
    with (args.output/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)
    audit_log.append(f'TOTAL coordinates={count} max_error/envelope={ratio:.9g} PASSED')
    (args.output/'extraction-audit.txt').write_text('\n'.join(audit_log)+'\n')
    (args.output/'runs.log').write_text('\n'.join(logs).rstrip()+'\n', encoding='utf-8')
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(audit_log[-1])

if __name__ == '__main__':
    main()
