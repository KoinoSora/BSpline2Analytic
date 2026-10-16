"""Diagnose only the 240 historical main-suite failures; CURRENT validator only.

Run after the uninstrumented three-run benchmark. No diagnostic timing is
reported as production performance. Oracle/frame swaps never feed candidates.
"""
import argparse
import collections
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, rows):
    if not rows:
        raise ValueError(f'No diagnostic rows: {path}')
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'sdk', 'data', 'failures', 'baseline', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--cases', nargs='+', type=int)
    p.add_argument('--timeout', type=float, default=120)
    p.add_argument('--resume', action='store_true')
    args = p.parse_args()
    ids = [int(r['id']) for r in read(args.failures) if r['suite'] == 'main']
    if len(ids) != 240 or len(set(ids)) != 240:
        raise ValueError('Expected exactly 240 unique main failures')
    if args.cases:
        if len(set(args.cases)) != len(args.cases) or not set(args.cases) <= set(ids):
            p.error('Cases must be unique members of the 240 failures')
        ids = args.cases
    exe, sdk, data, output = (getattr(args, n).resolve() for n in ('exe', 'sdk', 'data', 'output'))
    config = dict(exe_sha256=sha(exe), failures_sha256=sha(args.failures), cases=ids,
                  inputs_sha256={str(i): sha(data/f'general-{i:03}.b2a') for i in ids},
                  baseline_results_sha256=sha(args.baseline/'results.csv'),
                  baseline_attempts_sha256=sha(args.baseline/'attempts.csv'), timeout=args.timeout,
                  production_timing=False, parallel_processes=1)
    output.mkdir(parents=True, exist_ok=True)
    stamp = output/'run.json'
    if args.resume:
        if not stamp.exists() or json.loads(stamp.read_text()) != config:
            raise ValueError('Resume needs identical binary, inputs, baseline and configuration')
    elif stamp.exists():
        raise ValueError('Output already used')
    stamp.write_text(json.dumps(config, indent=2), encoding='utf-8')
    parts = output/'cases'
    parts.mkdir(exist_ok=True)
    env = dict(os.environ)
    key = next(k for k in env if k.lower() == 'path')
    env[key] = str(sdk/'x64-win/bin')+os.pathsep+env[key]
    executions = []
    tables = {name: [] for name in ('candidates', 'charts', 'witnesses')}
    for index, case in enumerate(ids):
        prefix = parts/f'main-{case:03}'
        done, log = prefix.with_suffix('.json'), prefix.with_suffix('.log')
        if args.resume and done.exists():
            execution = json.loads(done.read_text())
        else:
            start = time.perf_counter()
            with log.open('wb') as capture:
                try:
                    result = subprocess.run([str(exe), '--diagnose-main', str(data), str(prefix), str(case), '--rational-only'],
                                            cwd=sdk, env=env, stdout=capture, stderr=subprocess.STDOUT,
                                            timeout=args.timeout, creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
                    code, timeout = result.returncode, False
                except subprocess.TimeoutExpired:
                    code, timeout = None, True
            execution = dict(id=case, exit_code=code, timeout=timeout, wall_seconds=time.perf_counter()-start,
                             log_tail=log.read_text(encoding='utf-8', errors='replace')[-1200:])
            done.write_text(json.dumps(execution, indent=2), encoding='utf-8')
        executions.append(execution)
        if execution['timeout'] or execution['exit_code'] or f'COMPLETE {case} ' not in execution['log_tail']:
            (output/'executions.json').write_text(json.dumps(executions, indent=2))
            raise RuntimeError(f'Incomplete diagnostic case {case}: {execution}')
        for name in tables:
            rows = read(Path(str(prefix)+f'-{name}.csv'))
            if any(r['id'] != str(case) for r in rows):
                raise ValueError('Diagnostic identity mismatch')
            tables[name].extend(rows)
        if index % 10 == 0 or index + 1 == len(ids):
            print(f'{index+1}/{len(ids)} main:{case}: {execution["wall_seconds"]:.2f}s', flush=True)
    baseline = {(r['id'], r['attempt']): r for r in read(args.baseline/'attempts.csv') if r['run'] == '0' and int(r['id']) in ids}
    current = {(r['id'], r['candidate']): r for r in tables['candidates'] if r['mode'] == 'classified'}
    if set(baseline) != set(current):
        raise ValueError('Candidate sequence changed')
    for key, row in current.items():
        for name in ('type', 'status', 'accepted', 'charts', 'nodes', 'message'):
            if row[name] != baseline[key][name]:
                raise ValueError(f'Instrumented result differs from baseline: {key} {name}')
    for name, rows in tables.items():
        write(output/(name+'.csv'), rows)
    (output/'executions.json').write_text(json.dumps(executions, indent=2), encoding='utf-8')
    print('All candidate decisions match uninstrumented baseline', flush=True)
    print('Recovery branches:', dict(collections.Counter(r['reason'] for r in tables['charts'] if r['mode'] == 'classified')), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
