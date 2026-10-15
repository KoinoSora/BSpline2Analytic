"""Compare two revisions of the CURRENT NURBS boundary validator only.

The rational-only command rejects retired hybrid-comparison executables. Five
rounds alternate process order; source certificates and SDK output audits remain
outside algorithm timers. This is a selected regression, not a dataset sweep.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess

CASES = [(0, i) for i in (0, 100, 200, 218, 318, 254, 500)] + [(1, i) for i in range(5)]
TIMES = ('pipeline_ms', 'validate_ms', 'classify_ms', 'sampling_ms', 'recovery_ms', 'prepare_ms', 'boundary_ms')
COUNTS = ('nodes', 'refinements', 'max_degree', 'fused_errors', 'position_only', 'fast_checks',
          'fast_accepted', 'tangent_fallbacks', 'source_degree_u', 'source_degree_v', 'source_controls')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'optimized', 'sdk', 'data', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--runs', type=int, default=5)
    args = parser.parse_args()
    if args.runs < 1:
        parser.error('--runs must be positive')
    exe = {name: getattr(args, name).resolve() for name in ('baseline', 'optimized')}
    sdk, data, output = (getattr(args, name).resolve() for name in ('sdk', 'data', 'output'))
    work = output.parent / (output.stem + '-cases')
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    key = next(k for k in env if k.lower() == 'path')
    env[key] = str(sdk / 'x64-win/bin') + os.pathsep + env[key]
    rows, failures = [], []
    for synthetic, case in CASES:
        for round_id in range(args.runs):
            order = ('baseline', 'optimized') if round_id % 2 == 0 else ('optimized', 'baseline')
            for variant in order:
                target = work / f'{synthetic}-{case:03}-{round_id}-{variant}.csv'
                if target.exists():
                    target.unlink()
                try:
                    result = subprocess.run([str(exe[variant]), '--rational-case', str(data), str(target),
                                             str(case), '1', str(synthetic), '--rational-only'],
                                            cwd=sdk, env=env, capture_output=True, text=True,
                                            encoding='utf-8', errors='replace', timeout=60)
                    target.with_suffix('.log').write_text(result.stdout + result.stderr, encoding='utf-8')
                    if result.returncode:
                        failures.append(dict(case=case, synthetic=synthetic, round=round_id,
                                             variant=variant, exit_code=result.returncode))
                except subprocess.TimeoutExpired:
                    failures.append(dict(case=case, synthetic=synthetic, round=round_id, variant=variant, timeout=True))
                with target.open(encoding='utf-8') if target.exists() else open(os.devnull) as stream:
                    current = list(csv.DictReader(stream))
                if len(current) != 1 or any(r['mode'] != 'rational' for r in current):
                    failures.append(dict(case=case, synthetic=synthetic, round=round_id, variant=variant, invalid_rows=True))
                for row in current:
                    row.update(variant=variant, round=str(round_id))
                rows.extend(current)
        print(f'Completed synthetic={synthetic}, id={case}', flush=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with output.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary = dict(validator='rational_boundary', runs=args.runs, rows=len(rows), failures=failures,
                   executables_sha256={name: sha(path) for name, path in exe.items()},
                   inputs_sha256={f'general-{i:03}.b2a': sha(data / f'general-{i:03}.b2a') for s, i in CASES if not s},
                   cpu_affinity_fixed=False, exclusive_machine=False, cases=[])
    for synthetic, case in CASES:
        entry = dict(synthetic=synthetic, id=case)
        for variant in exe:
            selected = [r for r in rows if r['synthetic'] == str(synthetic) and r['id'] == str(case) and r['variant'] == variant]
            entry[variant] = dict(accepted_runs=sum(r['accepted'] == '1' for r in selected),
                                  outcomes=sorted(set(r['status'] for r in selected)),
                                  medians={k: statistics.median(float(r[k]) for r in selected) for k in TIMES + COUNTS
                                           if selected and all(k in r for r in selected)},
                                  max_pipeline_ms=max((float(r['pipeline_ms']) for r in selected), default=None))
        if entry['baseline']['outcomes'] != entry['optimized']['outcomes'] or entry['baseline']['accepted_runs'] != entry['optimized']['accepted_runs']:
            failures.append(dict(case=case, synthetic=synthetic, changed_outcome=True))
        expected = 0 if not synthetic and case in (254, 500) else args.runs
        if entry['optimized']['accepted_runs'] != expected:
            failures.append(dict(case=case, synthetic=synthetic, unexpected_acceptance=True))
        summary['cases'].append(entry)
    output.with_suffix('.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(f'Finished {len(CASES)} cases, {len(rows)} rows, {len(failures)} failures', flush=True)
    return int(bool(failures))


if __name__ == '__main__':
    raise SystemExit(main())
