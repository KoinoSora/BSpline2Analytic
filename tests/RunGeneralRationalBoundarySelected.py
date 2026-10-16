"""Rational-only validation: prior slow fixtures and five 200-control fixtures.

ValidateGeneralHybrid is retired and must not run. Synthetic cases benchmark a
supplied mother and fixed chart; rebuilt V2 fixtures use classified candidates.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess

SELECTED = (0, 100, 188, 200, 218, 230, 231, 254, 255, 284, 286, 288, 318, 388, 486, 488, 500)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', type=Path, required=True)
    p.add_argument('--sdk', type=Path, required=True)
    p.add_argument('--data', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--runs', type=int, default=5)
    a = p.parse_args()
    from RationalDatasetSupport import require_v2_data
    require_v2_data(a.data)
    if a.runs < 1:
        p.error('--runs must be positive')
    exe, sdk, data, output = [x.resolve() for x in (a.exe, a.sdk, a.data, a.output)]
    work = output.parent / (output.stem + '-cases')
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    key = next(k for k in env if k.lower() == 'path')
    env[key] = str(sdk / 'x64-win/bin') + os.pathsep + env[key]
    rows, failures, hashes = [], [], {}
    for synthetic, case in [(0, i) for i in SELECTED] + [(1, i) for i in range(5)]:
        suite = 'synthetic' if synthetic else 'main'
        path = data / suite / f'{suite}-{case:03}.b2a'
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        target = work / f'{synthetic}-{case:03}.csv'
        if target.exists():
            target.unlink()
        try:
            r = subprocess.run([str(exe), '--rational-case', str(data), str(target), str(case), str(a.runs), str(synthetic), '--rational-v2-only'],
                               cwd=sdk, env=env, capture_output=True, text=True, encoding='utf-8',
                               errors='replace', timeout=60)
            target.with_suffix('.log').write_text(r.stdout + r.stderr, encoding='utf-8')
            if r.returncode:
                failures.append(dict(synthetic=synthetic, case=case, exit_code=r.returncode))
        except subprocess.TimeoutExpired:
            failures.append(dict(synthetic=synthetic, case=case, timeout=True))
        current = list(csv.DictReader(target.open(encoding='utf-8'))) if target.exists() else []
        rows.extend(current)
        if len(current) != a.runs:
            failures.append(dict(synthetic=synthetic, case=case, incomplete_rows=len(current)))
        if any(r['mode'] != 'rational' for r in current):
            failures.append(dict(synthetic=synthetic, case=case, unexpected_validator=True))
        print(f'Completed synthetic={synthetic} id={case}: rows={len(current)}', flush=True)
    with output.open('w', newline='', encoding='utf-8') as f:
        if rows:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), quoting=csv.QUOTE_ALL)
            writer.writeheader(); writer.writerows(rows)
    summary = dict(validator='rational_boundary', dataset_contract='rational-domain-v2', failures=failures, runs=a.runs, rows=len(rows),
                   exe_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(), inputs_sha256=hashes,
                   cpu_affinity_fixed=False, exclusive_machine=False, cases=[])
    for synthetic, case in [(0, i) for i in SELECTED] + [(1, i) for i in range(5)]:
        entry = dict(synthetic=synthetic, id=case)
        for mode in ('rational',):
            selected = [r for r in rows if r['synthetic'] == str(synthetic) and r['id'] == str(case) and r['mode'] == mode]
            entry[mode] = dict(accepted_runs=sum(r['accepted'] == '1' for r in selected),
                outcomes=sorted(set(r['status'] for r in selected)),
                medians={k: statistics.median(float(r[k]) for r in selected) if selected else None
                         for k in ('pipeline_ms', 'validate_ms', 'classify_ms', 'sampling_ms', 'recovery_ms',
                                   'prepare_ms', 'boundary_ms', 'attempts', 'chart_attempts', 'nodes', 'refinements', 'max_degree')})
        summary['cases'].append(entry)
    output.with_suffix('.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(f'Finished {len(summary["cases"])} selected cases, {len(rows)} rows, {len(failures)} failures', flush=True)
    return int(bool(failures))


if __name__ == '__main__':
    raise SystemExit(main())
