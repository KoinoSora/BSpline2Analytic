"""Compare local Release binaries on 12 fixed fixtures, without a full dataset run.

The first --variant is the baseline. Each executable retains its independent
input/output certificates. Timing excludes these audits, as in RunGeneralHybrid.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess

CASES = [(0, i) for i in (0, 100, 101, 218, 254, 318, 400, 401, 500)] + [
    (1, i) for i in (0, 50, 55)]
INVARIANTS = ('status', 'accepted', 'output_type', 'output_audit', 'attempts',
              'samples', 'boundary_nodes', 'boundary_segments', 'mapping_refinements',
              'collapsed_edges', 'chart_attempts', 'chart_iterations')
BOUNDS = ('interior_max', 'finite_sample_max', 'boundary_bound', 'tangent_deg')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant', action='append', required=True, help='name=executable')
    parser.add_argument('--sdk', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--runs', type=int, default=5)
    args = parser.parse_args()
    if args.runs < 1:
        parser.error('--runs must be positive')
    variants = [(name, Path(path).resolve()) for name, path in
                (v.split('=', 1) for v in args.variant)]
    if len(set(name for name, _ in variants)) != len(variants):
        parser.error('variant names must be unique')
    sdk, data, output = args.sdk.resolve(), args.data.resolve(), args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    work = output.parent / (output.stem + '-cases')
    work.mkdir(exist_ok=True)
    env = dict(os.environ)
    key = next(k for k in env if k.lower() == 'path')
    env[key] = str(sdk / 'x64-win/bin') + os.pathsep + env[key]
    rows, inputs = [], {}
    for index, (extra, case_id) in enumerate(CASES):
        filename = f'general-{case_id:03}.b2a' if not extra else f'general-extra-{case_id:03}.b2a'
        inputs[filename] = sha(data / filename)
        # Rotate variant order between cases to reduce systematic order bias.
        offset = index % len(variants)
        for name, exe in variants[offset:] + variants[:offset]:
            target = work / f'{index:02}-{variants.index((name, exe))}.csv'
            if target.exists():
                target.unlink()
            run = subprocess.run([str(exe), '--case', str(data), str(target), str(case_id),
                                  str(extra), str(args.runs)], cwd=sdk, env=env,
                                 capture_output=True, text=True, encoding='utf-8',
                                 errors='replace', timeout=30)
            target.with_suffix('.log').write_text(run.stdout + run.stderr, encoding='utf-8')
            if run.returncode:
                raise RuntimeError(f'{name} {filename}: exit {run.returncode}; see {target.with_suffix(".log")}')
            current = list(csv.DictReader(target.open(encoding='utf-8')))
            if len(current) != args.runs:
                raise RuntimeError(f'{name} {filename}: missing measurements')
            for row in current:
                if row['fixture_certified'] != '1' or (row['accepted'] == '1' and row['output_audit'] != '1'):
                    raise RuntimeError(f'{name} {filename}: certificate failure')
                rows.append(dict(variant=name, **row))
        print(f'Completed {index+1}/{len(CASES)} selected cases', flush=True)
    with output.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), quoting=csv.QUOTE_ALL)
        writer.writeheader()
        writer.writerows(rows)
    baseline = variants[0][0]
    reference = {(r['extra'], r['id'], r['run']): r for r in rows if r['variant'] == baseline}
    differences, bound_changes = [], []
    for row in rows:
        ref = reference[row['extra'], row['id'], row['run']]
        for field, target in [(f, differences) for f in INVARIANTS] + [(f, bound_changes) for f in BOUNDS]:
            if row[field] != ref[field]:
                target.append(dict(variant=row['variant'], extra=row['extra'], id=row['id'],
                                   run=row['run'], field=field, before=ref[field], after=row[field]))
    summary = dict(cases=len(CASES), runs=args.runs, rows=len(rows),
                   invariant_differences=differences, bound_changes=bound_changes,
                   inputs_sha256=inputs, exe_sha256={n: sha(p) for n, p in variants},
                   runner_sha256=sha(Path(__file__).resolve()), cpu_affinity_fixed=False,
                   exclusive_machine=False, variants={})
    for name, _ in variants:
        positive = [r for r in rows if r['variant'] == name and r['group'] == 'feasible']
        ids = sorted({r['id'] for r in positive}, key=int)
        per_case = {i: {k: statistics.median(float(r[k]) for r in positive if r['id'] == i)
                       for k in ('prepare_ms', 'boundary_ms', 'validate_ms', 'pipeline_ms')}
                    for i in ids}
        summary['variants'][name] = dict(positive_case_medians=per_case,
            mean_of_case_medians={k: statistics.mean(v[k] for v in per_case.values())
                                  for k in ('prepare_ms', 'boundary_ms', 'validate_ms', 'pipeline_ms')})
    output.with_suffix('.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k: v['mean_of_case_medians'] for k, v in summary['variants'].items()}, indent=2))
    print(f'Invariant differences: {len(differences)}; bound changes: {len(bound_changes)}')
    return int(bool(differences))


if __name__ == '__main__':
    raise SystemExit(main())
