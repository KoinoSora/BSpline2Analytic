"""Check the 50 fixed polynomial sphere fixtures and summarize exclusive validation timing."""
import argparse
import collections
import csv
import json
import math
from pathlib import Path

STEPS = ('prepare', 'surface_bounds', 'subdivide', 'boundary', 'report')
PHASES = ('initial', 'refined')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def summarize(rows, runs):
    total = math.fsum(r['validation_total_ms'] for r in rows)
    result = {'cases_per_run': len(rows) // runs, 'measurements': len(rows),
              'validation_calls': sum(r['validation_calls'] for r in rows),
              'validation_ms_per_run': total / runs,
              'validation_mean_ms_per_case': total / len(rows), 'steps': {}}
    for step in STEPS:
        value = math.fsum(r[step + '_ms'] for r in rows)
        result['steps'][step] = {'total_ms': value, 'ms_per_run': value / runs,
                                 'mean_ms_per_case': value / len(rows),
                                 'percent': 100 * value / total if total else 0}
    for phase in PHASES:
        result[phase + '_validation_ms_per_run'] = math.fsum(
            r[phase + '_validation_total_ms'] for r in rows) / runs
    result['outer_validation_ms_per_run'] = math.fsum(
        r['classify_validate_ms'] + r['optimize_validate_ms'] for r in rows) / runs
    result['classification_compute_ms_per_run'] = math.fsum(r['classify_compute_ms'] for r in rows) / runs
    result['optimization_compute_ms_per_run'] = math.fsum(r['optimize_compute_ms'] for r in rows) / runs
    result['optimization_iterations'] = sum(r['opt_iterations'] for r in rows)
    result['surface_nodes'] = sum(r['initial_surface_bounds_calls'] + r['refined_surface_bounds_calls']
                                  for r in rows)
    result['validation_pieces'] = sum(r['validation_pieces'] for r in rows)
    result['input_surface_cells'] = sum(int(r['spans_u']) * int(r['spans_v']) for r in rows)
    result['input_boundary_cells'] = sum(2 * (int(r['spans_u']) + int(r['spans_v'])) for r in rows)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv', nargs='+', type=Path)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    baseline = {int(r['id']): r for r in read(args.baseline)}
    all_rows, per_run, mismatches = [], [], []
    for index, path in enumerate(args.csv, 1):
        rows = read(path)
        require(sorted(int(r['id']) for r in rows) == list(range(301, 400, 2)),
                'Expected exactly the 50 polynomial sphere IDs')
        for r in rows:
            require(r['type'] == 'Sphere' and r['path'] == 'polynomial', 'Wrong fixture selection')
            require(r['validation_timing_schema'] == 'exclusive-validation-v1', 'Wrong timing schema')
            for key in list(r):
                if key.endswith('_ms'):
                    r[key] = float(r[key])
                    require(math.isfinite(r[key]) and r[key] >= -1e-8, f'Invalid timing: {key}')
                elif key.endswith('_calls') or key in ('id', 'collapsed_mask', 'opt_iterations',
                                                       'validation_pieces', 'controls'):
                    r[key] = int(r[key])
            r['run'] = index
            for phase, outer in (('initial', 'classify_validate'), ('refined', 'optimize_validate')):
                total = r[phase + '_validation_total_ms']
                require(math.isclose(total, math.fsum(r[f'{phase}_{s}_ms'] for s in STEPS),
                                     rel_tol=1e-8, abs_tol=1e-8), 'Step partition mismatch')
                require(total <= r[outer + '_ms'] + 1e-8, 'Inner timing exceeds outer timing')
                require(r[phase + '_prepare_calls'] == r[outer + '_calls'], 'Lost validation calls')
                for step in STEPS:
                    require(r[f'{phase}_{step}_calls'] or r[f'{phase}_{step}_ms'] == 0,
                            'Unused stage has elapsed time')
            r['validation_total_ms'] = sum(r[p + '_validation_total_ms'] for p in PHASES)
            r['validation_calls'] = sum(r[p + '_prepare_calls'] for p in PHASES)
            for step in STEPS:
                r[step + '_ms'] = sum(r[f'{p}_{step}_ms'] for p in PHASES)
            old = baseline[r['id']]
            # Timing must not change statuses, certificates, work count or iteration count.
            for key in ('status', 'geometry_bound', 'tangent_deg', 'normal_deg',
                        'validation_pieces', 'opt_iterations', 'classify_validate_calls',
                        'optimize_validate_calls'):
                equal = r[key] == old[key] if key == 'status' else float(r[key]) == float(old[key])
                if not equal:
                    mismatches.append({'run': index, 'id': r['id'], 'field': key,
                                       'baseline': old[key], 'measured': r[key]})
        require(sum(not r['collapsed_mask'] for r in rows) == 40, 'Wrong regular case count')
        per_run.append({'run': index, 'file': path.name, 'all': summarize(rows, 1),
                        'regular': summarize([r for r in rows if not r['collapsed_mask']], 1)})
        all_rows.extend(rows)
    runs = len(args.csv)
    regular = [r for r in all_rows if not r['collapsed_mask']]
    collapsed = [r for r in all_rows if r['collapsed_mask']]
    unexpected = sorted({r['id'] for r in all_rows
                         if r['status'] != ('Degenerate' if r['collapsed_mask'] else 'Success')})
    cases = []
    for case_id in range(301, 400, 2):
        rows = [r for r in all_rows if r['id'] == case_id]
        item = {'id': case_id, 'status': rows[0]['status'], 'controls': rows[0]['controls'],
                'collapsed': bool(rows[0]['collapsed_mask'])}
        item.update({step + '_mean_ms': math.fsum(r[step + '_ms'] for r in rows) / runs
                     for step in STEPS})
        item['validation_mean_ms'] = math.fsum(r['validation_total_ms'] for r in rows) / runs
        cases.append(item)
    summary = {'runs': runs, 'cases': 50, 'baseline': args.baseline.name,
               'statuses_per_run': dict(collections.Counter(r['status'] for r in all_rows[:50])),
               'unexpected_ids': unexpected, 'baseline_mismatches': mismatches,
               'all': summarize(all_rows, runs), 'regular': summarize(regular, runs),
               'collapsed': summarize(collapsed, runs), 'per_run': per_run,
               'regular_buckets': {str(b): summarize([r for r in regular if (r['id'] % 100) // 20 == b], runs)
                                   for b in range(5)},
               'slowest_regular': max((r for r in cases if not r['collapsed']),
                                      key=lambda r: r['validation_mean_ms'])}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    with (args.output / 'case-means.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)
    print(json.dumps({k: summary[k] for k in ('runs', 'cases', 'statuses_per_run', 'unexpected_ids',
                                             'baseline_mismatches', 'regular', 'collapsed')}, indent=2))
    return bool(unexpected or mismatches)


if __name__ == '__main__':
    raise SystemExit(main())
