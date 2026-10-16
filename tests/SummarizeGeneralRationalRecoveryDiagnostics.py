"""Validate and summarize the current-validator V2 diagnostic evidence."""
import argparse
import collections
import csv
import json
import math
from pathlib import Path
import statistics


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--diagnostics', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    args = p.parse_args()
    d = args.diagnostics
    candidates, charts, witnesses = [read(d/(name+'.csv')) for name in ('candidates', 'charts', 'witnesses')]
    config = json.loads((d/'run.json').read_text())
    baseline_config = json.loads((args.baseline/'run.json').read_text())
    if config.get('dataset_contract') != 'rational-domain-v2' or baseline_config.get('dataset_contract') != 'rational-domain-v2':
        raise ValueError('Diagnostic summary requires V2 evidence')
    ids = sorted({int(r['id']) for r in candidates})
    baseline = [r for r in read(args.baseline/'results.csv') if r['suite'] == 'main' and int(r['id']) in ids]
    assert ids == sorted(config['cases']) and len(baseline) == len(ids)*baseline_config['runs']
    assert all(r['accepted'] == '0' and r['execution_status'] == 'completed' for r in baseline)
    for r in witnesses:
        if r['outside_certified'] == '1':
            assert float(r['global_lower']) > float(r['proof_tolerance']) >= 1.1e-5
            assert float(r['global_lower']) <= float(r['nearby_upper'])
        assert math.isfinite(float(r['source_box_width']))
    for r in candidates:
        if r['accepted'] == '1':
            assert float(r['boundary_bound']) < 1e-5 and float(r['tangent_deg']) < .1
    proof_keys = {(r['id'], r['mode'], r['candidate'], r['chart']) for r in witnesses if r['outside_certified'] == '1'}
    unresolved = [r for r in charts if r['reason'] == 'RationalLineWitnessExceeded' and
                  (r['id'], r['mode'], r['candidate'], r['chart']) not in proof_keys]
    with (d/'unresolved-chart-exclusions.csv').open('w', encoding='utf-8', newline='') as f:
        with (d/'charts.csv').open(encoding='utf-8-sig', newline='') as source:
            fields = csv.DictReader(source).fieldnames
        w = csv.DictWriter(f, fieldnames=fields);w.writeheader();w.writerows(unresolved)
    cases = []
    for case in ids:
        key = str(case)
        cr = [r for r in candidates if r['id'] == key]
        br = [r for r in baseline if r['id'] == key]
        assert sorted(int(r['run']) for r in br) == list(range(baseline_config['runs']))
        known = next(r for r in cr if r['mode'] == 'known_mother')
        source_type = known['type']
        classified = [r for r in cr if r['mode'] == 'classified']
        matching = [r for r in classified if r['type'] == source_type and int(r['charts'])]
        primary = [r for r in charts if r['id'] == key and r['mode'] == 'classified' and r['type'] == source_type]
        all_charts = [r for r in charts if r['id'] == key and r['mode'] == 'classified']
        excluded = sum((r['id'], r['mode'], r['candidate'], r['chart']) in proof_keys for r in primary)
        swapped = {r['mode']: r['status'] for r in cr if 'frame' in r['mode']}
        proofs = [r for r in witnesses if r['id'] == key and r['mode'] == 'classified' and r['outside_certified'] == '1']
        row = dict(id=case, type=source_type, rational=br[0]['rational'], degree_u=br[0]['degree_u'], degree_v=br[0]['degree_v'],
                   median_ms=statistics.median(float(r['total_ms']) for r in br), max_ms=max(float(r['total_ms']) for r in br),
                   candidate_attempts=len(classified), matching_mothers_reaching_chart=len(matching),
                   matching_chart_attempts=len(primary), matching_charts_with_distance_exclusion=excluded,
                   all_chart_attempts=len(all_charts), all_charts_with_distance_exclusion=sum(
                       (r['id'], r['mode'], r['candidate'], r['chart']) in proof_keys for r in all_charts),
                   classified_reasons=';'.join(sorted(set(r['reason'] for r in primary))),
                   known_status=known['status'], candidate_with_truth_frame=swapped.get('candidate_with_truth_frame', ''),
                   truth_with_candidate_frame=swapped.get('truth_with_candidate_frame', ''),
                   min_exclusion_lower=min((float(r['global_lower']) for r in proofs), default=math.nan),
                   min_nearby_upper=min((float(r['nearby_upper']) for r in proofs), default=math.nan))
        cases.append(row)
    with (d/'case-diagnosis.csv').open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(cases[0]));w.writeheader();w.writerows(cases)
    groups = {}
    for kind in sorted({r['type'] for r in cases}):
        rows = [r for r in cases if r['type'] == kind]
        times = sorted(r['median_ms'] for r in rows)
        groups[kind] = dict(cases=len(rows), median_ms=statistics.median(times), p95_case_median_ms=times[math.ceil(.95*len(times))-1],
                            max_run_ms=max(r['max_ms'] for r in rows), matching_charts=sum(r['matching_chart_attempts'] for r in rows),
                            excluded_matching_charts=sum(r['matching_charts_with_distance_exclusion'] for r in rows),
                            known_pass_ids=[r['id'] for r in rows if r['known_status'] == 'Passed'])
    summary = dict(dataset_contract='rational-domain-v2', known_mother_uses_reference_chart=True, cases=len(cases), benchmark_runs=len(baseline), diagnostic_timing_excluded=True,
                   classification_attempts=sum(r['mode'] == 'classified' for r in candidates), groups=groups,
                   classified_recovery_reasons=dict(collections.Counter(r['reason'] for r in charts if r['mode'] == 'classified')),
                   chart_rows=len(charts), witness_rows=len(witnesses),
                   excluded_chart_keys=len(proof_keys),
                   exclusions_by_mode=dict(collections.Counter(k[1] for k in proof_keys)),
                   unresolved_exclusions_by_mode=dict(collections.Counter(r['mode'] for r in unresolved)),
                   unresolved_classified_exclusion_case_ids=sorted({int(r['id']) for r in unresolved if r['mode']=='classified'}),
                   known_unresolved_exclusion_case_ids=sorted({int(r['id']) for r in unresolved if r['mode']=='known_mother'}),
                   cases_all_matching_charts_excluded=sum(r['matching_chart_attempts']>0 and r['matching_chart_attempts']==r['matching_charts_with_distance_exclusion'] for r in cases),
                   candidate_truth_frame_pass_ids=[r['id'] for r in cases if r['candidate_with_truth_frame'] == 'Passed'],
                   truth_candidate_frame_pass_ids=[r['id'] for r in cases if r['truth_with_candidate_frame'] == 'Passed'],
                   interpretation='Exclusion covers only each tested fixed chart and recovered endpoints, including export allowance; it is not a proof that no other chart or boundary is feasible.')
    (d/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
