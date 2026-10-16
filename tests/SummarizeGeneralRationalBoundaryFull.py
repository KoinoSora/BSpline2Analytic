"""Summarize full current-validator coverage, timings and every unresolved case."""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import statistics
from RunGeneralRationalBoundaryFull import read, write

TIMES = ('total_ms', 'classify_ms', 'validate_ms', 'sampling_ms', 'recovery_ms', 'prepare_ms', 'boundary_ms')
COUNTS = ('attempts', 'chart_attempts', 'nodes', 'refinements', 'fast_checks', 'fast_accepted', 'tangent_fallbacks', 'position_only', 'fused_errors',
          'tangent_residual_checks', 'tangent_residual_accepted', 'tangent_residual_interval_fallbacks')


def numbers(rows, key):
    values = []
    for row in rows:
        try:
            value = float(row.get(key, 'nan'))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
    return values


def quantile(values, q):
    return sorted(values)[max(0, math.ceil(q * len(values)) - 1)] if values else None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', type=Path)
    args = p.parse_args()
    directory = args.directory
    config = json.loads((directory / 'run.json').read_text())
    if config.get('dataset_contract') != 'rational-domain-v2':
        raise ValueError('V2 summary requires a V2 run; historical reports remain archived')
    provenance_path = directory / 'frozen-source.json'
    provenance = json.loads(provenance_path.read_text()) if provenance_path.exists() else {}
    if provenance and provenance.get('exe_sha256') != config['exe_sha256']:
        raise ValueError('Frozen-source provenance does not match measured executable')
    rows, attempts, known_rows = (read(directory / name) for name in ('results.csv', 'attempts.csv', 'known-mother.csv'))
    grouped, by_attempt = defaultdict(list), defaultdict(list)
    for row in rows:
        grouped[(row['suite'], int(row['id']))].append(row)
        timed = numbers([row], 'total_ms')
        if timed:
            if float(row['classify_ms']) + float(row['validate_ms']) > timed[0] + 1e-6:
                raise ValueError('Classification/validation timing overlap')
            if sum(float(row[k]) for k in ('sampling_ms', 'recovery_ms', 'prepare_ms', 'boundary_ms')) > float(row['validate_ms']) + 1e-6:
                raise ValueError('Validation phase timing overlap')
            if row.get('source_certified') != '1':
                raise ValueError('Measured an uncertified fixture')
            if row.get('accepted') == '1':
                for key, limit in (('boundary_bound', 1e-5), ('tangent_deg', .1), ('mother_bound', 1e-6)):
                    value = float(row[key])
                    if not math.isfinite(value) or not 0 <= value < limit:
                        raise ValueError('Accepted result violates its reported certificate bounds')
    for row in attempts:
        by_attempt[(row['suite'], int(row['id']))].append(row)
    known = {(r['suite'], int(r['id'])): r for r in known_rows}
    if set(grouped) != {tuple(x) for x in config['cases']}:
        raise ValueError('Missing/extra cases')
    cases = []
    for key, current in grouped.items():
        suite, case = key
        if sorted(int(r['run']) for r in current) != list(range(config['runs'])):
            raise ValueError(f'Incomplete/duplicate run keys: {key}')
        first = current[0]
        accepted = sum(r['accepted'] == '1' for r in current)
        audited = sum(r['accepted'] == '1' and r['audit_status'] == 'passed' for r in current)
        measured = numbers(current, 'total_ms')
        cancelled = sum(r.get('time_budget_exceeded') == '1' for r in current)
        overruns = sum(r.get('within_time_limit') == '0' for r in current)
        expectation = config["expectations"][f"{suite}:{case}"]
        if len(measured) != config['runs'] or any(r['execution_status'] != 'completed' for r in current):
            outcome = 'execution_incomplete'
        elif overruns:
            outcome = 'wall_time_exceeded'
        elif cancelled:
            outcome = 'time_budget_exceeded'
        elif accepted and accepted != config['runs']:
            outcome = 'unstable_acceptance'
        elif accepted != audited:
            outcome = 'audit_unresolved'
        elif expectation == 'accept':
            outcome = 'accepted' if accepted else 'not_certified_positive'
        elif expectation in ('reject_infeasible', 'reject_known_mother', 'reject_budget'):
            outcome = 'unexpected_acceptance' if accepted else 'expected_rejection'
        elif expectation == 'invalid_options':
            outcome = 'expected_invalid' if all(r['status'] == 'InvalidInput' for r in current) else 'invalid_option_mismatch'
        else:
            outcome = 'probe_accepted' if accepted else 'probe_rejected'
        item = {k: first.get(k, '') for k in ('suite', 'id', 'group', 'source_type', 'rational', 'degree_u', 'degree_v', 'controls', 'entry')}
        item.update(expectation=expectation, outcome=outcome, measured_runs=len(measured), accepted_runs=accepted, audited_runs=audited,
                    statuses='|'.join(sorted({r['status'] for r in current})),
                    median_total_ms=statistics.median(measured) if measured else None,
                    min_total_ms=min(measured) if measured else None, max_total_ms=max(measured) if measured else None,
                    runs_under_10ms=sum(t < 10 for t in measured))
        item.update(time_budget_runs=cancelled, over_time_limit_runs=overruns,
                    timeout_stages='|'.join(sorted({r.get('timeout_stage', '') for r in current} - {''})))
        for k in TIMES[1:] + COUNTS:
            values = numbers(current, k)
            item['median_' + k] = statistics.median(values) if values else None
        ka = known.get(key, {})
        item.update(known_status=ka.get('status', ''), known_accepted=ka.get('accepted', ''), known_message=ka.get('message', ''))
        candidates = by_attempt[key]
        def progress(a):
            return (float(a['nodes']) > 0) * 4 + (float(a['prepare_ms']) > 0) * 2 + (float(a['charts']) > 0)
        deepest = max(candidates, key=progress) if candidates else {}
        item['deepest_candidate_message'] = deepest.get('message', '')
        item['last_message'] = current[-1]['message']
        item['audit_message'] = '|'.join(sorted({r.get('audit_message', '') for r in current} - {''}))
        item['output_types'] = '|'.join(sorted({r.get('output_type', '') for r in current} - {'Unknown', ''}))
        if cancelled:
            item['diagnosis'] = 'time_budget/' + item['timeout_stages']
        elif outcome == 'not_certified_positive':
            if ka.get('accepted') == '1':
                item['diagnosis'] = 'known_mother_passes_candidate_path_fails'
            elif 'recovery unresolved' in ka.get('message', ''):
                item['diagnosis'] = 'known_mother_chart_recovery_unresolved'
            elif 'correspondence not certified' in ka.get('message', ''):
                item['diagnosis'] = 'known_mother_boundary_unresolved'
            else:
                item['diagnosis'] = ka.get('message', '') or item['deepest_candidate_message']
        else:
            item['diagnosis'] = ''
        cases.append(item)
    write(directory / 'case-summary.csv', cases)
    write(directory / 'not-accepted-cases.csv', [c for c in cases if c['audited_runs'] != config['runs']])
    failures = [c for c in cases if c['outcome'] not in ('accepted', 'expected_rejection', 'expected_invalid', 'probe_accepted', 'probe_rejected')]
    write(directory / 'failures.csv', failures)
    slow = sorted([c for c in cases if c['max_total_ms'] is not None and c['max_total_ms'] >= 10], key=lambda c: c['max_total_ms'], reverse=True)
    write(directory / 'over-10ms.csv', slow)
    write(directory / 'over-time-limit.csv', [c for c in cases if c['over_time_limit_runs']])
    write(directory / 'time-budget-cases.csv', [c for c in cases if c['time_budget_runs']])
    def aggregate(selected):
        timed = [c['median_total_ms'] for c in selected if c['median_total_ms'] is not None]
        successful = [c for c in selected if c['outcome'] == 'accepted']
        return dict(cases=len(selected), outcomes=dict(Counter(c['outcome'] for c in selected)),
                    time_budget_cases=sum(c['time_budget_runs'] > 0 for c in selected),
                    over_time_limit_cases=sum(c['over_time_limit_runs'] > 0 for c in selected),
                    measured_cases=len(timed), median_case_ms=statistics.median(timed) if timed else None,
                    p95_case_median_ms=quantile(timed, .95), p99_case_median_ms=quantile(timed, .99),
                    max_case_median_ms=max(timed) if timed else None,
                    max_run_ms=max((c['max_total_ms'] for c in selected if c['max_total_ms'] is not None), default=None),
                    median_under_10_cases=sum(c['median_total_ms'] is not None and c['median_total_ms'] < 10 for c in selected),
                    all_runs_under_10_cases=sum(c['measured_runs'] == config['runs'] and c['runs_under_10ms'] == config['runs'] for c in selected),
                    accepted_and_all_runs_under_10=sum(c['runs_under_10ms'] == config['runs'] for c in successful),
                    stages_ms={k: statistics.median([c['median_' + k] for c in selected if c['median_' + k] is not None])
                               for k in TIMES[1:] if any(c['median_' + k] is not None for c in selected)})
    groups = {}
    for suite in ('main', 'extra', 'extreme', 'synthetic'):
        subset = [c for c in cases if c['suite'] == suite]
        groups[suite] = aggregate(subset)
        if suite == 'main':
            groups['main_positive'] = aggregate([c for c in subset if c['expectation'] == 'accept'])
            groups['main_negative'] = aggregate([c for c in subset if c['expectation'] == 'reject_infeasible'])
        for name in sorted({c['group'] for c in subset}):
            groups[suite + '/' + name] = aggregate([c for c in subset if c['group'] == name])
        for name in sorted({c['source_type'] for c in subset if c['expectation'] == 'accept'}):
            groups[suite + '/positive/' + name] = aggregate([c for c in subset if c['source_type'] == name and c['expectation'] == 'accept'])
        for name in ('accepted', 'not_certified_positive'):
            groups[suite + '/' + name] = aggregate([c for c in subset if c['outcome'] == name])
    summary = dict(dataset_contract='rational-domain-v2', fixed_chart=config['fixed_chart'], algorithm_commit=provenance.get('algorithm_commit', 'unspecified'), cases=len(cases), rows=len(rows),
                   runs=config['runs'], quantiles='nearest rank over per-case medians; not production latency percentiles',
                   outcomes=dict(Counter(c['outcome'] for c in cases)), groups=groups,
                   failure_diagnoses=dict(Counter(c['diagnosis'] for c in failures)),
                   slowest_cases=sorted([c for c in cases if c['median_total_ms'] is not None], key=lambda c: c['median_total_ms'], reverse=True)[:30],
                   over_10ms_cases=len(slow),
                   max_case_ms=config.get('max_case_ms'),
                   time_budget_cases=sum(c['time_budget_runs'] > 0 for c in cases),
                   over_time_limit_cases=sum(c['over_time_limit_runs'] > 0 for c in cases),
                   input_certificate_errors=sum(r.get('source_certified') == '0' for r in rows),
                   independent_audit_unresolved_runs=sum(r.get('audit_status') == 'unresolved' for r in rows),
                   interior_blindspot=read(directory / 'cases/extra-049-interior.csv'))
    (directory / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    print(json.dumps({k: summary[k] for k in ('cases', 'rows', 'outcomes', 'over_10ms_cases', 'failure_diagnoses')}, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
