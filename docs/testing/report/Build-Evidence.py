"""Recalculate report statistics from committed records; never run geometry."""
from pathlib import Path
import collections
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
A = ROOT / 'docs/testing/results/sequential-non-extreme-20261016'
B = ROOT / 'docs/testing/general-rational-985-50ms-2026-10-16'
D = ROOT / 'docs/testing/introduction-cases-2026-10-16'


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    raw, cases = rows(A / 'runs.csv'), rows(A / 'cases.csv')
    assert len(raw) == 3480 and len(cases) == 1160
    grouped = collections.defaultdict(list)
    for r in raw:
        grouped[(r['suite'], int(r['id']))].append(r)
        assert r['source_certified'] == r['iso_attempted'] == '1'
        assert math.isfinite(float(r['total_ms'])) and float(r['total_ms']) >= float(r['api_ms'])
    output = []
    for c in cases:
        g = grouped[c['suite'], int(c['id'])]
        assert {int(r['run']) for r in g} == {0, 1, 2}
        maximum = max(float(r['total_ms']) for r in g)
        accepted = all(r['accepted'] == '1' for r in g)
        audited = accepted and all(r['audit'] == 'passed' for r in g)
        timeout = any(r['time_budget_exceeded'] == '1' or float(r['total_ms']) > 50 for r in g)
        audit_failed = any(r['audit'] == 'failed' for r in g)
        false_accept = c['category'] == 'negative' and any(r['accepted'] == '1' for r in g)
        correct = not timeout and not audit_failed and (
            audited if c['category'] == 'positive' else not false_accept if c['category'] == 'negative' else False)
        assert str(correct) == c['correct']
        assert math.isclose(maximum, float(c['max_ms']), rel_tol=1e-12)
        assert str(audited) == c['audited_success'] and str(timeout) == c['timeout']
        assert (maximum < 10) == (c['under_10ms'] == 'True')
        output.append(dict(suite=c['suite'], id=int(c['id']), category=c['category'],
                           accepted=accepted, audited=audited, timeout=timeout,
                           max_ms=maximum, under10=maximum < 10, correct=correct,
                           fallback=c['fallback'] == 'True', path=c['path']))
    def stats(selected):
        keys = {(x['suite'], x['id']) for x in selected}
        rr = [r for r in raw if (r['suite'], int(r['id'])) in keys]
        return dict(cases=len(selected), calls=len(rr), audited_success=sum(x['audited'] for x in selected),
                    under10=sum(x['under10'] for x in selected),
                    audited_under10=sum(x['audited'] and x['under10'] for x in selected),
                    ordinary_correct_under10=sum(x['correct'] and x['under10'] for x in selected),
                    budget_failures=sum(x['timeout'] for x in selected),
                    mean_ms=sum(float(r['total_ms']) for r in rr)/len(rr),
                    max_ms=max(float(r['total_ms']) for r in rr),
                    total_ms=sum(float(r['total_ms']) for r in rr))
    summary = {s: stats([x for x in output if x['suite'] == s]) for s in ('iso', 'main', 'extra')}
    summary['all'] = stats(output)
    assert summary['all']['audited_success'] == 973
    assert summary['all']['audited_under10'] == 943
    assert summary['all']['ordinary_correct_under10'] == 894
    run = json.loads((A / 'run.json').read_text())
    for r in run['cases']:
        p = ROOT / 'tests/data' / ('dataset500' if r['suite'] == 'iso' else 'general-rational-v2/' + r['suite'])
        p /= f"case-{r['id']:03}.b2a" if r['suite'] == 'iso' else f"{r['suite']}-{r['id']:03}.b2a"
        assert sha(p) == r['sha256'], p
    provenance = json.loads((A / 'provenance.json').read_text())
    differences = [dict(path=p, recorded_sha256=h, current_sha256=sha(ROOT/p))
                   for p, h in provenance['source_sha256'].items() if sha(ROOT/p) != h]
    br = rows(B/'results.csv')
    bc = rows(B/'case-summary.csv')
    assert len(br) == 2955 and len(bc) == 985
    assert sum(r['accepted'] == '1' for r in br) == 2133
    positive = [r for r in bc if r['expectation'] == 'accept']
    assert len(positive) == 850 and sum(int(r['accepted_runs']) == 3 for r in positive) == 710
    assert sum(float(r['max_total_ms']) < 10 for r in bc) == 951
    dr = rows(D/'runs.csv')
    selected_cases = []
    for name in ('sphere-317', 'cone-217', 'torus-triangle'):
        rr = [r for r in dr if r['case'] == name]
        assert len(rr) == 3 and all(r['accepted'] == r['audit'] == '1' for r in rr)
        assert all(r['iso_attempted'] == r['fallback_attempted'] == '1' for r in rr)
        selected_cases.append(dict(case=name, max_ms=max(float(r['total_ms']) for r in rr),
            boundary_bound=max(float(r['boundary_bound']) for r in rr),
            tangent_bound_deg=max(float(r['tangent_bound_deg']) for r in rr),
            audit_position_max=max(float(r['audit_position_max']) for r in rr),
            audit_edge_max=max(float(r['audit_edge_max']) for r in rr),
            output_mother_bound=max(float(r['output_mother_bound']) for r in rr)))
    dp = json.loads((D/'provenance.json').read_text())
    assert all(sha(ROOT/p) == h for p,h in dp['source_sha256'].items())
    records = [A/'runs.csv', A/'cases.csv', A/'run.json', A/'provenance.json', A/'sequential-unit.log',
               B/'results.csv', B/'case-summary.csv', B/'summary.json', B/'frozen-source.json',
               B/'extra-049-interior.csv',
               ROOT/'docs/testing/results/iso-full-20261016/manifest.json',
               ROOT/'docs/testing/ISO_SPECIALIZED_FULL_2026-10-16.md',
               D/'runs.csv', D/'summary.json', D/'provenance.json',
               ROOT/'docs/requirements/中船海舟样条转二次曲面.pdf']
    report = dict(geometry_executed=False, primary_batch=summary,
        successful_path_counts=dict(collections.Counter(x['path'] for x in output if x['audited'])),
        budget_exit_stages=dict(collections.Counter(r['timeout_stage'] for r in raw if r['time_budget_exceeded']=='1')),
        actual_over50_calls=sum(float(r['total_ms'])>50 for r in raw),
        accepted_calls_with_optimization=sum(r['accepted']=='1' and int(r['optimization_iterations'])>0 for r in raw),
        source_differences_from_primary_batch=differences, primary_input_hashes_verified=1160,
        general_batch=dict(cases=985,calls=2955,positive=850,positive_passed=710,under10=951,
                           max_ms=max(float(r['total_ms']) for r in br)),
        selected_current_source_hashes_match=True, selected_cases=selected_cases,
        evidence_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in records})
    (OUT/'report-evidence.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    failures=[x for x in output if not x['audited']]
    assert len(failures)==187
    with (OUT/'primary-not-accepted.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(failures[0]));writer.writeheader();writer.writerows(failures)
    print(json.dumps(dict(primary=summary['all'],source_differences=len(differences),selected=selected_cases),indent=2))


if __name__ == '__main__':
    main()
