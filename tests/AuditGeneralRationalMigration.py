"""Audit V2 file identity, provenance, representation and SOURCE certificates.

Only reads files. Does not invoke classification, conversion or a validator.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
from RationalDatasetSupport import require_v2_data


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def fixture(path):
    lines = path.read_text().splitlines()
    v2 = lines[0] == 'B2A_RATIONAL_DOMAIN_V2'
    if lines[0] not in ('B2A_RATIONAL_DOMAIN_V2', 'B2A_GENERAL_V1'):
        raise ValueError('Invalid fixture')
    meta = list(map(int, lines[1].split()))
    p, q, nu, nv, rational = map(int, lines[11 + v2].split())
    return dict(meta=meta, degree_u=p, degree_v=q, controls=nu*nv,
                mother=[float(x) for line in lines[2:6] for x in line.split()],
                bend=float(lines[10 + v2].split()[0]),
                tolerance=float(lines[10 + v2].split()[4]), rational=rational,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('data', 'legacy-main', 'legacy-extreme'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    root = require_v2_data(a.data)
    stamp = json.loads((root/'dataset.json').read_text())
    extreme_meta = {int(r['id']): r for r in read(a.legacy_extreme/'manifest.csv')}
    rows, stats = [], {}
    for suite, size in stamp['sizes'].items():
        audit_path = root/suite/f'audit-0-{size}.csv'
        if hashlib.sha256(audit_path.read_bytes()).hexdigest() != stamp['qualification_csv_sha256'][suite]:
            raise ValueError('Source qualification hash mismatch')
        audits = {int(r['id']): r for r in read(audit_path)}
        if set(audits) != set(range(size)):
            raise ValueError('Missing source audit rows')
        for i in range(size):
            current = fixture(root/suite/f'{suite}-{i:03}.b2a')
            legacy_id = i*100+4 if suite == 'synthetic' else i
            old_path = (a.legacy_extreme/f'extreme-{i:03}.b2a' if suite == 'extreme' else
                        a.legacy_main/f'general-{"extra-" if suite == "extra" else ""}{legacy_id:03}.b2a')
            old = fixture(old_path)
            if suite != 'synthetic' and current['meta'] != old['meta']:
                raise ValueError(f'Case identity/role/degeneracy changed: {suite}:{i}')
            if any(not math.isclose(x, y, rel_tol=1e-14, abs_tol=1e-15) for x, y in zip(current['mother'], old['mother'])):
                raise ValueError('Mother geometry changed')
            if not 1 <= current['controls'] <= 200 or current['rational'] != current['meta'][6]:
                raise ValueError('Representation/control budget violated')
            nonlinear = suite == 'extra' and i < 20
            if bool(current['bend']) != nonlinear:
                raise ValueError('Unexpected nonlinear reference boundary')
            audit = audits[i]
            if audit['source_certified'] != '1':
                raise ValueError('Uncertified source')
            if current['meta'][1] == 1:
                if not float(audit['rank_margin']) > 0:
                    raise ValueError('Missing global infeasibility proof')
                if not (0 <= float(audit['boundary_source_bound']) < current['tolerance'] and
                        0 <= float(audit['source_tangent_deg']) < .01):
                    raise ValueError('Infeasible input boundary missed the straight NURBS reference')
            elif not (0 <= float(audit['source_bound']) < current['tolerance'] and
                      0 <= float(audit['source_tangent_deg']) < .01):
                raise ValueError('Source bound does not meet construction tolerance')
            scenario = extreme_meta[i]['group'] if suite == 'extreme' else ''
            # These are inherited SCENARIO names; measured sizes live in audit CSV.
            row = dict(suite=suite, id=i, legacy_file=old_path.name, legacy_sha256=old['sha256'],
                       input_sha256=current['sha256'], scenario=scenario, type=current['meta'][2],
                       rational=current['rational'], collapsed_mask=current['meta'][8],
                       old_degree_u=old['degree_u'], old_degree_v=old['degree_v'], old_controls=old['controls'],
                       degree_u=current['degree_u'], degree_v=current['degree_v'], controls=current['controls'],
                       reference_boundary='nonlinear_probe' if nonlinear else 'straight_nurbs_domain',
                       source_bound=audit['source_bound'], source_tangent_deg=audit['source_tangent_deg'])
            rows.append(row)
        part = [r for r in rows if r['suite'] == suite]
        stats[suite] = dict(cases=size, rational=sum(r['rational'] for r in part),
                           collapsed_cases=sum(bool(r['collapsed_mask']) for r in part),
                           min_degree=min(min(r['degree_u'], r['degree_v']) for r in part),
                           max_degree=max(max(r['degree_u'], r['degree_v']) for r in part),
                           over_biquintic=sum(max(r['degree_u'], r['degree_v']) > 5 for r in part),
                           control_range=[min(r['controls'] for r in part), max(r['controls'] for r in part)],
                           degree_increased=sum(r['degree_u'] > r['old_degree_u'] or r['degree_v'] > r['old_degree_v'] for r in part),
                           controls_decreased=sum(r['controls'] < r['old_controls'] for r in part),
                           scenarios=dict(Counter(r['scenario'] for r in part if r['scenario'])))
    with (root/'migration.csv').open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary = dict(contract='rational-domain-v2', cases=len(rows), source_certified=len(rows),
                   representation_policy=stamp['representation_policy'],
                   fallback_tests_run=False, suites=stats,
                   qualification_exe_sha256=stamp['qualification_exe_sha256'],
                   migration_csv_sha256=hashlib.sha256((root/'migration.csv').read_bytes()).hexdigest())
    (root/'qualification.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
