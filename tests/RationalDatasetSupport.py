"""Shared version gate: retired analytic-UV inputs must never enter current runs."""
import json
import hashlib
from pathlib import Path

REPRESENTATION_POLICY = dict(degree_min=3, degree_max=5, max_controls=200,
                             geometry_tolerance=1e-5)


def require_representation(content):
    lines = content.decode('utf-8').splitlines()
    p, q, nu, nv, _ = map(int, lines[12].split())
    if not (3 <= p <= 5 and 3 <= q <= 5 and 0 < nu*nv <= 200):
        raise ValueError('V2 representation must use degrees 3-5 and at most 200 controls')
    meta = list(map(int, lines[1].split()))
    tolerance = float(lines[11].split()[4])
    # Deliberately out-of-tolerance probes describe their construction distance;
    # they never relax the algorithm's 1e-5 acceptance threshold.
    if not (meta[7] and meta[1] == 5 and 40 <= meta[0] <= 49) and not 0 < tolerance <= 1e-5:
        raise ValueError('Positive V2 construction tolerance exceeds 1e-5')
    return p, q, nu*nv


def require_v2_data(path):
    root = Path(path)
    stamp = root / 'dataset.json'
    if not stamp.exists():
        raise ValueError('Retired or incomplete fallback dataset; use tests/data/general-rational-v2')
    data = json.loads(stamp.read_text(encoding='utf-8'))
    sizes = dict(main=600, extra=60, extreme=320, synthetic=5)
    if data.get('contract') != 'rational-domain-v2' or data.get('sizes') != sizes:
        raise ValueError('Invalid V2 corpus contract/size')
    if data.get('representation_policy') != REPRESENTATION_POLICY:
        raise ValueError('V2 corpus does not declare the degree/control/tolerance policy')
    records = {(r['suite'], int(r['id'])): r for r in data['cases']}
    if len(data['cases']) != 985 or set(records) != {(s, i) for s, n in sizes.items() for i in range(n)}:
        raise ValueError('Incomplete/duplicate V2 case identities')
    for (suite, i), record in records.items():
        name = f'{suite}-{i:03}.b2a'
        if record['file'] != name:
            raise ValueError('Unexpected V2 filename')
        content = (root / suite / name).read_bytes()
        if not content.startswith(b'B2A_RATIONAL_DOMAIN_V2') or hashlib.sha256(content).hexdigest() != record['sha256']:
            raise ValueError('V2 input changed or retired: ' + name)
        actual = require_representation(content)
        if actual != (int(record['degree_u']), int(record['degree_v']), int(record['controls'])):
            raise ValueError('V2 representation metadata mismatch: ' + name)
    return root
