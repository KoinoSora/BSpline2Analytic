"""Exact rational audit of double error envelopes (no SDK/validator calls).

Treat every exported binary double endpoint as exact. Compute independent
interval Bernstein products with fractions.Fraction and exact binomial ratios.
This checks coefficients, not only sampled values or double-precision results.
"""
import json
import math
import sys
from fractions import Fraction as F


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def sub(a, b):
    return a[0] - b[1], a[1] - b[0]


def mul(a, b):
    x = [u * v for u in a for v in b]
    return min(x), max(x)


def product(a, b):
    p, q = len(a)-1, len(b)-1
    out = [(F(0), F(0)) for _ in range(p+q+1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            f = F(math.comb(p, i)*math.comb(q, j), math.comb(p+q, i+j))
            out[i+j] = add(out[i+j], mul(mul(x, y), (f, f)))
    return out


def difference(a, b):
    n = max(len(a), len(b))
    a = product(a, [(F(1), F(1))]*(n-len(a)+1))
    b = product(b, [(F(1), F(1))]*(n-len(b)+1))
    return [sub(x, y) for x, y in zip(a, b)]


def derivative(a):
    n = len(a)-1
    return [mul(sub(a[i+1], a[i]), (F(n), F(n))) for i in range(n)] or [(F(0), F(0))]


def main(path):
    cases = json.load(open(path, encoding='utf-8'))
    count = 0
    for case_id, raw in enumerate(cases):
        if 'angle_sine_squared' not in raw:
            raise ValueError('Residual envelope export required; historical dot/cross exports are archived')
        sine_squared = F(raw.pop('angle_sine_squared'))
        scale = [tuple(F(x) for x in pair) for pair in raw.pop('angle_scale')]
        reference = [F(x) for x in raw.pop('angle_reference')]
        case = {key: [[tuple(F(v) for v in pair) for pair in axis] for axis in cell]
                for key, cell in raw.items()}
        s, t = case['source'], case['target']
        expected = [difference(product(t[k], s[3]), product(s[k], t[3])) for k in range(3)]
        expected.append(product(s[3], t[3]))
        tangent = [difference(product(derivative(s[k]), s[3]), product(s[k], derivative(s[3])))
                   if len(s[3]) > 1 else product(derivative(s[k]), s[3]) for k in range(3)]
        for label, axes in [('error', expected), ('tangent', tangent)]:
            for k, axis in enumerate(axes):
                actual = case[label][k]
                assert len(actual) == len(axis), (case_id, label, k, 'degree')
                for i, ((lo, hi), (a, b)) in enumerate(zip(actual, axis)):
                    assert lo <= a <= b <= hi, (case_id, label, k, i, 'outside envelope')
                    count += 1
        a, b = case['angle_source'], case['angle_target']
        assert all(0 < lo <= hi for lo, hi in scale), (case_id, 'nonpositive scale')
        residual = case['angle_residual'][:3]
        for k in range(3):
            expected = difference(a[k], product(b[k], scale))
            assert len(expected) == len(residual[k]), (case_id, 'residual degree')
            for (lo, hi), (x, y) in zip(residual[k], expected):
                assert lo <= x <= y <= hi, (case_id, 'residual coefficient outside envelope')
                count += 1
        assert len({len(axis) for axis in a[:3]}) == 1
        projection_lo = min(sum(mul(coefficient[k], (reference[k], reference[k]))[0]
                                for k in range(3)) for coefficient in zip(*a[:3]))
        reference_squared = sum(x*x for x in reference)
        assert projection_lo > 0 and reference_squared > 0, (case_id, 'source length unproved')
        length_squared_lo = projection_lo**2 / reference_squared
        if len({len(axis) for axis in residual}) == 1:
            error_squared = max(sum(max(abs(x) for x in component)**2 for component in coefficient)
                                for coefficient in zip(*residual))
        else:
            error_squared = sum(max(abs(x) for pair in axis for x in pair)**2 for axis in residual)
        assert 0 <= sine_squared < 1 and sine_squared*length_squared_lo >= error_squared, (case_id, 'angle understated')
    print(f'PASS: {len(cases)} cases; {count} exact rational coefficient enclosures; '
          f'{len(cases)} positive-scale residual bounds; no failures')


if __name__ == '__main__':
    main(sys.argv[1])
