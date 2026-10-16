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
        ratio = F(raw.pop('angle_ratio'))
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
        dot_axes = [product(a[k], b[k]) for k in range(3)]
        dot_lo = min(sum(x[0] for x in coefficient) for coefficient in zip(*dot_axes))
        cross = [difference(product(a[(k+1)%3], b[(k+2)%3]), product(a[(k+2)%3], b[(k+1)%3]))
                 for k in range(3)]
        upper_square = sum(max(abs(x) for pair in axis for x in pair)**2 for axis in cross)
        assert dot_lo > 0 and ratio*ratio*dot_lo*dot_lo >= upper_square, (case_id, 'angle understated')
    print(f'PASS: {len(cases)} cases; {count} exact rational coefficient enclosures; '
          f'{len(cases)} dot/cross bounds; no failures')


if __name__ == '__main__':
    main(sys.argv[1])
