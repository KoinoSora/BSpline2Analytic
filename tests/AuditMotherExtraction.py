"""Independent 80-digit blossom audit of exported binary64 source envelopes.

This is a numerical regression audit, not a replacement for the analytic error
bound. Decimal inputs are the exact values of the exported binary64 numbers.
"""
import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path

def exact(x):
    return Decimal.from_float(float(x))

def axis(knots, mults, degree, first, lo, hi):
    full = [exact(k) for k, m in zip(knots, mults) for _ in range(m)]
    a, n = first + degree, degree + 1
    result = []
    # Blossom directly on the actual (possibly trimmed) domain, rather than
    # reproduce the implementation's full-span matrix followed by restriction.
    for j in range(n):
        work = [[Decimal(int(i == k)) for k in range(n)] for i in range(n)]
        for r in range(1, n):
            t = exact(lo if r <= degree-j else hi)
            for i in range(degree, r-1, -1):
                alpha = (t-full[a-degree+i]) / (full[a+1+i-r]-full[a-degree+i])
                work[i] = [(1-alpha)*x + alpha*y for x, y in zip(work[i-1], work[i])]
        result.append(work[degree])
    return result

def audit(path):
    data = json.loads(path.read_text())
    p, q = data['p'], data['q']
    origin = list(map(exact, data['origin']))
    controls = [[[exact(x)-o for x, o in zip(point, origin)] for point in row] for row in data['controls']]
    count, maximum = 0, Decimal(0)
    for cell in data['cells']:
        u0, u1, v0, v1 = cell['domain']
        fu, fv = cell['first_u'], cell['first_v']
        u = axis(data['knots_u'], data['mults_u'], p, fu, u0, u1)
        v = axis(data['knots_v'], data['mults_v'], q, fv, v0, v1)
        # Contract V first; production contracts U first and reuses strips.
        strip = [[[sum(v[j][l]*controls[fu+k][fv+l][d] for l in range(q+1))
                   for d in range(3)] for j in range(q+1)] for k in range(p+1)]
        for i in range(p+1):
            for j in range(q+1):
                for d in range(3):
                    reference = sum(u[i][k]*strip[k][j][d] for k in range(p+1))
                    observed = exact(cell['coefficients'][i*(q+1)+j][d])
                    error = exact(cell['error'][d])
                    delta = abs(reference-observed)
                    if delta > error:
                        raise AssertionError(f'{path.name}: coefficient {(i,j,d)} error {delta} > {error}')
                    if error:
                        maximum = max(maximum, delta/error)
                    count += 1
    return count, float(maximum)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', type=Path)
    args = parser.parse_args()
    with localcontext() as ctx:
        ctx.prec = 80
        total, maximum = 0, 0.0
        for path in args.files:
            count, ratio = audit(path)
            total += count
            maximum = max(maximum, ratio)
            print(f'{path.name}: coordinates={count} max_error/envelope={ratio:.6g} PASSED')
        print(f'TOTAL coordinates={total} max_error/envelope={maximum:.6g} PASSED')

if __name__ == '__main__':
    main()
