# 文件职责：源 Bézier 提取包络的独立高精度审计。
# 采用不同收缩顺序和直接裁剪域提取核对双精度误差包络，不以生产实现自身作参考答案。

"""以 80 位精度独立审计导出的 binary64 源面提取包络。
这是数值回归审计，不替代解析误差上界；Decimal 输入精确表示导出的二进制双精度数。"""
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
    # 在实际参数域（可含裁剪）上直接进行多仿射提取，
    # 避免复刻生产实现先构造整区间矩阵再限制的计算路径。
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
        # 审计先沿 V 收缩；生产实现先沿 U 收缩并复用条带，以保持实现独立性。
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
