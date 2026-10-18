# 文件职责：业务极端矩阵的数据结构、尺寸及配额审计。
# 读取保存的输入和清单，验证数据构造合同，不运行转换器。

"""独立审计已保存业务矩阵的结构、尺寸与配额。
不调用 Convert、不重新生成数据，也不替代 SDK 侧连续配对距离检查。
保存的证书与采样角度仍按各自保证范围解释；可选输出包含每个输入的 SHA-256。"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


def read_case(path):
    tokens = iter(path.read_text(encoding="utf-8").split())
    assert next(tokens) == "B2A_CASE_V1", path
    meta = [next(tokens) for _ in range(9)]
    frame = [float(next(tokens)) for _ in range(9)]
    radius, major, angle = [float(next(tokens)) for _ in range(3)]
    domain = [float(next(tokens)) for _ in range(4)]
    surfaces = []
    for _ in range(2):
        pu, pv, nu, nv, weighted = [int(next(tokens)) for _ in range(5)]
        for degree, count in ((pu, nu), (pv, nv)):
            n = int(next(tokens))
            knots = [(float(next(tokens)), int(next(tokens))) for _ in range(n)]
            assert sum(m for _, m in knots) == count + degree + 1
            assert all(math.isfinite(t) and m > 0 for t, m in knots)
            assert all(a[0] < b[0] for a, b in zip(knots, knots[1:]))
        points = []
        weights = []
        for _ in range(nu * nv):
            points.append(tuple(float(next(tokens)) for _ in range(3)))
            weights.append(float(next(tokens)) if weighted else 1.0)
        assert all(math.isfinite(x) for p in points for x in p)
        assert all(math.isfinite(w) and w > 0 for w in weights)
        surfaces.append((nu, nv, points, weights))
    assert next(tokens, None) is None
    return meta, frame, (radius, major, angle), domain, surfaces


def audit(directory):
    directory = Path(directory)
    assert (directory / "profile.txt").read_text().strip() == "business-extreme200-v1"
    with (directory / "manifest.csv").open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    assert [int(r["id"]) for r in rows] == list(range(200))
    assert len(list(directory.glob("extreme-*.b2a"))) == 200
    def size(x):
        return math.isfinite(x) and .01 * (1 - 1e-10) <= x <= 1000 * (1 + 1e-10)
    actual_counts, bands = Counter(), Counter()
    for row in rows:
        i = int(row["id"])
        meta, frame, params, d, surfaces = read_case(directory / f"extreme-{i:03}.b2a")
        t = int(meta[2]); mask = int(meta[7]); r, R, alpha = params
        assert int(meta[0]) == i and int(meta[1]) == i // 2 and int(meta[6]) == i % 2
        assert int(row["collapsed_mask"]) == mask
        assert (row["group"] == "collapsed") == (mask != 0)
        assert row["type"] == ("Plane", "Cylinder", "Cone", "Sphere", "Torus")[t]
        if t == 0:
            assert size(d[1]-d[0]) and size(d[3]-d[2])
        elif t == 1:
            assert size(r) and size(d[3]-d[2])
        elif t == 2:
            assert (d[2] == 0 if mask else size(d[2]*math.tan(alpha)))
            assert size(d[3]*math.tan(alpha)) and size(d[3]-d[2]) and d[3] <= 1000*(1+1e-10)
        elif t == 3:
            assert size(r)
        else:
            assert size(r) and size(R) and R > r
        nu, nv, points, weights = surfaces[1]
        count = nu * nv
        assert count == int(row["controls"]) and 50 <= count <= 200
        if i % 2:
            assert all(w == 1 for w in weights) and row["rational"] == "0"
        actual_counts[count] += 1
        bands[next(k for k, (lo, hi) in enumerate(((50,79),(80,109),(110,139),(140,169),(170,200))) if lo <= count <= hi)] += 1
        for e in range(4):
            if mask & (1 << e):
                for u, v, pts, _ in surfaces:
                    edge = [pts[(0 if e == 0 else u-1)*v+j] for j in range(v)] if e < 2 else [pts[j*v+(0 if e == 2 else v-1)] for j in range(u)]
                    assert all(p == edge[0] for p in edge)
            else:
                assert size(float(row[f"edge{e}"]))
        assert row["qualified"] == row["dimensions_passed"] == "1"
        assert float(row["pair_bound"]) < 1e-6
        for key in ("source_tangent_sampled_deg", "source_normal_sampled_deg"):
            value = float(row[key])
            assert math.isnan(value) if mask else math.isfinite(value) and value < .05
    assert bands == Counter({0:40, 1:40, 2:40, 3:40, 4:40})
    assert actual_counts[50] and actual_counts[200]
    assert sum(r["collapsed_mask"] != "0" for r in rows) == 8
    files = sorted(directory.glob("*.b2a")) + [directory/"manifest.csv", directory/"profile.txt"]
    finite_max = lambda field: max(float(r[field]) for r in rows if math.isfinite(float(r[field])))
    return {
        "profile": "business-extreme200-v1", "count": len(rows),
        "source_counts": dict(Counter(r["source_kind"] for r in rows)),
        "type_counts": dict(Counter(r["type"] for r in rows)),
        "scenario_counts": dict(Counter(r["group"] for r in rows)),
        "control_counts": dict(sorted(actual_counts.items())), "control_bands": dict(sorted(bands.items())),
        "degenerate_inputs": 8, "pair_bound_max": finite_max("pair_bound"),
        "source_tangent_sampled_deg_max": finite_max("source_tangent_sampled_deg"),
        "source_normal_sampled_deg_max": finite_max("source_normal_sampled_deg"),
        "sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.directory)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k != "sha256"}, indent=2))
