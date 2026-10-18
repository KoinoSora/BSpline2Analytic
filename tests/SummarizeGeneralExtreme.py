# 文件职责：历史一般极端记录的结构、尺寸与性能汇总。
# 控制凸包证明与弦长测量分别表述，不重跑归档输入。

"""只审计历史 GeneralExtreme320 成员、尺寸、哈希与测试记录。"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

QUOTAS = dict(radius_arc=64, cylinder_helix=32, slender_cylinder_helix=32,
              sphere_triangle=48, cone_helix=32, slender_cone=32,
              torus_extreme=32, world_coordinates=48)
SUCCESS = {"same_type_success", "alternate_type_success"}
WORLD_COORDINATE_LIMIT = 1e6


def read(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def fixture(path):
    tokens = iter(path.read_text().split())
    assert next(tokens) == "B2A_GENERAL_V1"
    metadata = [int(next(tokens)) for _ in range(9)]
    origin, axis, xdir = ([float(next(tokens)) for _ in range(3)] for _ in range(3))
    kind = int(next(tokens))
    radius, major, angle = (float(next(tokens)) for _ in range(3))
    uv = [[float(next(tokens)), float(next(tokens))] for _ in range(4)]
    extra = [float(next(tokens)) for _ in range(5)]
    p, q, nu, nv, rational = (int(next(tokens)) for _ in range(5))
    for degree, count in ((p, nu), (q, nv)):
        count_knots = int(next(tokens))
        knots = [(float(next(tokens)), int(next(tokens))) for _ in range(count_knots)]
        assert knots[0][0] == 0 and knots[-1][0] == 1
        assert all(a[0] < b[0] for a, b in zip(knots, knots[1:]))
        assert sum(m for _, m in knots) == count+degree+1
        assert all(0 < m <= degree for _, m in knots[1:-1])
    points, weights = [], []
    for _ in range(nu*nv):
        points.append([float(next(tokens)) for _ in range(3)])
        weights.append(float(next(tokens)) if rational else 1.0)
    assert next(tokens, None) is None
    assert all(math.isfinite(v) for point in points for v in point)
    assert all(math.isfinite(w) and w > 0 for w in weights)
    assert bool(metadata[6]) == (max(weights) != min(weights))
    assert 50 <= nu*nv <= 200 and 0 < p < nu and 0 < q < nv
    return dict(metadata=metadata, origin=origin, axis=axis, xdir=xdir, kind=kind,
                radius=radius, major=major, angle=angle, uv=uv, extra=extra,
                count=nu*nv, ulp=max(math.ulp(v) for point in points for v in point),
                coordinate_min=[min(point[d] for point in points) for d in range(3)],
                coordinate_max=[max(point[d] for point in points) for d in range(3)])


def audit_data(data):
    rows = read(data/"manifest.csv")
    assert [int(r["id"]) for r in rows] == list(range(320)), "Incomplete or reordered matrix"
    assert dict(Counter(r["group"] for r in rows)) == QUOTAS
    assert Counter(r["rational"] for r in rows) == {"0": 160, "1": 160}
    hashes, previous = {}, None
    for row in rows:
        path = data/row["file"]
        assert path.name == f"extreme-{int(row['id']):03}.b2a"
        f = fixture(path)
        # fixture() 已检查正权重，因此控制凸包直接包围整个曲面，
        # 包括四边，无需依赖采样推断覆盖。
        assert all(-WORLD_COORDINATE_LIMIT <= lo <= hi <= WORLD_COORDINATE_LIMIT
                   for lo, hi in zip(f["coordinate_min"], f["coordinate_max"])), (row["id"], "world bounds")
        if row["group"] == "world_coordinates":
            kind, slot = divmod(int(row["id"])-272, 8)
            near_limit = (1 if kind % 2 == 0 else -1)*(WORLD_COORDINATE_LIMIT-2000)
            translation = (0.0, 1e3, 1e5, near_limit)[slot//2]
            assert float(row["translation"]) == translation, (row["id"], "outdated translation")
            assert f["origin"] == [translation, -.73*translation, .41*translation], (row["id"], "origin mismatch")
        assert f["metadata"][0] == int(row["id"]) and f["count"] == int(row["controls"])
        assert f["metadata"][8] == int(row["collapsed_mask"])
        if row["type"] in {"Cylinder", "Sphere", "Torus"}:
            assert .01 <= f["radius"] <= 1000
        if row["type"] == "Torus":
            assert f["radius"] < f["major"] <= 1000
        if row["type"] == "Cone":
            assert 0 < f["angle"] < math.pi/2
            assert all(0 <= p[1] <= 1000 and (p[1] == 0 or .01-1e-12 <= p[1]*math.tan(f["angle"]) <= 1000+1e-10) for p in f["uv"])
        for key in ("length", "width"):
            assert .01-1e-12 <= float(row[key]) <= 1000+1e-10, (row["id"], key)
        # 边长采用 1024/2048 段弦长测量，不是区间证明。
        allowance = 2*float(row["edge_measure_change"])+1e-8
        assert float(row["edge_min"]) >= .01-allowance, (row["id"], "short boundary")
        assert float(row["edge_max"]) <= 1000+allowance, (row["id"], "long boundary")
        assert float(row["coordinate_ulp"]) == f["ulp"]
        assert float(row["source_bound"]) == f["extra"][1]
        if row["qualified"] == "1":
            tolerance = min(1e-6, .1*float(row["sag"])) if float(row["sag"]) > 0 else 1e-6
            assert 0 <= f["extra"][1] < tolerance
        if int(row["id"]) % 2:
            for key in ("origin", "axis", "xdir", "kind", "radius", "major", "angle", "uv"):
                assert f[key] == previous[key], (row["id"], "pair mismatch", key)
        previous = f
        hashes[row["id"]] = hashlib.sha256(path.read_bytes()).hexdigest()
    return rows, hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    manifest, hashes = audit_data(args.data)
    report = {"count": 320, "groups": QUOTAS, "qualified": sum(r["qualified"] == "1" for r in manifest),
              "controls_min": min(int(r["controls"]) for r in manifest),
              "controls_max": max(int(r["controls"]) for r in manifest), "input_sha256": hashes}
    print(f"320 fixtures; rational/polynomial 160/160; dimensions and membership passed; qualified={report['qualified']}")
    failures = 320-report["qualified"]
    if args.results:
        rows = read(args.results)
        assert [int(r["id"]) for r in rows] == list(range(320)), "Incomplete results"
        assert all(r["input_sha256"] == hashes[r["id"]] for r in rows), "Input/result hash mismatch"
        for r in rows:
            if r["status"] in SUCCESS:
                assert r["qualified"] == "1"
                assert all(0 <= float(r[k]) < 1e-5 for k in ("output_bound", "geometry_bound", "boundary_bound"))
            if math.isfinite(float(r["total_ms"])):
                assert abs(float(r["total_ms"])-sum(float(r[k]) for k in ("classify_ms", "optimize_ms", "validate_ms"))) < max(1, float(r["total_ms"])*.01)
        report["statuses"] = dict(Counter(r["status"] for r in rows))
        report["known_stage_seconds"] = {k: sum(float(r[k]) for r in rows if math.isfinite(float(r[k])))/1000
                                         for k in ("classify_ms", "optimize_ms", "validate_ms", "total_ms")}
        report["missing_stage_count"] = sum(not math.isfinite(float(r["total_ms"])) for r in rows)
        groups = defaultdict(list)
        for row in rows:
            groups[row["group"]].append(row)
        report["results_by_group"] = {}
        print("| Group | Inputs | Success | Other type | Unfinished timing | Known max classify+optimize ms |")
        print("| --- | ---: | ---: | ---: | ---: | ---: |")
        for name, group in groups.items():
            finite = [r for r in group if math.isfinite(float(r["total_ms"]))]
            stats = dict(count=len(group), success=sum(r["status"] in SUCCESS for r in group),
                         alternate=sum(r["status"] == "alternate_type_success" for r in group),
                         missing=len(group)-len(finite),
                         known_max_solve_ms=max((float(r["classify_ms"])+float(r["optimize_ms"]) for r in finite), default=None),
                         statuses=dict(Counter(r["status"] for r in group)))
            report["results_by_group"][name] = stats
            print(f"| {name} | {stats['count']} | {stats['success']} | {stats['alternate']} | {stats['missing']} | {stats['known_max_solve_ms']} |")
        failures = sum(r["status"] not in SUCCESS for r in rows)
        print("Statuses:", report["statuses"])
        print("Known stage totals (s):", report["known_stage_seconds"])
        print(f"Passed {320-failures}/320; missing timing={report['missing_stage_count']}")
    if args.json:
        args.json.write_text(json.dumps(report, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    return int(failures != 0)


if __name__ == "__main__":
    raise SystemExit(main())
