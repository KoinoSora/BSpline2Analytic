# 文件职责：当前母面联合切向验证的定向计时与审计。
# 串行收集冷/热缓存结果，校验联合条件及阶段计时，失败进程也保留部分证据。

"""串行测量当前母面验证器的指定案例，分别记录冷/热缓存。"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
from datetime import datetime, timezone

CASES = [1, 101, 201, 301, 351, 381, 383, 391, 395, 401, 451, 481, 485, 489, 491, 497]
STAGES = ["preparation_ms", "position_ms", "subdivision_ms", "boundary_ms", "reporting_ms"]
METRICS = ["classification_ms", "validation_ms", "algorithm_ms"] + STAGES + [
    "source_extraction_ms", "source_boundary_ms", "target_preparation_ms", "mapping_ms", "sample_ms",
    "residual_build_ms", "residual_bound_ms", "boundary_position_ms", "boundary_angle_ms",
    "extraction_plan_ms", "source_contraction_ms", "source_envelope_ms"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def audit(rows, case, runs):
    assert len(rows) == 2 * runs, "Incomplete measurements"
    assert {(int(r["run"]), r["cache"]) for r in rows} == {
        (i, cache) for i in range(runs) for cache in ("cold", "warm")}, "Duplicate/missing run"
    for r in rows:
        assert int(r["id"]) == case and r["status"] == "Passed"
        assert int(r["degree_u"]) <= 5 and int(r["degree_v"]) <= 5 and int(r["controls"]) < 200
        t = {k: float(r[k]) for k in METRICS if r[k]}
        assert all(math.isfinite(x) and x >= 0 for x in t.values()), "Invalid timing"
        assert abs(sum(t[k] for k in STAGES) - t["validation_ms"]) < 1e-7
        if r["cache"] == "cold":
            assert abs(t["classification_ms"] + t["validation_ms"] - t["algorithm_ms"]) < 1e-7
            assert int(r["source_cache_hits"]) == 0
        else:
            assert not r["classification_ms"] and not r["algorithm_ms"]
            assert int(r["source_cache_hits"]) == 1
        assert float(r["mother_bound"]) < 1e-5 and float(r["edge_bound"]) < 1e-5
        for field in ("tangent_u_bound_deg", "tangent_v_bound_deg", "joint_tangent_bound_deg"):
            assert 0 <= float(r[field]) < .1
        assert 0 <= float(r["joint_sine_squared_bound"]) < math.sin(math.radians(.1)) ** 2
        assert not r["joint_tangent_failure"]
        assert float(r["audit_position_max"]) <= float(r["mother_bound"]) + 1e-11
        assert float(r["audit_edge_max"]) <= float(r["edge_bound"]) + 1e-10
        for axis in ("u", "v"):
            assert float(r[f"audit_tangent_{axis}_max_deg"]) <= float(r[f"tangent_{axis}_bound_deg"]) + 1e-7
        assert float(r["audit_joint_max_deg"]) <= float(r["joint_tangent_bound_deg"]) + 1e-7
        assert t["boundary_position_ms"] + t["boundary_angle_ms"] <= t["boundary_ms"] + 1e-7


def averages(group):
    out = {}
    for metric in METRICS:
        values = [float(r[metric]) for r in group if r[metric]]
        out[metric] = statistics.mean(values) if values else ""
    return out


def main():
    parser = argparse.ArgumentParser(__doc__)
    for name in ("exe", "sdk", "data", "work", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--runs", type=int, default=5)
    parser.add_argument("--cases", type=int, nargs="+", default=CASES)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--unit-log", type=Path)
    args = parser.parse_args()
    assert 1 <= args.runs <= 100 and len(set(args.cases)) == len(args.cases)
    for name in ("exe", "sdk", "data", "work", "output"):
        setattr(args, name, getattr(args, name).resolve())
    args.work.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parent.parent
    source_files = sorted(list((root / "src").glob("*.cpp")) + list((root / "src").glob("*.h")) +
                          list((root / "include").glob("*.h")))
    source_files += [root / "CMakeLists.txt", Path(__file__).resolve(),
                     root / "tests/MotherDistanceSelectedTests.cpp", root / "tests/Dataset500.cpp",
                     root / "tests/PairedSurfaceBounds.cpp"]
    git = ["git", "-c", "safe.directory=" + root.as_posix()]
    manifest = dict(started_utc=datetime.now(timezone.utc).isoformat(),
                    source_commit=subprocess.check_output(git + ["rev-parse", "HEAD"], cwd=root, text=True).strip(),
                    cases=args.cases, runs=args.runs, executable_sha256=sha(args.exe),
                    source_sha256={p.relative_to(root).as_posix(): sha(p) for p in source_files},
                    input_sha256={str(i): sha(args.data / f"case-{i:03}.b2a") for i in args.cases},
                    commands=[], failures=[], contract="joint U/V tangent budget strictly below 0.1 degree; no normal fallback")
    env = os.environ.copy()
    env["PATH"] = str(args.sdk / "x64-win/bin") + os.pathsep + env["PATH"]
    rows, logs = [], []
    for case in args.cases:
        dest = args.work / f"joint-{case:03}.csv"
        command = [str(args.exe), str(args.data), str(case), str(dest), str(args.runs)]
        manifest["commands"].append(command)
        if dest.exists():
            raise RuntimeError(f"Refusing to mix an existing measurement: {dest}")
        try:
            result = subprocess.run(command, cwd=args.sdk, env=env, capture_output=True,
                                    text=True, errors="replace", timeout=args.timeout)
            logs.append(f"case {case}, exit {result.returncode}\n{result.stdout}\n{result.stderr}")
            current = list(csv.DictReader(dest.open(encoding="utf-8", newline=""))) if dest.exists() else []
            rows.extend(current)  # 即使子进程失败，也保留已获得的部分证据。
            if result.returncode:
                raise RuntimeError(f"Process exited {result.returncode}; see runs.log")
            audit(current, case, args.runs)
            print(f"case {case}: {len(current)} measurements and independent audits passed", flush=True)
        except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError, KeyError) as exc:
            manifest["failures"].append(dict(id=case, reason=str(exc)))
            logs.append(f"case {case}: {exc}")
            print(f"case {case}: FAILED {exc}", flush=True)
        write_csv(args.output / "measurements.csv", rows)
        (args.output / "runs.log").write_text("\n".join(logs), encoding="utf-8")
        (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if args.unit_log:
        (args.output / "unit.log").write_bytes(args.unit_log.read_bytes())
    if manifest["failures"]:
        return 2
    per_case, by_type = [], []
    for case in args.cases:
        for cache in ("cold", "warm"):
            group = [r for r in rows if int(r["id"]) == case and r["cache"] == cache]
            metric = "algorithm_ms" if cache == "cold" else "validation_ms"
            values = [float(r[metric]) for r in group]
            per_case.append(dict(id=case, type=group[0]["type"], cache=cache, runs=len(group),
                degree_u=int(group[0]["degree_u"]), degree_v=int(group[0]["degree_v"]), controls=int(group[0]["controls"]),
                **averages(group), measured_total_median_ms=statistics.median(values),
                measured_total_min_ms=min(values), measured_total_max_ms=max(values),
                runs_under_10ms=sum(x < 10 for x in values)))
    for surface in ("Plane", "Cylinder", "Cone", "Sphere", "Torus"):
        for cache in ("cold", "warm"):
            group = [r for r in rows if r["type"] == surface and r["cache"] == cache]
            if group:
                by_type.append(dict(type=surface, cache=cache, cases=len({r["id"] for r in group}),
                                    measurements=len(group), **averages(group)))
    cold = [r for r in rows if r["cache"] == "cold"]
    summary = dict(cases=len(args.cases), measurements=len(rows), passed=len(rows),
        cold_runs_under_10ms=sum(float(r["algorithm_ms"]) < 10 for r in cold), cold_runs=len(cold),
        cold_case_means_under_10ms=sum(r["algorithm_ms"] < 10 for r in per_case if r["cache"] == "cold"),
        cold_max_ms=max(float(r["algorithm_ms"]) for r in cold),
        max_joint_bound_deg=max(float(r["joint_tangent_bound_deg"]) for r in rows),
        max_audit_joint_deg=max(float(r["audit_joint_max_deg"]) for r in rows),
        boundary_fallbacks=sum(int(r["boundary_fallbacks"]) for r in rows),
        maximum_mapping_attempts=max(int(r["mapping_attempts"]) for r in rows),
        maximum_surface_nodes=max(int(r["surface_nodes"]) for r in rows),
        scope="classification compute plus first candidate validation; no optimizer run; audits excluded",
        historical_comparison="not a paired performance comparison; older acceptance contracts differ")
    write_csv(args.output / "per-case.csv", per_case)
    write_csv(args.output / "modules-by-type.csv", by_type)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
