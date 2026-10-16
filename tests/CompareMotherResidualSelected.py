"""Compare the current residual certificate with the saved joint-budget executable."""
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import statistics
import subprocess

from RunMotherJointSelected import CASES, METRICS, audit, sha, write_csv

METRICS = METRICS + ["boundary_nodes", "tangent_residual_fits", "tangent_residual_splits"]


def means(rows):
    return {key: statistics.mean(float(r[key]) for r in rows if r.get(key))
            if any(r.get(key) for r in rows) else "" for key in METRICS}


def save_rows(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    write_csv(path, [{k: row.get(k, "") for k in keys} for row in rows])


def main():
    parser = argparse.ArgumentParser(__doc__)
    for name in ("exe", "baseline-exe", "baseline-manifest", "sdk", "data", "work", "output", "unit-log"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--runs", default=5, type=int)
    args = parser.parse_args()
    assert 1 <= args.runs <= 100
    for name, value in vars(args).items():
        if isinstance(value, Path): setattr(args, name, value.resolve())
    if args.output.exists() or args.work.exists():
        raise RuntimeError("Use new output/work directories; never overwrite a measurement")
    baseline = json.loads(args.baseline_manifest.read_text(encoding="utf-8"))
    assert sha(args.baseline_exe) == baseline["executable_sha256"], "Unknown baseline executable"
    inputs = {str(i): sha(args.data / f"case-{i:03}.b2a") for i in CASES}
    assert inputs == baseline["input_sha256"], "Baseline inputs changed"
    root = Path(__file__).resolve().parent.parent
    files = sorted(list((root / "src").glob("*.cpp")) + list((root / "src").glob("*.h")) +
                   list((root / "include").glob("*.h")))
    allowed = {"src/MotherTangentValidation.cpp", "src/MotherTangentValidation.h",
               "include/MotherDistanceValidation.h"}
    for p in files:
        key = p.relative_to(root).as_posix()
        if key not in allowed:
            assert sha(p) == baseline["source_sha256"][key], "Unrelated library changed: " + key
    files += [root / p for p in ["CMakeLists.txt", "tests/MotherDistanceSelectedTests.cpp",
              "tests/MotherDistanceTests.cpp", "tests/Dataset500.cpp", "tests/PairedSurfaceBounds.cpp",
              "tests/RunMotherJointSelected.py", "tests/CompareMotherResidualSelected.py"]]
    hashes = {p.relative_to(root).as_posix(): sha(p) for p in files}
    git = ["git", "-c", "safe.directory=" + root.as_posix()]
    manifest = dict(started_utc=datetime.now(timezone.utc).isoformat(), cases=CASES, runs=args.runs,
        source_base_commit=subprocess.check_output(git + ["rev-parse", "HEAD"], cwd=root, text=True).strip(),
        source_sha256=hashes, input_sha256=inputs, residual_executable_sha256=sha(args.exe),
        baseline_executable_sha256=sha(args.baseline_exe), baseline_manifest_sha256=sha(args.baseline_manifest),
        baseline_source_commit=baseline["source_commit"], commands=[], failures=[],
        contract="same strict 0.1-degree joint tangent condition and 1e-5 positions; no normal fallback",
        order="serial case blocks, alternating baseline/residual first; five cold/warm pairs per executable block")
    args.work.mkdir(parents=True)
    args.output.mkdir(parents=True)
    env = os.environ.copy()
    env["PATH"] = str(args.sdk / "x64-win/bin") + os.pathsep + env["PATH"]
    rows, logs = [], []
    for index, case in enumerate(CASES):
        order = ("baseline", "residual") if index % 2 == 0 else ("residual", "baseline")
        for method in order:
            dest = args.work / f"{method}-{case:03}.csv"
            exe = args.baseline_exe if method == "baseline" else args.exe
            command = [str(exe), str(args.data), str(case), str(dest), str(args.runs)]
            manifest["commands"].append(dict(method=method, case=case, command=command))
            try:
                result = subprocess.run(command, cwd=args.sdk, env=env, capture_output=True,
                                        text=True, errors="replace", timeout=120)
                logs.append(f"{method} case {case}, exit {result.returncode}\n{result.stdout}\n{result.stderr}")
                current = list(csv.DictReader(dest.open(encoding="utf-8", newline=""))) if dest.exists() else []
                rows.extend(dict(method=method, **r) for r in current)
                if result.returncode: raise RuntimeError(f"Exit code {result.returncode}")
                audit(current, case, args.runs)
                if method == "residual":
                    assert all(int(r["tangent_residual_fits"]) > 0 for r in current)
                print(f"{method} case {case}: {len(current)} measurements and audits passed", flush=True)
            except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError, KeyError) as exc:
                manifest["failures"].append(dict(method=method, case=case, reason=str(exc)))
                logs.append(f"{method} case {case}: FAILED {exc}")
                print(logs[-1], flush=True)
            save_rows(args.output / "measurements.csv", rows)
            (args.output / "runs.log").write_text("\n".join(logs), encoding="utf-8")
            (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (args.output / "unit.log").write_bytes(args.unit_log.read_bytes())
    assert all(sha(root / key) == value for key, value in hashes.items()), "Sources changed during run"
    if manifest["failures"]: return 2
    per_case, by_type, comparison, summary = [], [], [], {}
    for method in ("baseline", "residual"):
        group = [r for r in rows if r["method"] == method]
        cold = [r for r in group if r["cache"] == "cold"]
        for case in CASES:
            for cache in ("cold", "warm"):
                selected = [r for r in group if int(r["id"]) == case and r["cache"] == cache]
                times = [float(r["algorithm_ms" if cache == "cold" else "validation_ms"]) for r in selected]
                per_case.append(dict(method=method, id=case, type=selected[0]["type"], cache=cache,
                    runs=len(selected), controls=int(selected[0]["controls"]), **means(selected),
                    measured_total_median_ms=statistics.median(times), measured_total_min_ms=min(times),
                    measured_total_max_ms=max(times), runs_under_10ms=sum(t < 10 for t in times)))
        for kind in ("Plane", "Cylinder", "Cone", "Sphere", "Torus"):
            for cache in ("cold", "warm"):
                selected = [r for r in group if r["type"] == kind and r["cache"] == cache]
                by_type.append(dict(method=method, type=kind, cache=cache,
                    cases=len({r["id"] for r in selected}), measurements=len(selected), **means(selected)))
        summary[method] = dict(cases=len(CASES), measurements=len(group), passed=len(group),
            cold_case_means_under_10ms=sum(r["algorithm_ms"] < 10 for r in per_case if r["method"] == method and r["cache"] == "cold"),
            cold_runs_under_10ms=sum(float(r["algorithm_ms"]) < 10 for r in cold),
            cold_max_ms=max(float(r["algorithm_ms"]) for r in cold),
            max_joint_bound_deg=max(float(r["joint_tangent_bound_deg"]) for r in group),
            max_audit_joint_deg=max(float(r["audit_joint_max_deg"]) for r in group),
            tangent_residual_splits=sum(int(r["tangent_residual_splits"]) for r in group)
                if method == "residual" else None,
            boundary_fallbacks=sum(int(r["boundary_fallbacks"]) for r in group))
    for case in CASES:
        pair = {r["method"]: r for r in per_case if r["id"] == case and r["cache"] == "cold"}
        out = dict(id=case, type=pair["baseline"]["type"])
        for metric in ("algorithm_ms", "validation_ms", "boundary_angle_ms"):
            a, b = pair["baseline"][metric], pair["residual"][metric]
            out.update({"baseline_" + metric: a, "residual_" + metric: b,
                        metric + "_reduction_percent": 100 * (a - b) / a})
        comparison.append(out)
    write_csv(args.output / "per-case.csv", per_case)
    write_csv(args.output / "modules-by-type.csv", by_type)
    write_csv(args.output / "comparison.csv", comparison)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
