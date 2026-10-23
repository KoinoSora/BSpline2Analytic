"""串行比较保存的母面基线与圆锥因子版本；近锥顶探针另计位置模块耗时。"""
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import statistics
import subprocess

from RunMotherJointSelected import CASES, METRICS, audit, sha, write_csv


def main():
    parser = argparse.ArgumentParser(__doc__)
    for name in ("baseline", "exe", "probe", "sdk", "data", "work", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--baseline-sha256", required=True)
    parser.add_argument("--baseline-commit", required=True)
    parser.add_argument("--runs", default=5, type=int)
    args = parser.parse_args()
    assert 1 <= args.runs <= 100
    for name, value in vars(args).items():
        if isinstance(value, Path): setattr(args, name, value.resolve())
    assert sha(args.baseline).lower() == args.baseline_sha256.lower(), "Baseline executable changed"
    if args.work.exists() or args.output.exists():
        raise RuntimeError("Use new work/output directories; measurements must not be overwritten")
    args.work.mkdir(parents=True); args.output.mkdir(parents=True)
    env = os.environ.copy()
    env["PATH"] = str(args.sdk / "x64-win/bin") + os.pathsep + env["PATH"]
    # Dataset500 每组最后四例为坍缩边；仅选择 40 个常规非有理圆锥。
    cases = sorted(set(CASES + [i for i in range(201, 300, 2) if i % 20 < 16]))
    root = Path(__file__).resolve().parent.parent
    files = [root / "CMakeLists.txt", Path(__file__)]
    files += sorted((root / "src").glob("*.*")) + sorted((root / "include").glob("*.h"))
    files += [root / "tests/MotherConeFactorTests.cpp", root / "tests/MotherDistanceSelectedTests.cpp"]
    manifest = dict(started_utc=datetime.now(timezone.utc).isoformat(), cases=cases, runs=args.runs,
        baseline_commit=args.baseline_commit, baseline_sha256=sha(args.baseline), exe_sha256=sha(args.exe),
        probe_sha256=sha(args.probe), source_sha256={str(f.relative_to(root)): sha(f) for f in files if f.is_file()},
        input_sha256={str(i): sha(args.data / f"case-{i:03}.b2a") for i in cases},
        commands=[], failures=[], scope="55 fixed cases: classification + initial validation, cold/warm separated; synthetic probes: prepared-source position only")
    rows, logs = [], []
    def save():
        write_csv(args.output / "measurements.csv", rows)
        (args.output / "runs.log").write_text("\n".join(logs), encoding="utf-8")
        (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    def run(command):
        manifest["commands"].append(command)
        result = subprocess.run(command, cwd=args.sdk, env=env, capture_output=True, text=True,
                                errors="replace", timeout=120)
        logs.append(" ".join(command) + f"\nexit={result.returncode}\n{result.stdout}\n{result.stderr}")
        if result.returncode: raise RuntimeError(f"Command failed: {command[0]}")
    try:
        for index, case in enumerate(cases):
            for method in (("baseline", "factor") if index % 2 == 0 else ("factor", "baseline")):
                path = args.work / f"{method}-{case}.csv"
                exe = args.baseline if method == "baseline" else args.exe
                run([str(exe), str(args.data), str(case), str(path), str(args.runs)])
                current = list(csv.DictReader(path.open(encoding="utf-8")))
                audit(current, case, args.runs)
                rows.extend(dict(method=method, **r) for r in current)
            print(f"case {case}: baseline/factor validated", flush=True)
            save()
        run([str(args.probe)])
        run([str(args.probe), "--benchmark", str(args.output / "position-probes.csv")])
    except Exception as error:
        manifest["failures"].append(str(error)); save(); raise
    summary = []
    for method in ("baseline", "factor"):
        for cache in ("cold", "warm"):
            for kind in sorted({r["type"] for r in rows}):
                group = [r for r in rows if (r["method"], r["cache"], r["type"]) == (method, cache, kind)]
                item = dict(method=method, cache=cache, type=kind, cases=len({r["id"] for r in group}), measurements=len(group))
                for key in METRICS:
                    values = [float(r[key]) for r in group if r[key]]
                    item[key] = statistics.mean(values) if values else ""
                summary.append(item)
    write_csv(args.output / "summary.csv", summary)
    probes = list(csv.DictReader((args.output / "position-probes.csv").open()))
    probe_summary = []
    for case in range(4):
        for enabled in ("0", "1"):
            group = [r for r in probes if (r["case"], r["enabled"]) == (str(case), enabled)]
            expected = case == 0 or (enabled == "1" and case in (1, 2))
            assert len(group) == 100 and all(bool(int(r["passed"])) == expected for r in group)
            item = dict(case=case, enabled=enabled, measurements=len(group), passed=sum(int(r["passed"]) for r in group))
            for key in ("total_ms", "position_ms", "subdivision_ms", "build_ms", "factor_ms", "nodes", "attempts", "accepted_nodes"):
                item[key] = statistics.mean(float(r[key]) for r in group)
            item["max_ms"] = max(float(r["total_ms"]) for r in group)
            probe_summary.append(item)
    write_csv(args.output / "position-summary.csv", probe_summary)
    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    save()
    print(f"PASS: {len(rows)} fixed-case validations and {len(probes)} position probes")


if __name__ == "__main__":
    main()
