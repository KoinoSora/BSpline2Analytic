# 文件职责：等参增量分类的代表例串行版本对照。
# 检查证书与候选尝试统计，完整保留来源和输出审计，区别冷上下文与共享源缓存。

"""只在 16 个代表等参案例上串行比较版本。"""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import statistics
import subprocess

from RunMotherJointSelected import CASES, METRICS as BASE_METRICS, audit, sha, write_csv

CLASS_TIMES = ["class_sampling_ms", "class_algebraic_ms", "class_geometric_ms", "class_scoring_ms", "class_recovery_ms", "class_other_ms"]
COUNTERS = ["validation_attempts", "curve_extractions", "curve_point_evaluations", "score_point_evaluations", "point_cache_hits", "fit_passes", "geometric_iterations", "fast_circle_fits", "sampling_refinements", "classification_stages", "raw_candidates", "distance_evaluations", "shared_mother_scores", "duplicate_candidates", "patch_recovery_attempts"]
METRICS = BASE_METRICS + CLASS_TIMES + COUNTERS


def means(rows):
    return {k: statistics.mean(float(r[k]) for r in rows if r.get(k))
            if any(r.get(k) for r in rows) else "" for k in METRICS}


def save(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    write_csv(path, [{k: r.get(k, "") for k in keys} for r in rows])


def check_incremental(rows, attempts, case, runs):
    # 保留原有证书、源数据合同、计时与密集审计检查。
    # 只调整冷缓存命中次数的预期：同一调用中的失败候选也共享源面准备。
    normalized = [dict(r) for r in rows]
    for r in normalized:
        if r["cache"] == "cold": r["source_cache_hits"] = "0"
    audit(normalized, case, runs)
    for r in rows:
        if r["cache"] != "cold": continue
        parts = [a for a in attempts if int(a["run"]) == int(r["run"])]
        assert len(parts) == int(r["validation_attempts"]) >= 1
        assert [int(a["attempt"]) for a in parts] == list(range(1,len(parts)+1))
        assert parts[-1]["status"] == "Passed"
        assert sum(a["status"] == "Passed" for a in parts) == 1
        assert int(r["source_cache_hits"]) == len(parts)-1
        for k in ["validation_ms", "preparation_ms", "position_ms", "subdivision_ms", "boundary_ms", "reporting_ms"]:
            assert abs(sum(float(a[k]) for a in parts)-float(r[k])) < 1e-7
        assert abs(sum(float(r[k]) for k in CLASS_TIMES)-float(r["classification_ms"])) < 1e-7
        assert all(math.isfinite(float(r[k])) and float(r[k]) >= 0 for k in CLASS_TIMES+COUNTERS)


def main():
    parser=argparse.ArgumentParser(__doc__)
    for name in ("exe", "baseline-exe", "baseline-manifest", "source", "sdk", "data", "work", "output"):
        parser.add_argument("--"+name, required=True, type=Path)
    parser.add_argument("--runs",type=int,default=5)
    args=parser.parse_args()
    for k,v in vars(args).items():
        if isinstance(v,Path): setattr(args,k,v.resolve())
    assert 1<=args.runs<=100
    if args.work.exists() or args.output.exists(): raise RuntimeError("Use new output/work directories")
    base=json.loads(args.baseline_manifest.read_text(encoding="utf-8"))
    assert sha(args.baseline_exe)==base["residual_executable_sha256"]
    inputs={str(i):sha(args.data/f"case-{i:03d}.b2a") for i in CASES}
    assert inputs==base["input_sha256"]
    files=sorted(list((args.source/"src").glob("*.cpp"))+list((args.source/"src").glob("*.h"))+list((args.source/"include").glob("*.h")))
    files += [args.source/p for p in ["CMakeLists.txt", "tests/MotherDistanceSelectedTests.cpp", "tests/Dataset500.cpp", "tests/PairedSurfaceBounds.cpp"]]
    hashes={p.relative_to(args.source).as_posix():sha(p) for p in files}
    # 所有实际连续认证实现文件必须与保存的基线一致。
    for key,value in hashes.items():
        if key.startswith(("src/Mother", "src/Paired", "include/Mother", "include/Paired", "src/Interval")):
            assert value==base["source_sha256"][key], "Validator changed: "+key
    manifest=dict(started_utc=datetime.now(timezone.utc).isoformat(), cases=CASES, runs=args.runs,
        source_sha256=hashes, input_sha256=inputs, executable_sha256=sha(args.exe),
        baseline_executable_sha256=sha(args.baseline_exe), baseline_manifest_sha256=sha(args.baseline_manifest),
        source_directory=str(args.source), commands=[], failures=[],
        runner_sha256={str(Path(__file__).name):sha(Path(__file__)), "RunMotherJointSelected.py":sha(Path(__file__).with_name("RunMotherJointSelected.py"))},
        contract="1e-5 positions and strict 0.1-degree U/V joint condition; all failed candidate validation costs included",
        order="serial case blocks alternating baseline/incremental first; process warmup then five cold/warm pairs")
    args.work.mkdir(parents=True); args.output.mkdir(parents=True)
    env=os.environ.copy(); env["PATH"]=str(args.sdk/"x64-win/bin")+os.pathsep+env["PATH"]
    rows=[]; all_attempts=[]; logs=[]
    for ix,case in enumerate(CASES):
        for method in (("baseline","incremental") if ix%2==0 else ("incremental","baseline")):
            dest=args.work/f"{method}-{case:03d}.csv"
            exe=args.baseline_exe if method=="baseline" else args.exe
            command=[str(exe),str(args.data),str(case),str(dest),str(args.runs)]
            manifest["commands"].append(dict(method=method,case=case,command=command))
            try:
                result=subprocess.run(command,cwd=args.sdk,env=env,capture_output=True,text=True,errors="replace",timeout=120)
                logs.append(f"{method} {case}, exit {result.returncode}\n{result.stdout}\n{result.stderr}")
                current=list(csv.DictReader(dest.open(encoding="utf-8",newline=""))) if dest.exists() else []
                rows.extend(dict(method=method,**r) for r in current)
                parts=[]
                if method=="incremental" and Path(str(dest)+".attempts.csv").exists():
                    parts=list(csv.DictReader(Path(str(dest)+".attempts.csv").open(encoding="utf-8",newline="")))
                    all_attempts.extend(parts)
                if result.returncode: raise RuntimeError(f"Exit {result.returncode}")
                if method=="baseline": audit(current,case,args.runs)
                else: check_incremental(current,parts,case,args.runs)
                print(f"{method} {case}: {len(current)} measurements passed",flush=True)
            except (subprocess.TimeoutExpired,RuntimeError,AssertionError,ValueError,KeyError) as exc:
                manifest["failures"].append(dict(method=method,case=case,reason=str(exc)))
                print(f"FAILED {method} {case}: {exc}",flush=True)
            save(args.output/"measurements.csv",rows); save(args.output/"attempts.csv",all_attempts)
            (args.output/"runs.log").write_text("\n".join(logs),encoding="utf-8")
            (args.output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    assert all(sha(args.source/k)==v for k,v in hashes.items()), "Build sources changed during measurement"
    if manifest["failures"]: return 2
    per_case=[]; by_type=[]; summary={}
    for method in ("baseline","incremental"):
        group=[r for r in rows if r["method"]==method]
        cold=[r for r in group if r["cache"]=="cold"]
        for case in CASES:
            for cache in ("cold","warm"):
                selected=[r for r in group if int(r["id"])==case and r["cache"]==cache]
                times=[float(r["algorithm_ms" if cache=="cold" else "validation_ms"]) for r in selected]
                per_case.append(dict(method=method,id=case,type=selected[0]["type"],cache=cache,runs=len(selected),**means(selected),
                    total_min_ms=min(times),total_max_ms=max(times),total_median_ms=statistics.median(times),runs_under_10ms=sum(t<10 for t in times)))
        for kind in ("Plane","Cylinder","Cone","Sphere","Torus"):
            selected=[r for r in cold if r["type"]==kind]
            by_type.append(dict(method=method,type=kind,cases=len({r["id"] for r in selected}),**means(selected)))
        summary[method]=dict(cases=len(CASES),measurements=len(group),passed=len(group),
            cold_case_means_under_10ms=sum(r["algorithm_ms"]<10 for r in per_case if r["method"]==method and r["cache"]=="cold"),
            cold_runs_under_10ms=sum(float(r["algorithm_ms"])<10 for r in cold),cold_max_ms=max(float(r["algorithm_ms"]) for r in cold),
            cold_mean_ms=statistics.mean(float(r["algorithm_ms"]) for r in cold),
            classification_mean_ms=statistics.mean(float(r["classification_ms"]) for r in cold),
            validation_mean_ms=statistics.mean(float(r["validation_ms"]) for r in cold))
    save(args.output/"per-case.csv",per_case); save(args.output/"modules-by-type.csv",by_type)
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    manifest["completed_utc"]=datetime.now(timezone.utc).isoformat()
    (args.output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2),flush=True)
    return 0

if __name__=="__main__": raise SystemExit(main())
