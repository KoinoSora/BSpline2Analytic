# 文件职责：配对验证器准备阶段的隔离基线对照。
# 校验共享几何依赖，比较冷缓存及复用缓存的成本，避免把其他算法变化归因于准备优化。

"""比较保存的配对验证器与当前冷缓存、复用缓存模式。
基线数值代码保持原样，只重命名私有命名空间/公共符号并补齐准备阶段子计时。
生成文件保留在 build/；需先构建当前 Release 库，不执行仓库检出或重置。"""
import argparse
import csv
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--cmake", default="cmake")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = args.build_dir.resolve()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    work = build / "paired-preparation-comparison"
    work.mkdir(parents=True, exist_ok=True)
    env = {k.upper(): v for k, v in os.environ.items()}
    git = ["git", "-c", f"safe.directory={root.as_posix()}"]

    def read_git(spec):
        return subprocess.check_output(git + ["show", spec], cwd=root, env=env).decode("utf-8")

    sha = subprocess.check_output(git + ["rev-parse", args.baseline + "^{commit}"],
                                  cwd=root, env=env, text=True).strip()
    # 共享依赖变化不能静默破坏新旧实现对照的有效性。
    for name in ("src/IntervalBounds.cpp", "src/IntervalBounds.h", "src/Numerics.cpp", "src/Numerics.h",
                 "src/SurfaceSampleContext.cpp", "src/SurfaceSampleContext.h"):
        if read_git(f"{sha}:{name}").replace("\r\n", "\n") != (root / name).read_text(encoding="utf-8"):
            raise ValueError(f"Baseline shared dependency differs: {name}")
    # 恢复接口现在可选保留轨迹；对隔离基线实际使用的几何辅助函数逐项校验，
    # 避免仅因无关的接口扩展拒绝对照。
    def function(text, marker):
        start = text.index(marker)
        brace = text.index("{", start)
        depth = 1
        end = brace + 1
        while depth:
            depth += (text[end] == "{") - (text[end] == "}")
            end += 1
        return re.sub(r"\s+", "", text[start:end])
    old_geometry = read_git(f"{sha}:src/ApproximateGeometry.cpp")
    new_geometry = (root / "src/ApproximateGeometry.cpp").read_text(encoding="utf-8")
    for marker in ("Frame3 FrameOf(", "std::array<double, 2> MotherParameters(",
                   "ClassificationSamples BuildClassificationSamples(",
                   "PatchRecoveryTrace TraceApproximateBoundaries(", "void CheckApproximateInput("):
        if function(old_geometry, marker) != function(new_geometry, marker):
            raise ValueError(f"Baseline geometry helper changed: {marker}")
    originals = {}
    files = ["PairedValidationInternal.h", "PairedRationalChart.cpp", "PairedConversionValidation.cpp"]
    if subprocess.run(git + ["cat-file", "-e", f"{sha}:src/PairedSourcePreparation.cpp"],
                      cwd=root, env=env, stderr=subprocess.DEVNULL).returncode == 0:
        files.append("PairedSourcePreparation.cpp")
    if subprocess.run(git + ["cat-file", "-e", f"{sha}:src/PairedPositionBounds.cpp"],
                      cwd=root, env=env, stderr=subprocess.DEVNULL).returncode == 0:
        files.append("PairedPositionBounds.cpp")
    for name in files:
        original = read_git(f"{sha}:src/{name}")
        originals[name] = hashlib.sha256(original.encode()).hexdigest()
        code = original.replace("detail::paired", "detail::paired_baseline")
        code = code.replace('"PairedValidationInternal.h"', '"PairedBaselineInternal.h"')
        code = code.replace("ValidatePairedConversion(", "ValidatePairedConversionBaseline(")
        code = code.replace("PairedValidationStatusName(", "PairedValidationStatusNameBaseline(")
        if name == "PairedConversionValidation.cpp":
            old = "        sourceCells = ExtractCells(source);"
            if old in code:
                assert code.count(old) == 1
                code = code.replace(old, "        auto extractionStart = Clock::now();\n" + old +
                                    "\n        out.sourceExtractionMs += Milliseconds(extractionStart);")
                old = "        for (auto &c : sourceCells)\n            c.p = LocalSource(c.p, patch.parameters);"
                assert code.count(old) == 1
                code = code.replace(old, "        auto localStart = Clock::now();\n" + old +
                                    "\n        out.localTransformMs += Milliseconds(localStart);")
            if "struct PairedValidationContext::Impl" in code:
                # 仅运行冷缓存基线；避免产生重复的公共上下文符号，
                # 同时保留原有冷启动辅助函数和数值实现。
                code = code[:code.index("struct PairedValidationContext::Impl")] + '''
PairedValidationReport ValidatePairedConversionBaseline(const sggk::BSplineSurface &s,
 const AnalyticPatch &p, const PairedValidationOptions &o) { return ValidatePaired(s,p,o,nullptr); }
std::shared_ptr<void> MakePairedBaselineCache()
{ return std::make_shared<detail::paired_baseline::SourceCache>(); }
PairedValidationReport ValidatePairedConversionBaselineCached(const sggk::BSplineSurface &s,
 const AnalyticPatch &p, const PairedValidationOptions &o, void *cache)
{ return ValidatePaired(s,p,o,static_cast<detail::paired_baseline::SourceCache *>(cache)); }
} // namespace bspline2analytic
'''
        if name == "PairedConversionValidation.cpp" and "MakePairedBaselineCache" not in code:
            raise ValueError("This comparison requires a baseline with SourceCache support")
        target = "PairedBaselineInternal.h" if name.endswith(".h") else name
        (work / target).write_text(code, encoding="utf-8")
    cmake = r'''cmake_minimum_required(VERSION 3.22)
project(PairedPreparationComparison LANGUAGES CXX)
add_executable(PairedPreparationComparison
 "${REPO}/tests/PairedPreparationBenchmark.cpp" "${REPO}/tests/Dataset500.cpp"
 "${REPO}/tests/PairedSurfaceBounds.cpp" PairedRationalChart.cpp PairedConversionValidation.cpp)
target_compile_features(PairedPreparationComparison PRIVATE cxx_std_17)
target_compile_options(PairedPreparationComparison PRIVATE /W4 /utf-8 /fp:strict /permissive-)
target_compile_definitions(PairedPreparationComparison PRIVATE NOMINMAX)
target_include_directories(PairedPreparationComparison PRIVATE "${REPO}/include" "${REPO}/src" "${REPO}/tests")
target_include_directories(PairedPreparationComparison SYSTEM PRIVATE "${SDK}/include")
target_link_directories(PairedPreparationComparison PRIVATE "${SDK}/x64-win/lib")
target_link_libraries(PairedPreparationComparison PRIVATE "${LIBRARY}" SGGK_Foundation SGGK_Math SGGK_GeomBase SGGK_Geometry SGGK_BSplineAlgo)
'''
    if "PairedSourcePreparation.cpp" in files:
        cmake += '\ntarget_sources(PairedPreparationComparison PRIVATE PairedSourcePreparation.cpp)\n'
    if "PairedPositionBounds.cpp" in files:
        cmake += '\ntarget_sources(PairedPreparationComparison PRIVATE PairedPositionBounds.cpp)\n'
    (work / "CMakeLists.txt").write_text(cmake, encoding="utf-8")
    # 读取本次库构建实际选用的 SDK 路径，不假设固定安装位置。
    cache = (build / "CMakeCache.txt").read_text(encoding="utf-8")
    sdk = Path(next(line.split("=", 1)[1] for line in cache.splitlines()
                    if line.startswith("SGGK_SDK_DIR:PATH="))).resolve()
    library = build / "Release/BSpline2Analytic.lib"
    if not library.is_file():
        raise FileNotFoundError("Build the current Release library first")
    env["PATH"] = str(sdk / "x64-win/bin") + os.pathsep + env.get("PATH", "")
    commands = [
        ("configure", [args.cmake, "-S", str(work), "-B", str(work / "out"), "-G",
                       "Visual Studio 17 2022", "-A", "x64", f"-DREPO={root.as_posix()}",
                       f"-DSDK={sdk.as_posix()}", f"-DLIBRARY={library.as_posix()}"], root),
        ("build", [args.cmake, "--build", str(work / "out"), "--config", "Release", "--parallel", "2"], root),
        ("run", [str(work / "out/Release/PairedPreparationComparison.exe"),
                 str(root / "tests/data/dataset500"), str(output), str(args.runs)], sdk),
    ]
    for name, command, cwd in commands:
        log = output.with_suffix(f".{name}.log")
        with log.open("wb") as stream:
            result = subprocess.run(command, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT)
        print(f"{name}: exit {result.returncode}; {log}", flush=True)
        if result.returncode:
            raise RuntimeError(log.read_text(encoding="utf-8", errors="replace")[-4000:])
    with output.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 40 * args.runs and len({r["id"] for r in rows}) == 40
    mean = lambda key: sum(float(r[key]) for r in rows) / len(rows)
    summary = {"baseline": sha, "cases": 40, "repetitions": args.runs,
               "baseline_original_sha256": originals,
               "current_library_sha256": hashlib.sha256(library.read_bytes()).hexdigest(),
               "configuration": "Release /fp:strict", "cpu_affinity_fixed": False,
               "exclusive_machine": False, "cache_fill_mean_ms": mean("cache_fill_ms"),
               "baseline_cache_fill_mean_ms": mean("baseline_cache_fill_ms"),
               "means": {mode: {field: mean(f"{mode}_{field}") for field in
                                ("wall_ms", "prepare_ms", "extract_ms", "local_ms", "cache_ms",
                                 "surface_ms", "boundary_ms", "compose_ms", "product_ms",
                                 "difference_ms", "fast_bound_ms", "tight_bound_ms",
                                 "source_center_ms", "source_boundary_ms", "target_transform_ms",
                                 "target_boundary_ms", "edge_compose_ms", "edge_algebra_ms",
                                 "edge_bounds_ms", "edge_subdivision_ms")}
                         for mode in ("baseline", "new", "cached", "baseline_cached")}}
    summary["speedup"] = {"cold": mean("baseline_wall_ms") / mean("new_wall_ms"),
                          "cached": mean("baseline_cached_wall_ms") / mean("cached_wall_ms")}
    summary["position_speedup"] = mean("baseline_surface_ms") / mean("new_surface_ms")
    summary["new_fast_bound_passes"] = sum(int(r["new_fast_passes"]) for r in rows)
    summary["new_tight_bound_checks"] = sum(int(r["new_tight_checks"]) for r in rows)
    output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
