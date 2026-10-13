"""Compare a Git baseline with the built Release library on the 500 positives."""
import argparse
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, help="Commit with the original two-argument classifier")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = root / "build"
    library = build / "Release/BSpline2Analytic.lib"
    if not library.is_file():
        raise SystemExit("Build the current Release library first")
    cache = (build / "CMakeCache.txt").read_text(encoding="utf-8")
    sdk = Path(next(line.split("=", 1)[1] for line in cache.splitlines()
                    if line.startswith("SGGK_SDK_DIR:PATH="))).resolve()
    work = build / "general-classification-benchmark"
    work.mkdir(parents=True, exist_ok=True)
    source = subprocess.run(["git", "show", f"{args.baseline}:src/GeneralClassification.cpp"],
                            cwd=root, check=True, capture_output=True, encoding="utf-8").stdout
    if source.count("ClassifyGeneral(") != 1:
        raise SystemExit("Baseline classifier signature cannot be identified uniquely")
    (work / "Baseline.cpp").write_text(source.replace("ClassifyGeneral(", "ClassifyGeneralBeforeDedup("),
                                      encoding="utf-8")
    (work / "CMakeLists.txt").write_text(f'''cmake_minimum_required(VERSION 3.22)
project(GeneralClassificationBenchmark LANGUAGES CXX)
add_executable(GeneralClassificationBenchmark Baseline.cpp
    "{root.as_posix()}/tests/GeneralClassificationBenchmark.cpp"
    "{root.as_posix()}/tests/GeneralDataset.cpp"
    "{root.as_posix()}/tests/GeneralDatasetCertificate.cpp")
target_compile_features(GeneralClassificationBenchmark PRIVATE cxx_std_17)
target_compile_options(GeneralClassificationBenchmark PRIVATE /W4 /utf-8 /permissive-)
target_compile_definitions(GeneralClassificationBenchmark PRIVATE NOMINMAX)
target_include_directories(GeneralClassificationBenchmark PRIVATE
    "{root.as_posix()}/include" "{root.as_posix()}/src" "{root.as_posix()}/tests")
target_include_directories(GeneralClassificationBenchmark SYSTEM PRIVATE "{sdk.as_posix()}/include")
set_source_files_properties("{root.as_posix()}/tests/GeneralDatasetCertificate.cpp" PROPERTIES COMPILE_OPTIONS "/fp:strict")
target_link_directories(GeneralClassificationBenchmark PRIVATE "{sdk.as_posix()}/x64-win/lib")
target_link_libraries(GeneralClassificationBenchmark PRIVATE "{library.as_posix()}"
    SGGK_Foundation SGGK_Math SGGK_GeomBase SGGK_Geometry SGGK_BSplineAlgo)
''', encoding="utf-8")
    output = work / "out"
    subprocess.run(["cmake", "-S", str(work), "-B", str(output), "-G", "Visual Studio 17 2022", "-A", "x64"], check=True)
    subprocess.run(["cmake", "--build", str(output), "--config", "Release", "--parallel", "2"], check=True)
    env = os.environ.copy()
    env["PATH"] = str(sdk / "x64-win/bin") + os.pathsep + env.get("PATH", "")
    subprocess.run([str(output / "Release/GeneralClassificationBenchmark.exe"),
                    str(root / "tests/data/dataset500"), str(build / "general-classification-paired.csv")],
                   cwd=sdk, env=env, check=True)


if __name__ == "__main__":
    main()
