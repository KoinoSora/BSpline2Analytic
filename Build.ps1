# 文件职责：模块构建与受支持回归的 PowerShell 入口。
# 校验互斥开关和废弃路径，配置本地 SDK、构建目标并检查退出码；许可证和 SDK 仅从本机读取。

param(
    [string]$SdkDirectory = (Join-Path $PSScriptRoot '..\..\..\SGGK'),
    [ValidateSet('Debug','Release','RelWithDebInfo')][string]$Configuration = 'Release',
    [switch]$Benchmark,
    [switch]$ApproximateBenchmark,
    [switch]$PlateauOnly,
    [switch]$JacobianOnly,
    [switch]$Dataset500Only,
    [switch]$GeneralOnly,
    [switch]$GeneralPositiveOnly,
    [switch]$ExtremeBoundaryOnly,
    [switch]$ConeOnly,
    [switch]$SpherePolynomialOnly,
    [ValidateRange(1,64)][int]$BuildJobs = 1
)
$ErrorActionPreference = 'Stop'
if ($GeneralOnly -or $GeneralPositiveOnly) {
    throw 'The analytic-UV fallback datasets are retired. Use RunGeneralRationalBoundaryFull.py with general-rational-v2.'
}
if ($ApproximateBenchmark -or $PlateauOnly -or $JacobianOnly -or $ConeOnly -or $ExtremeBoundaryOnly) {
    throw 'This test path invokes retired ValidateConversion. Its code is retained, but testing is disabled pending migration.'
}
$focusedCount = [int]$PlateauOnly.IsPresent + [int]$JacobianOnly.IsPresent + [int]$Dataset500Only.IsPresent + [int]$GeneralOnly.IsPresent + [int]$GeneralPositiveOnly.IsPresent + [int]$ExtremeBoundaryOnly.IsPresent + [int]$ConeOnly.IsPresent + [int]$SpherePolynomialOnly.IsPresent
if ($focusedCount -gt 1 -or ($focusedCount -gt 0 -and ($Benchmark -or $ApproximateBenchmark))) {
    throw 'Choose only one focused test switch and do not combine it with broader benchmarks'
}
$sdk = (Resolve-Path -LiteralPath $SdkDirectory).Path
$build = Join-Path $PSScriptRoot 'build'
cmake -S $PSScriptRoot -B $build -G 'Visual Studio 17 2022' -A x64 "-DSGGK_SDK_DIR=$sdk"
if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed' }
if ($ExtremeBoundaryOnly) {
    cmake --build $build --config $Configuration --target BSpline2AnalyticExtremeTests --parallel $BuildJobs
} elseif ($GeneralOnly -or $GeneralPositiveOnly) {
    cmake --build $build --config $Configuration --target BSpline2AnalyticGeneralDatasetTests --parallel $BuildJobs
} elseif ($ConeOnly) {
    cmake --build $build --config $Configuration --target BSpline2AnalyticDataset500Tests BSpline2AnalyticConeTests --parallel $BuildJobs
} elseif ($Dataset500Only -or $SpherePolynomialOnly) {
    cmake --build $build --config $Configuration --target BSpline2AnalyticPairedValidationTests --parallel $BuildJobs
} elseif ($JacobianOnly) {
    cmake --build $build --config $Configuration --target BSpline2AnalyticJacobianTests --parallel $BuildJobs
} elseif ($PlateauOnly) {
    cmake --build $build --config $Configuration --target BSpline2AnalyticPlateauTests --parallel $BuildJobs
} else {
    cmake --build $build --config $Configuration --parallel $BuildJobs
}
if ($LASTEXITCODE -ne 0) { throw 'C++ build failed' }
if ($ExtremeBoundaryOnly) {
    ctest --test-dir $build -C $Configuration -R '^BSpline2Analytic\.Extreme(200|Modules|Data)$' --output-on-failure
} elseif ($GeneralPositiveOnly) {
    ctest --test-dir $build -C $Configuration -R '^BSpline2Analytic\.GeneralPositive500$' --output-on-failure
} elseif ($GeneralOnly) {
    ctest --test-dir $build -C $Configuration -R '^BSpline2Analytic\.General(600|Extra60)$' --output-on-failure
} elseif ($ConeOnly) {
    ctest --test-dir $build -C $Configuration -R '^BSpline2Analytic\.(ConeDataset|Cone)$' --output-on-failure
} elseif ($SpherePolynomialOnly) {
    ctest --test-dir $build -C $Configuration -R '^BSpline2Analytic\.PairedSphere$' --output-on-failure
} elseif ($Dataset500Only) {
    ctest --test-dir $build -C $Configuration -R '^BSpline2Analytic\.PairedDataset500$' --output-on-failure
} elseif ($JacobianOnly) {
    ctest --test-dir $build -C $Configuration -R '^BSpline2Analytic\.Jacobian$' --output-on-failure
} elseif ($PlateauOnly) {
    ctest --test-dir $build -C $Configuration -R '^BSpline2Analytic\.Plateau$' --output-on-failure
} else {
    ctest --test-dir $build -C $Configuration --output-on-failure
}
if ($LASTEXITCODE -ne 0) { throw 'Regression tests failed' }
if ($Benchmark -or $ApproximateBenchmark) {
    $suffix = if ($Configuration -eq 'Debug') {'d'} elseif ($Configuration -eq 'RelWithDebInfo') {'w'} else {''}
    $previousPath = $env:PATH
    try {
        $env:PATH = (Join-Path $sdk ('x64-win\bin' + $suffix)) + ';' + $env:PATH
        Push-Location $sdk
        try {
            if ($Benchmark) {
                & (Join-Path $build "$Configuration\BSpline2AnalyticTests.exe") --benchmark (Join-Path $build "benchmark-$Configuration.csv")
                if ($LASTEXITCODE -ne 0) { throw 'Exact benchmark failed' }
            }
            if ($ApproximateBenchmark) {
                & (Join-Path $build "$Configuration\BSpline2AnalyticApproximateBenchmark.exe") (Join-Path $build "approximate-benchmark-$Configuration.csv")
                if ($LASTEXITCODE -ne 0) { throw 'Approximate benchmark failed' }
            }
        }
        finally { Pop-Location }
        if ($LASTEXITCODE -ne 0) { throw 'Benchmark or batch consistency check failed' }
    } finally { $env:PATH = $previousPath }
}
