param(
    [string]$SdkDirectory = (Join-Path $PSScriptRoot '..\..\..\SGGK'),
    [ValidateSet('Debug','Release','RelWithDebInfo')][string]$Configuration = 'Release',
    [switch]$Benchmark,
    [switch]$ApproximateBenchmark,
    [ValidateRange(1,64)][int]$BuildJobs = 1
)
$ErrorActionPreference = 'Stop'
$sdk = (Resolve-Path -LiteralPath $SdkDirectory).Path
$build = Join-Path $PSScriptRoot 'build'
cmake -S $PSScriptRoot -B $build -G 'Visual Studio 17 2022' -A x64 "-DSGGK_SDK_DIR=$sdk"
if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed' }
cmake --build $build --config $Configuration --parallel $BuildJobs
if ($LASTEXITCODE -ne 0) { throw 'C++ build failed' }
ctest --test-dir $build -C $Configuration --output-on-failure
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
