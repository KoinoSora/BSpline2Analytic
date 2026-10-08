param(
    [string]$SdkDirectory = (Join-Path $PSScriptRoot '..\..\..\SGGK'),
    [ValidateSet('Debug','Release','RelWithDebInfo')][string]$Configuration = 'Release',
    [switch]$Benchmark
)
$ErrorActionPreference = 'Stop'
$sdk = (Resolve-Path -LiteralPath $SdkDirectory).Path
$build = Join-Path $PSScriptRoot 'build'
cmake -S $PSScriptRoot -B $build -G 'Visual Studio 17 2022' -A x64 "-DSGGK_SDK_DIR=$sdk"
if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed' }
cmake --build $build --config $Configuration --parallel
if ($LASTEXITCODE -ne 0) { throw 'C++ build failed' }
ctest --test-dir $build -C $Configuration --output-on-failure
if ($LASTEXITCODE -ne 0) { throw 'Regression tests failed' }
if ($Benchmark) {
    $suffix = if ($Configuration -eq 'Debug') {'d'} elseif ($Configuration -eq 'RelWithDebInfo') {'w'} else {''}
    $previousPath = $env:PATH
    try {
        $env:PATH = (Join-Path $sdk ('x64-win\bin' + $suffix)) + ';' + $env:PATH
        Push-Location $sdk
        try { & (Join-Path $build "$Configuration\BSpline2AnalyticTests.exe") --benchmark (Join-Path $build "benchmark-$Configuration.csv") }
        finally { Pop-Location }
        if ($LASTEXITCODE -ne 0) { throw 'Benchmark or batch consistency check failed' }
    } finally { $env:PATH = $previousPath }
}
