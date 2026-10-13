param([string]$OutputDirectory = (Join-Path $PSScriptRoot 'build'))
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force -Path $OutputDirectory | Out-Null
$destination = (Resolve-Path -LiteralPath $OutputDirectory).Path
Push-Location $PSScriptRoot
try {
    for ($pass = 0; $pass -lt 2; $pass++) {
        & xelatex -interaction=nonstopmode -halt-on-error -file-line-error "-output-directory=$destination" main.tex
        if ($LASTEXITCODE -ne 0) { throw 'XeLaTeX document build failed' }
    }
    Copy-Item -LiteralPath (Join-Path $destination 'main.pdf') -Destination (Join-Path $PSScriptRoot 'BSpline2Analytic-algorithm.pdf')
} finally { Pop-Location }
