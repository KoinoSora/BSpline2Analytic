# 文件职责：方案设计文件的 PDF 构建入口。
# 编译当前目录中的 LaTeX 源文件，不参与算法运行或性能测量。

param([string]$XeLaTeX = 'xelatex')
$ErrorActionPreference = 'Stop'
$documentName = '面向线性PCurve的一类通用样条面到二次曲面智能转化方案—方案设计文件'
$buildDirectory = Join-Path $PSScriptRoot 'build'
New-Item -ItemType Directory -Path $buildDirectory -Force | Out-Null
Push-Location $PSScriptRoot
try {
    for ($pass = 1; $pass -le 2; $pass++) {
        & $XeLaTeX '-interaction=nonstopmode' '-halt-on-error' "-output-directory=$buildDirectory" "$documentName.tex" > (Join-Path $buildDirectory "compile-$pass.txt")
        if ($LASTEXITCODE -ne 0) { throw "XeLaTeX pass $pass failed; see build/compile-$pass.txt" }
    }
    Copy-Item -LiteralPath (Join-Path $buildDirectory "$documentName.pdf") -Destination (Join-Path $PSScriptRoot "$documentName.pdf")
    Write-Output (Join-Path $PSScriptRoot "$documentName.pdf")
} finally { Pop-Location }
