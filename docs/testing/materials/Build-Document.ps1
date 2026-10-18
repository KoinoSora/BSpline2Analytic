# 文件职责：测试材料说明的 PDF 构建入口。
# 文档编译与算法测试分离，输出结果以编译器退出状态为准。

param([string]$XeLaTeX = 'xelatex')
$ErrorActionPreference = 'Stop'
$documentName = '测试材料说明'
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
