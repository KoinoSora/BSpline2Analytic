# 算法说明书

[阅读 PDF](BSpline2Analytic-algorithm.pdf) · [排版源文件](main.tex) · [正文源文件](content.tex)

正文为 A4 十页，含十幅矢量流程图、六张表；按总体架构、子模块数学方法、真实案例、工程化与缺口核对组织。本文反映已有实现，不是新的算法修改方案。PDF 是正式交付件，本文件只提供编译和来源信息。

## 实现与数据基线

- 编写日期：2026-10-13。
- 阅读的仓库版本：`3c772a716ffe9adc650a61b14c67d26780e6b303`。
- 生产转换算法最近变更：`c62044addd325e5cdcd5ffecedfb0682a83319c6`。
- 环境：SGK 1.7.11、Windows x64、MSVC 19.42、Release；正式批次每例运行一次，每类型/实际分支预热一次。
- 当前批次汇总见 [500 例重测报告](../../DATASET500_RETEST_2026-10-13.md)，保存的输入及生成流程见 [Dataset500](../../DATASET500.md)。

第九页引用以下已有本地记录，路径均相对于模块根目录。`build/` 按项目规则不提交；哈希用于区分当前文件与后续重跑结果，报告中的汇总已随仓库保存。

| 记录 | 用途 |
| --- | --- |
| `build/dataset500-20261013-current-Release.csv` | 当前 500 例正式结果；案例 100、249、385、317 的时间与误差界。 |
| `build/repeat385-20261013-release.txt` | 同实现下案例 385 的重复诊断；候选、采样、验收片计数。 |
| `build/verification-cone-focused-release.txt` | 当前生产算法圆锥专项；扰动例两种 Jacobian 模式均为四次通用迭代。 |
| `build/plateau-Release.csv` | 历史平面平台专项；只用于说明标量直接求解，不与当前批次混算时间。 |

上述文件的 SHA-256 依次为：

```text
EDCA4AAF82F536614EF1B8923EB07C8A5821B16E6AC294F3EEF8A0E4B20BCFE6
3BB536BF9064374280CB8E4E1B2B1701EF3136277AC1A203B0C9110E1A1E9C01
8305A3BF1F06AFE5EB97926824F99E7DD64AC14159C8AB2DCDEFD616DAB1DF0F
D7309E613B996DCA300EF325FD264D20A504C3BD44F104AA82116C000CBC3331
```

本轮未重新执行算法测试；复用既有结果并核对源码。全量普通案例均零优化，故正文另以扰动专项解释参数修正，未编造日志没有记录的逐次轨迹。

## 编译与版式核验

需要包含 `ctex`、Fandol 字体、TikZ、`xurl` 等包的 XeLaTeX。已使用 TeX Live 2022 验证；将 `xelatex` 加入 PATH 后，在 PowerShell 运行：

```powershell
.\Build-Document.ps1
```

脚本编译两遍以更新总页数和书签，中间产物放在本目录 `build/`，最终更新 `BSpline2Analytic-algorithm.pdf`。仅 PDF、TeX 源文件、构建脚本及本说明进入版本控制。

本次已用 `pdfinfo` 核验十页 A4，渲染并逐页检查全部页面；编译日志没有溢出版心、缺失字符或未定义引用。后续修改必须重新编译和渲染；若不能在保持可读性的前提下控制十页，依算法文档 Skill 先说明所需增页数量并取得用户同意。
