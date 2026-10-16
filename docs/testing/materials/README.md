# 测试材料说明

[PDF](测试材料说明.pdf) · [LaTeX 源文件](测试材料说明.tex)

依据 `70e03d1` 的实现和已提交批次，覆盖全部固定输入构造族、程序化单测、当前与停用入口、C++ API、单例／选定多例／整批命令，以及正确性和整体耗时判定。本文编制没有重跑几何测试。

## 附件

- [逐例索引](case-inventory.csv)：2,866 个固定输入（数据目录 2,865 个，展示圆环三角片 1 个），附路径、状态、生成器、原始元数据及 SHA-256。按文件身份计数，不声称独立几何形状数。
- [测试源文件索引](test-source-inventory.csv)：83 个测试源文件、构造器和脚本的哈希。
- [索引汇总](inventory-summary.json)：分组数量和附件哈希。
- [验证记录](verification.json)：本次文档及静态检查结果。

## 更新材料

从仓库根目录执行：

```powershell
python docs/testing/materials/Build-Inventory.py
powershell -NoProfile -ExecutionPolicy Bypass -File docs/testing/materials/Build-Document.ps1
```

索引脚本只读 Git 跟踪的数据及测试实现，不执行任何几何算法。增加输入或改变配额后，需要审查脚本的数量断言和正文；新增数据须先进入 Git 跟踪范围才能被收录。V2 文件哈希须符合其声明清单。

PDF 使用 Windows 中文字体及 XeLaTeX，全部流程图内嵌。编译中间文件位于忽略的 `build/` 目录。测试命令位于第八章，按照当前入口白名单执行；历史数据与旧验证器仍遵循根目录 `AGENTS.md` 的禁测要求。
