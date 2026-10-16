# 测试报告（实验报告）

[PDF](测试报告.pdf) · [LaTeX 源文件](测试报告.tex)

依据原始《中船海舟样条转二次曲面》PDF 的命题目标、测试材料及实验报告要求撰写。包含需求追踪、版本与环境、设计可行性、接口覆盖、正确性与整体时间、误差分析、退化案例、稳定性边界、不足及优化空间。

主证据为 1,160 例串联批次，另分别引用 985 例一般验证器批次、900 例专用压力批次和 3 个退化展示例。不同版本、入口、数据成员和计时合同不合并为同一回归。

**版本核验：**主批次记录的源码指纹有 8 个文件与文档整理基线 `d527849` 不同；1,160 个输入哈希均一致。退化展示批次清单内源码与当前一致。历史主批次不能作为当前版本的全量验收结论。此次只复算已存记录，未重新执行几何算法。

## 附件与复核

- [report-evidence.json](report-evidence.json)：原始数据复算、测试源码差异、输入校验和证据 SHA-256。
- [primary-not-accepted.csv](primary-not-accepted.csv)：主批次 187 个未接受案例，其中 183 个预算失败，4 个非超时拒绝。
- [verification.json](verification.json)：文档编译、页数、布局和静态验证记录。
- [测试材料说明](../materials/README.md)：详细构造、单例／多例调用及 2,866 个固定输入索引。

从仓库根目录重建材料：

```powershell
python docs/testing/report/Build-Evidence.py
powershell -NoProfile -ExecutionPolicy Bypass -File docs/testing/report/Build-Document.ps1
```

证据脚本只读取原始记录及源码，不启动 SDK 或转换器。若当前源码继续变化，应同步审查版本差异叙述；不会自动把旧批次改写成新测试结果。
