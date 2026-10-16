# rational-domain-v2

985 个固定输入，正例边界定义于母面 NURBS 参数域。两方向次数均为 **3–5**，控制顶点总数 **≤200**，算法几何容差为 **1e-5**。构造、合同、来源证书与次数统计见 [数据说明](../../../docs/testing/GENERAL_RATIONAL_DATASET_V2.md)。

- `dataset.json`：文件 SHA256、案例预期、来源认证程序及 CSV 哈希。
- `qualification.json`、`migration.csv`：配额、次数变化和历史规格来源。
- 各子集 `manifest.csv`：固定输入及来源构造界；`audit-0-N.csv`：最终复核界与采样几何尺寸。
- `generation-logs/`：数据生成和来源认证日志，未调用兜底模块。

2026-10-16：生成与来源认证已完成；兜底单测、分类、转换、性能测试均未运行。20 个非线性边界专项是单独的探针，不属于直线边界正例。
