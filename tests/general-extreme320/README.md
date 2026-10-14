# GeneralExtreme320 独立测试入口

从仓库根目录配置 `cmake -S tests/general-extreme320 -B build-general-extreme320`，构建 `BSpline2AnalyticGeneralExtremeTests`。该入口关闭其他测试目标，只注册 `BSpline2Analytic.GeneralExtreme320`。

输入位于 `../data/general-extreme320/`（相对于本目录），共 320 个固定样条；NURBS、非有理各 160 个。执行、数据资格及失败分类详见 [GeneralExtreme320 构造说明](../../docs/testing/GENERAL_EXTREME320.md)。

测试只调用一般兜底 `ConvertGeneral`，不通过统一入口调用前两个转换模块。
