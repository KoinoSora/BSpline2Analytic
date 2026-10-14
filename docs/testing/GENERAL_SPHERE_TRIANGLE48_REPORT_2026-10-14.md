# 球面三角片补丁验证（2026-10-14）

本次仅执行 GeneralExtreme320 中的 `sphere_triangle` 组（ID 128–175），共 48 个固定输入：NURBS 和非有理样条各 24 个。所有输入与[修复前结果](GENERAL_EXTREME320_RESULTS_2026-10-14.csv)的 SHA-256 一致，未重新生成或替换数据。其他极端组、旧 500 例及一般 600 例均未运行。

## 补丁范围

`src/GeneralPairedBounds.cpp` 保留原有 `RegularQuad` 判断，在它不通过时，补充 `CollapsedSideTriangle`：恰有一条 UV 相邻边的两端完全重合、四角恰有三个不同顶点，且三角形行列式的向外舍入区间严格为正或负，才允许使用已有连续配对误差界。

此时双线性映射的 Jacobian 为非零三角形行列式乘 `s`、`1-s`、`t` 或 `1-t`，只在指定退化边为零，内部不折叠。近似重合端点不会被吸附成同一点；共线、整体退化和折叠域仍被拒绝。三维极点处的边退化不等同于 UV 边退化：UV 仍为正常四边形的极点片继续经过原判断。

新增判断只比较四个角点，使用固定大小栈数组，不增加采样、优化迭代或细分预算。分类、球到圆环回退、参数坐标系恢复及既有投影验证回退均未修改。距离容差仍为 `1e-5`；配对验证仍要求误差界的两倍严格小于容差，并保留原有细分限制与失败处理。

## 结果

| 子组 | 输入数 | 修复前通过 / 失败 / 超时 | 修复后通过 / 失败 / 超时 | 转换累计耗时 | 单例最大转换耗时 |
| --- | ---: | --- | --- | ---: | ---: |
| 普通 UV 三角片 | 16 | 4 / 4 / 8 | 16 / 0 / 0 | 274.708 ms | 43.310 ms |
| 极点三角片 | 16 | 16 / 0 / 0 | 16 / 0 / 0 | 315.902 ms | 31.221 ms |
| 细长三角片 | 16 | 1 / 4 / 11 | 16 / 0 / 0 | 202.140 ms | 28.502 ms |
| 合计 | 48 | 21 / 8 / 19 | 48 / 0 / 0 | 792.750 ms | 43.310 ms |

48 个输出均为球面，全部通过生产连续验收和独立输出认证。SDK PCurve 均为一次曲线，共享端点完全一致，PCurve 的三维求值一致性检查全部通过。每例只尝试一个母面候选；母面拟合优化迭代为零。参数坐标系恢复仍可能迭代，本次累计 37 次，不包含在母面优化迭代数内。

- 来源有限面片独立误差上界最大值：`9.779933310177196e-7`。
- 输出有限面片独立误差上界最大值：`1.5145969596752404e-6`。
- 生产几何及边界 Hausdorff 上界最大值：`4.997626012388016e-6`；其两倍仍小于 `1e-5`。
- 转换平均耗时 `16.516 ms`，中位数 `15.366 ms`。分类累计 `31.403 ms`，母面优化 `0 ms`，验证累计 `760.324 ms`；剩余时间为流程开销。
- 运行器墙钟总耗时 `14.996 s`，包含数据检查、逐例启动、加载、输入/输出独立认证及转换。它与转换器内部累计 `0.793 s` 的统计口径不同。

修复前每例外部超时上限同为 30 秒；19 个超时样本没有完整转换计时，因此不据此计算整组转换加速倍数。耗时为本机一次 Release 运行结果。本轮没有其他类型的回归测试结果，不据此宣称全部测试集无退化。

## 退化域定向检查

新增 `--sphere-triangle-checks`，仅加载本组 ID 128、144、160。三个正确域均通过生产配对验证；在这些球面输入上分别构造点域、共线三角形、折叠四边形、非相邻重复顶点、近退化但折叠域，共 15 次非法域检查，全部在进入几何细分前被拒绝（`piecesVisited == 0`）。这验证了新增入口没有将无效参数域直接放行。

## 构建与复现

Release 目标 `BSpline2AnalyticGeneralExtremeTests` 构建成功。测试可执行文件 SHA-256：`154f1b47883730e956e783b57f259f7a6d1b7d9ce91b785606d738c5d0c1abff`。运行器和定向检查均返回 0。

从仓库根目录执行以下命令，仅运行本组 48 例，`--last` 不包含端点：

```powershell
cmake --build build-general-extreme320 --config Release --target BSpline2AnalyticGeneralExtremeTests --parallel 2
python tests/RunGeneralExtreme.py --exe build-general-extreme320/Release/BSpline2AnalyticGeneralExtremeTests.exe --sdk ../../../SGGK --data tests/data/general-extreme320 --output build-general-extreme320/sphere-triangle48.csv --first 128 --last 176 --timeout 30
```

定向检查使用同一可执行文件：`--sphere-triangle-checks <general-extreme320 数据绝对路径>`，从 SDK 目录启动并将 `SGGK/x64-win/bin` 加入进程 PATH。Windows 环境含大小写重复的 PATH/Path 时，用 Python `subprocess.run(..., env=dict(os.environ))` 启动构建器。本次采用此方式。

完整逐例记录见 [球面三角片 48 例 CSV](GENERAL_SPHERE_TRIANGLE48_RESULTS_2026-10-14.csv)。历史 320 例结果保持不变。
