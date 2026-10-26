---
name: closed-loop-validation
description: 为 BSpline2Analytic 面向 main 的 PR 整理闭环验证、证据和管理员审批；实验阶段只运行无 SDK 的合成协议演练，不运行几何或性能测试。
---

# 闭环提交验证

本项目的门禁规则见 [需求与运行规范](../../../docs/requirements/CLOSED_LOOP_VALIDATION.md)，实现和证据协议见 [算法路径](../../../docs/algorithm/CLOSED_LOOP_VALIDATION.md)。

## 当前模式

`tests/closed_loop/policy.json` 的 `real_tests_enabled=false`，基线注册表尚未获批。用户于 2026-10-26 要求只检查流程，暂不启动 AGENT 的真实测试。不得调用 SGK、构建几何测试、运行原数据集、基准测试或采集正式基线。可以执行 Python 协议单测和合成演练；不调用其他模型模拟算法结果。

`SIMULATION_PASS`、`closed-loop/protocol-smoke`、`closed-loop/dry-run` 仅表示协议连通。不得写入正式 `closed-loop/gate` 的 success，不得把合成耗时写成实测结果，不得自动合并。

## 工作流程

1. 查看 Git 状态，保留已有用户修改；面向 main 的修改使用工作分支和 PR。纯分析不产生提交。
2. 从管理员批准的控制版本取得策略、测试清单和基线，不能采用待测 PR 自行放宽的门槛。实验期使用当前框架演练，明确不具备正式信任来源。
3. 运行 `python -m unittest discover -s tests/closed_loop -p 'test_*.py' -v`。这只测试协议，不调用产品算法。
4. 获取 PR 的 head、当前 main 和预合并树。按规范执行 `run.py smoke`，核验包内哈希。数值全部合成；真实管理员批准在演练中只作状态机模拟。
5. 在用户授权的仓库交付任务中，可使用 `github_bridge.py publish-dry-run` 上传固定格式的合成包，并设置独立 dry-run 状态。普通分析不授权远端写入。GitHub 连接器能力不足时使用仓库现有认证；不得记录凭据。
6. 发布后读取远端证据和状态，报告 PR、提交、执行链接及剩余激活条件。遇到提交或 main 变化，旧结果失效，重新冻结上下文；网络错误最多重试一次，其他错误先查原因。
7. 管理员审阅后才能合并。不得以作者自身批准、旧提交的批准、普通用户批准或模拟批准替代真实指定管理员批准。

## 正式启用边界

正式基线、可信隔离执行器、数据适配器和远端强制规则均须单独批准。用户授权正式启动后，先完成规范中的激活检查；不能只把 JSON 开关改成 true。正式入口必须覆盖完整 API 与合法数据集，保留原集 5%、新增案例 <10 ms、对应几何合同和独立审计。任何缺测、超时、证据不完整或统计不确定均不能放行。

AGENT 负责解释失败、生成整改清单、核对算法和影响文档；确定性脚本负责数值判定。不得删除失败案例、修改预期或门槛以修复失败。整改产生新提交后重新验证，最多一次预先规定的环境复测，不得重复挑选最好批次。
