# 闭环提交协议：实验验证记录

本轮仅验证协议流程，不启动模型 AGENT，不运行 C++、SGK、几何数据集或真实性能测试。不存在本轮算法准确率、实际 5% 回退或实际 10 ms 达标结论。

## 本地验证

Windows / Python 3.11。`python -m unittest discover -s tests/closed_loop -p 'test_*.py' -v` 共 **34 项协议测试通过**。

覆盖合成成功流程；5% 精确边界与确定退化；固定基线防累计退化；单例退化不能被总体抵消；新增和既有 10 ms 边界；漏例、重例、空集、漏轮次、无效耗时；负例误接受、超时不算拒绝、拒绝证据；几何界超限与缺少审计；合同错配、退休入口与数据；历史未认证恢复仍受性能约束；主分支/提交/合并树失效；规则和环境变化；缺少文档；统计不确定；指定管理员、作者自批、撤销和要求修改；包篡改；正式模式保持关闭。

故障注入使用两个虚构案例及人工构造的数字。20 轮、4.00/4.04/6.00 ms 等均为测试常量，不能引用为产品性能。诊断也不能把 NOT_RUN、缺失状态或崩溃改标为预期成功；明确的超时诊断与普通负例超时分开判断。

## 远端演练路径

工作分支 push 和面向 main 的 PR 触发 `.github/workflows/closed-loop-protocol.yml`。GitHub-hosted Ubuntu 执行同一协议检查，生成合成报告、JSON、CSV 和文件哈希并上传运行附件。提交、base 与预合并树记录在证据内。

专用发布命令重新读取实时 PR，将包只追加到 `validation-evidence`，校验远端 Git blob 与分支引用，然后写 `closed-loop/dry-run` 状态。实际远端运行结果、日志和永久报告应通过对应 PR 的 Checks 与该状态链接查阅；本文件不预填尚未发生的远端成功结果。

`review-status` 只读取真实 review 并显示等待状态；单测中的管理员批准属于模拟事件，不是 GitHub 审批。本轮不执行 merge。

## 保留状态

`real_tests_enabled=false`；基线为空且未批准；没有部署正式采集器或修改 main 分支保护。`closed-loop/protocol-smoke` 和 `closed-loop/dry-run` 不等于正式 `closed-loop/gate`，不能作为真实合并授权。激活条件见[运行规范](../requirements/CLOSED_LOOP_VALIDATION.md)，实现路径见[算法说明](../algorithm/CLOSED_LOOP_VALIDATION.md)。

不涉及 `src/`、`include/`、CMake、原几何数据或算法容差修改。外部已有工作不包含在此提交中。
