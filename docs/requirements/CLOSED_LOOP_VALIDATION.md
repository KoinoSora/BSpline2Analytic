# 闭环提交验证：实验运行规范

日期：2026-10-26。当前为 **协议演练阶段**，真实测试和正式合并门禁未激活。

## 1. 当前交付范围

已提供 AGENT Skill、固定规则判定器、逐例证据协议、故障注入单测、GitHub Actions 协议演练、PR 模板、CODEOWNERS，以及 GitHub 证据发布桥接器。全过程无需 SGK、编译器、模型 API 或第三方 Python 包。

本轮不执行几何数据集、性能基准、C++ 构建或正式基线采集，不启动模型 AGENT。`baseline-registry.json` 故意为空且未批准；不把历史不完整报告提升为新基线。GitHub 工作流只执行合成协议，不能因其变绿就宣称几何认证或 10 ms 达标。

正式测试执行器的隔离部署、现有测试输出的合同适配及正式基线采集属于后续激活步骤，当前不能直接用于真实合并授权。这个限制由关闭开关、未批准注册表、可信执行器条件和独立状态名称共同表达，而非在报告中隐瞒。

## 2. 规则

| 对象 | 规则 |
| --- | --- |
| 原成功案例 | 每轮均成功，原合同检查和独立审计通过；不得以其他案例恢复抵消回退 |
| 原负例 | 正确拒绝；TIMEOUT、CRASH、NOT_RUN、INVALID 不计正确拒绝 |
| 历史未认证 | 注册表包含明确 debt_issue；保留 UNRESOLVED 或经审计恢复 ACCEPT，其他状态阻断 |
| 旧集耗时 | 对冻结基线和当前 main 分别比较，每例及分组累计中位耗时均不得增加超过 5% |
| 原低于 10 ms 案例 | 冻结基线的全部正式测量低于 10 ms 时，新版也保持该上限 |
| 新正式案例 | 每轮满足预期，完整 API 冷上下文的全部测量严格小于 10 ms，不允许新增未认证欠账 |
| 诊断 | 独立预期和状态，不混入普通性能与成功统计 |
| 几何 | 距离按既有严格 1e-5 合同；角度类型、范围和 0.1° 约束按各验证器合同，不互换 |
| 证据 | 完整案例/轮次、有限正耗时、版本/数据/参数/环境/文档一致，远端发布核验完成 |
| 审批 | 指定管理员、当前提交、非作者审批；最后一次有效 review 状态决定是否批准 |

真实基线应按“数据哈希＋入口＋参数＋合同”注册，不仅按编号注册；原集不得缩减。严格保留 AGENTS.md 中的退休验证器和数据禁测要求。

现有 50 ms 算法预算和历史超过 10 ms 的案例不因本规范变成合格。全量当前基线需要重新批准。快速失败变成成功认证时，5% 仍适用；超限进入独立能力扩展/基线审批，不自动豁免。基线更新不得与普通功能提交捆绑自我批准。

## 3. 计时与复测

20 轮冷上下文、固定机器/SDK/工具链/编译选项、串行交替 A/B，记录全部候选、失败回退、优化、准备及验证。IO、SDK 初始化、框架独立审计和报告生成单列，不包含在算法耗时。父阶段包含子阶段时不重复累加。

固定随机种子的 2000 次配对 bootstrap 给出每例中位耗时比的 95% 区间；上界 ≤1.05 才通过，下界 >1.05 为失败，其余为 INCONCLUSIVE。置信水平为逐例统计约定，不宣称全数据集的同时覆盖率。环境漂移整批无效，不能删除单个慢样本；最多一次预先规定的复测，保留两批。跨硬件数值不可直接比较，普通 Windows 测量不构成硬实时保证。

## 4. 演练命令

仓库根目录，Python 3.11 或以上，无第三方依赖：

```powershell
python -m unittest discover -s tests/closed_loop -p 'test_*.py' -v
$closedLoopHead = git rev-parse HEAD
$closedLoopBase = git rev-parse origin/main
$closedLoopTree = git rev-parse 'HEAD^{tree}'
python tests/closed_loop/run.py smoke --head $closedLoopHead --base $closedLoopBase --merge-tree $closedLoopTree --output out/protocol-local-001
python tests/closed_loop/run.py verify --bundle out/protocol-local-001
```

上述本地命令采用本地树，仅用于离线协议检查；发布到 PR 前必须换成 GitHub 当前预合并树。输出目录必须不存在，禁止覆盖旧批次。`evaluate` 的 FAIL/INCONCLUSIVE 返回 2，缺字段、非有限数据及封装错误也非零。

```powershell
python tests/closed_loop/github_bridge.py --git-credential context --pr <PR编号> --output out/pr-context.json
```

读取 context 后将其三个 SHA 传给 `smoke`，生成新目录，再发布：

```powershell
python tests/closed_loop/github_bridge.py --git-credential publish-dry-run --pr <PR编号> --bundle out/protocol-pr-001 --receipt out/protocol-pr-001-receipt.json --allow-network-writes
```

发布器只接受固定十个合成文件，重算结果并校验实时 PR，写入同仓库 `validation-evidence` 分支的唯一目录，按 Git blob SHA 校验远端内容；只写 `closed-loop/dry-run` 状态。上传后 main/head 变化则拒绝放行，已发布旧证据仍可审计。回执置于包外，防止改变封装哈希。令牌只在内存，不能进入日志或远端文件。

GitHub Actions 的 `closed-loop/protocol-smoke` 使用只读令牌，PR 代码不持有发布凭据。上传 artifact 只作运行附件，永久证据另由受信任操作者显式发布。没有自动评论、自动请求评审、自动批准、自动合并或修改仓库可见性的代码。

## 5. 后续激活清单

1. 确认仓库真实可见性及允许发布材料；当前远端读取为 public，与历史“保持私有”说明不一致，本轮没有改变可见性。SDK、许可证和凭据绝不上传。
2. 指定可独立批准 PR 的管理员。默认 CODEOWNERS 为 KoinoSora；若 PR 作者也是该账号，必须增加另一位有权限的管理员，不能自批。
3. 管理员批准完整基线注册表，包含冻结 anchor、每例哈希/参数/合同、历史欠账、环境和批准记录。空表不得启用。
4. 部署专用、可清理的 Windows 隔离执行器。PR 构建与持有发布凭据的控制器隔离；公开仓库不得直接将任意 PR 放到个人 SDK 工作机器执行。
5. 将当前合法测试适配为本文证据合同，独立审计覆盖各自证明范围。一般兜底的 interior_sample 不能冒称全域证明；需要统一全域证明时必须先补齐实现。
6. 控制器固定从批准的 main/受保护控制版本读取规则、注册表和判定器；不能信任 PR 自带策略或 self-reported trusted_collector。获取真实 PR/review、发布回执及可信执行器身份。
7. 单独评审正式工作流、GitHub App 身份和分支规则后才允许 real_tests_enabled=true。普通 PR 不得修改规则使自身通过。当前演练 workflow 会拒绝开关为 true，防止只改开关即意外启动。
8. 远端 main 要求 PR、指定代码所有者批准、旧批准失效、主分支更新失效、最终 `closed-loop/gate` 成功，限制管理员绕过；检查来源绑定受信任 App。所有分支任务必须由始终执行的最终检查核对，缺测或 skipped 不能作为成功。

本轮不安装 self-hosted runner、不配置正式必需状态、不改变 main 保护规则，不会把未批准基线锁成无法解除的生产门禁。当前是可审阅的 bootstrap PR；正式强制执行必须完成上述激活步骤。

## 6. 影响评估

没有修改 `src/`、`include/`、CMake、固定几何输入或算法参数，因此不产生算法行为或运行时开销；未测量任何真实性能收益。新增开销仅在协议检查和发布阶段，Python 标准库足够。

Git 工作流改为工作分支＋PR，原“修改后推送”的授权继续用于工作分支，不表示允许自行合并 main。证据保存在独立分支，报告提交不会改变被测代码 SHA。管理员可独立回退本 PR 的工作流和 Python 文件，产品算法无需回滚。

风险与限制：模拟审批不是真实审批；Git 文件哈希不认证生产者身份；20 轮统计不能证明硬实时；首次框架运行来自待审 PR，仅是 bootstrap；正式隔离执行器与基线尚未部署。以上未完成项必须继续出现在交付说明中。
