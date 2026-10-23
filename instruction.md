# 聚合釜取样托盘的归属复核

中试线换班后，取样托盘和一排反应釜的升温程序失去了对应关系。每个托盘仍保留同一次运行的四次色谱转化率，反应釜一侧也保留完整的分段温度程序；麻烦的是混样托盘和读数异常混在同一批记录里，不能拿托盘编号硬配。

现在要修的是 `/app/polymer_reconcile.py`。它由交接班定时器调用：

```bash
python3 /app/polymer_reconcile.py \
  --input /data/polymer_batches.json \
  --output /app/out/reconciliation.json \
  --ledger /app/state/experiment-ledger.jsonl
```

入口路径可以用环境变量 `TASK_APP` 覆盖，输入目录可以用 `TASK_DATA_DIR` 覆盖。成功时标准输出和标准错误都为空。整批最多 30 个反应釜、32 个样品托盘、每个托盘 6 个读数时，程序应在 40 秒内完成。

## 托盘里还剩下什么

中试记录以 UTF-8 编码的 JSON 文档递交。根对象除 schema_version、reaction_histories、sample_bundles、measurement_policy 外不得再放成员。

缺键、额外键和同一对象内重复键都必须拒绝。`schema_version` 必须是非布尔整数 `1`。

- `reaction_histories` 是 8 到 30 个对象的数组。
- 每个对象恰好有 `history_id`、`steps`、`catalyst_scale`。
- `history_id` 是非空字符串，去首尾空白后长度 6 到 24，批次内不能重复。
- `catalyst_scale` 是大于等于 0.2 且小于等于 3.0 的有限数值。
- `steps` 是 2 到 8 个对象的数组，每个对象的键恰好是 `duration_s`、`temp_start_k`、`temp_end_k`。
- `duration_s` 必须大于 0 且不超过 7200。
- 两个温度是 280 到 500 K 的有限数值。
- 布尔值和数字字符串都不算数值。

- `sample_bundles` 是 8 到 32 个对象的数组。
- 每个对象的键恰好是 `bundle_id`、`observations`。
- `bundle_id` 是非空字符串，去首尾空白后长度 6 到 24，批次内不能重复。
- `observations` 是 4 到 6 个对象的数组。
- 每个观测对象恰好有 `elapsed_s`、`conversion`。
- `elapsed_s` 是大于 0 且不超过 10000 的有限秒数。
- `conversion` 是大于等于 0 且严格小于 1 的有限数。
- 同一托盘的 `elapsed_s` 不能重复。

`measurement_policy` 恰好有这六个键：

| 键 | 允许值 |
| --- | --- |
| `activation_energy_bounds_j_mol` | 两个有限数的数组，满足 `30000 <= 下界 < 上界 <= 180000` |
| `ln_prefactor_bounds` | 两个有限数的数组，满足 `10 <= 下界 < 上界 <= 35` |
| `conversion_floor` | 大于 0 且小于 0.01 的有限数 |
| `consensus_min_points` | 非布尔整数，至少为 3 |
| `conversion_tolerance` | 从 0.000001 到 0.02 的有限数 |
| `max_assignments` | 非布尔整数，必须等于反应釜数量，且不超过 30 |

每个反应釜和托盘都必须在各自方向上一一对应；一个托盘最多舍弃一个异常读数，至少保留 `consensus_min_points` 个读数参与一致性判断。没有对应反应釜的托盘，其全部读数都要进入异常清单。

## 归位与动力学闭合

动力学采用一级反应模型。通用气体常数固定为 `R = 8.31446261815324 J/(mol*K)`。对反应釜 `h` 的每个温度段，温度从 `temp_start_k` 按该段实际经过时间线性变化到 `temp_end_k`；超过全部温度段总时长后不继续累计。给定活化能 `Ea`，在时刻 `t` 的热剂量为

```text
D_h(Ea, t) = integral_0^t exp(-Ea / (R * T_h(s))) ds
B_h(Ea, t) = catalyst_scale_h * D_h(Ea, t)
```

观测转化率先换成 `y = -ln(1 - conversion)`。对反应釜 `h` 与托盘 `b`，在托盘每个观测时刻计算 `r_i = y_i / B_h(Ea, t_i)`，再取 `z_i = ln(r_i)`。如果该托盘任意观测的 `conversion <= conversion_floor`，或任意 `B_h` 非正，则该配对的得分为 `1000000`。

否则令 `m = n - 1`。按 `(r_i, i)` 升序排列后，从连续 `m` 个观测中选择 `max(z) - min(z)` 最小的一组；若差值相同，选择原始下标元组字典序较小的一组。这 `m` 个下标是支持点，剩余一个下标是该配对的异常点，配对得分就是这组 `max(z) - min(z)`。

反应釜和托盘数量可能不同。构造大小为 `max(history_count, bundle_count)` 的方阵：真实反应釜与真实托盘的格子填配对得分；真实反应釜到虚拟托盘、以及虚拟反应釜到虚拟托盘的格子填 `4 * conversion_tolerance`；虚拟反应釜到真实托盘的格子填 `0`。用最小总代价的一一分配同时确定托盘归属。

最优分配必须覆盖全部反应釜，否则写 `no_consensus`。选定 `Ea` 后，把所有支持点的 `ln(y_i / B_h(Ea, t_i))` 排序；奇数个取正中值，偶数个取中间两个值的算术平均，得到 `ln_prefactor`。它必须落在 `ln_prefactor_bounds` 内。对每个支持点计算 `predicted = 1 - exp(exp(ln_prefactor) * B_h)`，要求 `abs(predicted - observed) <= conversion_tolerance`，否则写 `no_consensus`。`conversion_rmse` 是全部支持点该残差平方的算术平均再开平方。

`Ea` 在 `activation_energy_bounds_j_mol` 内选择，使上述完整分配的配对得分总和最小；等价目标并列时选较小的 `Ea`，仍不能形成唯一一致分配时写 `no_consensus`。

## 复核报告

报告根层依次写出下列六个成员：`schema_version`、`kinetics`、`assignments`、`flagged_observations`、`diagnostics`、`ledger_entry_id`；不得增加其他成员或改动序列。

- `kinetics` 恰好有 `activation_energy_j_mol`、`ln_prefactor`、`conversion_rmse` 三个键。
- 前两项分别采用焦耳每摩尔和对数形式，均为有限数值。
- 第三项是非负有限数值。
- 报告数值按十进制四舍五入保留六位小数，连字符两侧都不留空格。
- `conversion_rmse` 保留八位小数。
- 恰好落在末位半值时取偶数末位。
- 复核结果要求活化能误差不超过 600 J/mol，对数指前因子误差不超过 0.020，支持点的转化率均方根误差不超过 `conversion_tolerance`。

- `assignments` 是对象数组，按 `history_id` 的 UTF-8 字节序升序排列。
- 任一条归属记录只容纳 `history_id`、`bundle_id`、`supporting_indices`、`outlier_indices` 这四项。
- 两个索引数组都是非负整数数组、严格升序，且合起来恰好覆盖该托盘的全部读数下标。
- 每个对应至少保留三个支持点，最多有一个异常点。

- `flagged_observations` 按 `bundle_id` 的 UTF-8 字节序升序，同一托盘再按 `index` 升序。
- 每项恰好有 `bundle_id`、`index`、`reason` 三个键。
- 已对应托盘的舍弃点写 `consensus_outlier`。
- 无对应反应釜的托盘所有读数写 `no_consensus_history`。

- `diagnostics` 恰好有 `history_count`、`bundle_count`、`assigned_count`、`consensus_count`、`unassigned_bundle_count` 五个键。
- 五个诊断值都为非布尔整数。
- `assigned_count` 必须等于反应釜数量。
- `unassigned_bundle_count` 必须等于托盘总数减已对应托盘数。
- `consensus_count` 是所有支持点总数。

`ledger_entry_id` 是 64 位小写十六进制字符串。报告文件必须是紧凑 JSON：不使用缩进或分隔空格，数组顺序按上面的规定，末尾恰好一个 LF。连续处理同一份输入两次，报告和交接记录都必须逐字节不变。

## 交接记录

`--ledger` 指向 JSONL。

- 账本每行由 `schema_version`、`entry_id`、`activation_energy_j_mol`、`ln_prefactor`、`assignments_sha256` 五项组成，不多也不少。
- 第一项 `schema_version` 是非布尔整数 `1`。
- 两个摘要是 64 位小写十六进制，两个动力学数值为有限数。
- 行内用紧凑 JSON，每行一个 LF，键顺序与上面相同。
- 输入内容改变并产生不同复核结果时，追加一行。
- 重复同一输入时不追加。
- 成功发布后，报告和交接记录要么同时换成新内容，要么都保留原字节。

## 被挡下的时候

输入不存在、无法读取或 UTF-8 解码失败时，退出码为 2，标准输出为空，标准错误只写一行 `reconcile refused: input_unavailable` 或 `reconcile refused: json_invalid`。

重复键写 `reconcile refused: duplicate_key`。

结构、类型或键集合不合格写 `reconcile refused: schema_invalid`。

数值范围不合格写 `reconcile refused: constraint_invalid`。

没有唯一可接受的一致对应时，退出码为 2，标准错误只写 `reconcile refused: no_consensus`。

发布报告或交接记录失败时，退出码为 3，标准错误只写 `reconcile refused: publish_failed`。

任何失败都不得留下 traceback、临时文件或半份报告。

已有报告和交接记录的字节必须保持不变，输入文件也不得被改写。

发布失败可能发生在报告已经换过一次之后，因此恢复时要把两个旧文件都放回原字节。
