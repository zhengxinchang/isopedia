我需要你再理解src中代码的前提下，修改isopedia isoform输出的表格的格式。

RC_FSM_JC	与查询转录本 完整 junction chain 匹配 的 FSM read 数；代表主要的直接 splice-structure 支持。
RC_FSM_JC_TSS	RC_FSM_JC 中，同时满足查询转录本 TSS 边界 的 read 数。
RC_FSM_JC_TES	RC_FSM_JC 中，同时满足查询转录本 TES 边界 的 read 数。
RC_FSM_JC_TSS_TES	RC_FSM_JC 中，由同一条 read 同时支持 TSS 和 TES 的 read 数；代表最严格的 full-transcript direct evidence。
RC_EM_ISM	通过 EM 算法分配给查询转录本的 ISM read 期望数；仅用于 abundance estimation，不用于定义 transcript presence。
CPM_FSM_JC	基于 RC_FSM_JC 计算的 CPM；表示仅由完整 junction-chain 直接支持得到的 normalized abundance。
CPM_FSM_JC_TSS_TES	基于 RC_FSM_JC_TSS_TES 计算的 CPM；表示同时有完整 junction chain 和 TSS/TES 支持的严格 normalized abundance。
CPM_EST	基于 RC_FSM_JC + RC_EM_ISM 计算的 CPM；表示结合 FSM 和 EM-assigned ISM 的 estimated transcript abundance。
FRAC_FSM_JC_TSS	RC_FSM_JC_TSS / RC_FSM_JC；表示 FSM junction-chain reads 中支持 TSS 的比例。
FRAC_FSM_JC_TES	RC_FSM_JC_TES / RC_FSM_JC；表示 FSM junction-chain reads 中支持 TES 的比例。
FRAC_FSM_JC_TSS_TES	RC_FSM_JC_TSS_TES / RC_FSM_JC；表示 FSM junction-chain reads 中同时支持 TSS 和 TES 的比例。


其中最核心的三类信息可以理解为：
- Direct splice-chain evidence：RC_FSM_JC
- Strict full-transcript evidence：RC_FSM_JC_TSS_TES
- Estimated abundance：CPM_EST

多外显子 FSM 的 junction-chain 匹配由 `-f/--flank` 控制；单外显子 FSM 的两端匹配由 `-F/--mono-exon-wobble` 控制。


输出的格式也要从原来的：
CPM:COUNT:FSM_CPM:FSM_COUNT:EM_CPM:EM_COUNT 
进行对应的修改


具体修改计划（基于当前 src/ 数据流）
==================================

一、目标和统计口径

1. 每个 sample 列按以下固定顺序输出 11 个以冒号分隔的数值；FORMAT 声明和实际数值必须严格一一对应，不再带目前只声明、未写值的 INFO：
   RC_FSM_JC:RC_FSM_JC_TSS:RC_FSM_JC_TES:RC_FSM_JC_TSS_TES:RC_EM_ISM:CPM_FSM_JC:CPM_FSM_JC_TSS_TES:CPM_EST:FRAC_FSM_JC_TSS:FRAC_FSM_JC_TES:FRAC_FSM_JC_TSS_TES

2. 多外显子转录本：先由 `-f/--flank`（当前默认 0 bp）判定完整 junction chain。候选的 junction 数必须与注释相同，按顺序逐对比较 donor/acceptor 坐标，两个端点各自满足 `abs_diff <= flank`；不以 TSS/TES 是否匹配作为 RC_FSM_JC 的条件。现有 `find_fsm()` 依赖位置结果交集及 splice-site 数量，需要在加载 PTIR 后补一次有序 junction-pair 核验，避免较大 flank 下“同数目、但非同一完整链”的假阳性。

3. 对每个通过完整链核验的 PTIR，按 `sample_offset_arr` 和 `sample_evidence_arr` 遍历该样本的 `isoform_reads_slim_vec`，每条 read 恰好给该转录本的 RC_FSM_JC 加 1；同一条 read 独立计算 TSS 命中、TES 命中，并分别更新 RC_FSM_JC_TSS、RC_FSM_JC_TES。仅当这条 read 两端均命中时才更新 RC_FSM_JC_TSS_TES；不能用 TSS 数和 TES 数的最小值或乘积推断交集。四个 RC 字段使用整数计数。

4. 终点命中使用当前 CLI 的 `--tss-wob`、`--tes-wob`（默认各 50 bp），边界包含等号：正链 TSS=read.left 对 transcript.start、TES=read.right 对 transcript.end；负链 TSS=read.right 对 transcript.end、TES=read.left 对 transcript.start。分别使用 `abs_diff <= tss_wob/tes_wob`，方向由注释转录本的 strand 决定。任一端不命中不影响 RC_FSM_JC；这是将“直接 splice-chain 支持”和“严格全长支持”分开的关键。

5. 单外显子没有 junction chain，采用两端匹配的替代定义：`-F/--mono-exon-wobble`（默认 50 bp）既将注释区间两侧各扩展指定 bp 以检索候选 read，也作为 FSM 边界阈值。PTIR 必须是 `is_mono_exonic()`，且区间起点和终点分别与注释转录本对应边界相差不超过该 wobble（包含等号），才计入 RC_FSM_JC；`-f/--flank` 不参与单外显子 FSM 判定。随后仍逐条 read 用 strand-aware 的 tss_wob/tes_wob 计算三个终点字段。与注释区间内部相容但不满足上述 FSM 条件的单外显子 reads 沿用现有 mono-exon EM 路径。默认 mono_exon_wobble=50、tss_wob=tes_wob=50 时，单外显子 FSM 的三个终点比例通常都是 1；阈值不同时则不能作此假定，应写入用户文档。

6. RC_EM_ISM 取 EM 收敛后的 `abundance_cur`，应用现有 `min_em_abundance` 阈值后写出，可为小数。多外显子只把与转录本连续、真子集的 junction chain 送入 EM；完整链 reads 即使 TSS/TES 不命中，也属于 RC_FSM_JC，不能因终点不命中而丢弃或再进入 EM。RC_EM_ISM 不参与 transcript presence 判定。当前 `overall_fsm_ptrs` 随转录本顺序逐步建立，计划改为每个 group 先收集所有转录本的 FSM PTIR offset，再构建 ISM/MSJC 候选，消除先进入 EM、后被识别为 FSM 的顺序依赖。相同完整链可支持多个查询转录本，RC_FSM_JC 是“每个转录本的支持数”，不保证跨转录本可相加为唯一 read 数。

7. 同一 sample 的三个 CPM 共用 index 中 `DatasetInfo.sample_total_evidence_vec[sid]` 作为分母，遵循此前要求的 index read/evidence 总数，而不是当前 `GlobalStats.fsm_em_tx_abd_total[sid]`（它受输入 GTF 和 EM 阈值影响）：
   CPM_FSM_JC = RC_FSM_JC / index_total_evidence * 1,000,000
   CPM_FSM_JC_TSS_TES = RC_FSM_JC_TSS_TES / index_total_evidence * 1,000,000
   CPM_EST = (RC_FSM_JC + RC_EM_ISM) / index_total_evidence * 1,000,000
   计算时把整数分母直接转为 f64，不先转 f32，以免大样本整数读数被舍入。index_total_evidence 来自 profile 记录的 evidence 累加，不代表原始 BAM 中所有未过滤 reads；分母为 0 时三个 CPM 输出 0。

8. 三个比例分别用对应 RC 除以 RC_FSM_JC；RC_FSM_JC 为 0 时都输出 0，不输出 NaN/Inf。必须满足 `0 <= RC_FSM_JC_TSS_TES <= min(RC_FSM_JC_TSS, RC_FSM_JC_TES) <= RC_FSM_JC`，比例均在 [0,1]。`--no-check-tss-tes` 在新模型里不再用于跳过统计：RC_FSM_JC 本来就不依赖终点，而终点三个字段必须如实计算。保留此 CLI 参数作兼容入口，执行时提示它已无效，后续版本再移除，不能把未检查的 reads 当作 TSS/TES 支持。

9. 当前 `SingleRead::process()` 的聚合 signature 由染色体和 junction 坐标构成，不含 read strand；PTIR 的每条 read 虽保存 strand，现有 FSM/EM 候选分配并不按它过滤。本轮保持原有结构匹配口径，终点方向只按注释转录本 strand 解释；不要在直接计数侧单独增加 strand 过滤而让 EM 继续混合两条链。如将来要求“read 与注释必须同向”，需同时改 PTIR/MSJC 的分链计数和 EM 分配，作为独立行为变更验证。

二、按源码数据流修改

1. `src/cmd/profile.rs`、`src/reads.rs` -> `src/cmd/merge.rs`、`src/ptir.rs` -> `src/cmd/index.rs`、`src/bptree.rs` 已提供样本 evidence、每条 read 的 left/right、聚合 junction chain 及 offset；本次直接复用，不修改 profile、index 的磁盘格式，也不要求用户重建 index。`src/gtf.rs` 已提供 transcript.start/end、splice junctions、strand。

2. `src/grouped_tx.rs` 的 `GroupedTx::update_results()`：将完整链候选的核验和 FSM offset 收集放在 EM 候选建立之前；把 `TxAbundance::update_fsm_evidence_count()` 从“只有 TSS 与 TES 同时合格才增加一个 fsm_abundance”改为单次 read 遍历累积四个独立 RC 向量。mono-exon 的 `check_mono_exon_fsm()` 使用 mono_exon_wobble 控制两端匹配，并与多外显子共用 read-level 终点统计函数。`MSJC::new()`、`m_step()` 继续负责 ISM 的 EM 丰度，必要时仅调整候选过滤，不更改 EM 迭代公式。

3. `TxAbundance` 增加四个按 sample 排列的直接计数向量；旧 `fsm_abundance` 改为清晰的 RC_FSM_JC 含义。`TxAbundanceView::encode()/from_bytes()` 与 `TmpOutputManager` 的临时 shard 数据必须同步携带这四个向量，保证原始 GTF 顺序重排后仍可写出逐样本数据；临时格式是运行过程内部文件，无需迁移已建 index。EM 向量仍独立存放。

4. `src/cmd/isoform.rs` 修改 `FORMAT` 字符串与 sample 表头；`src/grouped_tx.rs` 的 `TxAbundanceView::write_line_directly()` 按上面 11 字段一次性写出真实数值。`src/results.rs`/`src/myio.rs` 的读写和 `isopedia-tool output --mode merge-isoform` 要按新 FORMAT 解析，不再假定旧样本字段 `[1]` 是整数 COUNT；合并时依据 RC_FSM_JC 重算检测/排名，遇到新旧格式混合时明确报错。

5. 固定列中现有 `detected(total:fsm:em)` 和 `n_pos_samples(total:fsm:em/sample_size)` 会把 EM 当作 presence。改为 `detected` 和 `n_pos_samples/sample_size`，分别表示是否有任一样本达到 `RC_FSM_JC >= min_read`、以及达到该阈值的样本数；校验 `min_read >= 1`，避免阈值为 0 时零 read 样本也被判阳性。`ranking_score` 也以 RC_FSM_JC 和 index 样本总 evidence 计算，避免 EM-only 转录本获得“直接证据”排名。排名公式沿用现有阳性样本比例、CPM 几何平均和 Gini 组合，但直接计数改用可容纳大样本的整数类型，并在计算中转成 f64、处理零总数。`src/global_stats.rs` 与命令末尾的统计日志同步区分直接检测和估计丰度，不用 EM-positive 数作为 transcript presence。

6. 更新 `README.md` 的 11 字段表、三个 CPM 的共同分母、单外显子特殊定义、`-f` 与 tss_wob/tes_wob 的不同职责，以及旧格式不兼容说明；顺手把 README 中过时的 isoform `-f` 默认值 10 改为当前源码的 0。`src/cmd/splice.rs`、`src/cmd/fusion.rs` 的输出不使用 isoform 的 FORMAT，保持原行为。

三、验证与完成标准

1. 为单个多外显子 PTIR 构造四类 reads：仅 TSS 命中、仅 TES 命中、两端都命中、两端都不命中；验证 RC_FSM_JC=4、各子计数及联合计数正确，终点失败的 FSM 仍在 RC_FSM_JC 中且不进入 EM。测试正负链使用不同的 tss_wob/tes_wob，并测试 `-f` 边界等号、超 1 bp、junction 数相同但顺序/配对不同。

2. 测试单外显子 mono_exon_wobble=0 与非零、两端刚好处于阈值及超出 1 bp、内部相容但非 FSM 的 EM 路径；验证更改 `-f` 不影响单外显子 FSM。测试同一 PTIR 对多个转录本的 FSM/ISM 分类与转录本遍历顺序无关。检查整数计数和三个比例的不变量，以及分母为 0、RC_FSM_JC 为 0 时输出全为有限数。

3. 测试临时 shard encode/decode、11 个 FORMAT 名与 11 个 sample 值逐项对齐、旧新格式混合合并报错、新格式多 sample 合并后检测状态/排名正确；运行 `cargo test` 和 release 构建。对同一 index/GTF 保存修改前后 isoform 输出，按 transcript/sample 对照 RC_FSM_JC、严格联合计数和 RC_EM_ISM；确认 splice/fusion 输出未受影响。成功标准是 EM-only 转录本可有 CPM_EST，但 `detected=no`，且终点失败的 FSM 仍贡献 CPM_FSM_JC。
