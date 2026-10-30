# v1.7.4：为 `isoform --info` 输出 read 坐标和方向

本文件记录基于 `v1.7.3`（`b2bc12a`）的修改计划及完成验证。

## 1. 输出约定

已按用户确认的原计划实现：`INFO` 只列出计入该 transcript、该 sample 的 `RC_FSM_JC` 的 reads。多外显子在 `-f 0` 时，这些 reads 的完整 junction 链均精确匹配；TSS/TES 可以不同。EM reads 不包含在内。

- 未指定 `--info`：保留现有 13 个 FORMAT 字段、样本数值及检测结果。
- 指定 `--info`：在现有 FORMAT 末尾追加 `:INFO`，在每个样本值末尾追加相应内容。
- 单条 read 为 `LEFT|RIGHT|STRAND`；多条 read 用逗号分隔；没有 FSM read 时为 `.`。
- 示例：`100|500|+,110|510|-`。只输出这三个值，无第四个字段、无尾部 `|`。
- `LEFT`、`RIGHT` 直接使用索引中 read 的主比对坐标，为 0-based、左闭右开 `[LEFT, RIGHT)`；负链也不交换左右端点。方向使用 read 自身的 `+` / `-`，不使用查询 transcript 的方向。
- 包含 exact 和 wobble-only FSM reads，也包含没有通过 TSS/TES 检查的 FSM reads；不额外使用 `min_read` 过滤明细。
- 一条索引 read evidence 对应一项，条目数应等于 `RC_FSM_JC`。坐标相同的不同 reads 保留重复项；沿用已有匹配去重逻辑，避免同一 PTIR 因重复搜索被添加多次。
- 同一 read 可以支持多个查询 transcript，分别出现在对应行中。沿用现有遍历顺序，不新增 read 排序或唯一归属算法。
- 不需要 `profile --rname`；不输出 read 名称、MAPQ、附加 tags、补充比对或 EM 分配权重。

## 2. 已有实现与复用位置

- `src/cmd/isoform.rs` 已定义 `AnnIsoCli.info`；`ISOFORM_FORMAT` 是当前 13 字段格式，`TableOutput::new` 接收最终格式字符串。
- `src/ptir.rs` 的 `ReadDiffSlim` 已包含 `left`、`right`、`strand`；`sample_offset_arr` 和 `sample_evidence_arr` 可定位每个样本的 read 切片。
- `src/grouped_tx.rs::update_fsm_evidence_count` 是多外显子、单外显子 FSM 计数的共同入口，可在已有 PTIR 已加载时收集明细，无需再次查询或读取索引。
- 数据现有路径为 `TxAbundance → TxAbundanceView::encode/from_bytes → TmpOutputManager → write_line_directly`，明细必须走完整路径，才能在临时分片排序后保留。
- 复用 `Strand` 的 `Display`、已有字节缓冲区和 `write!`。不调用 `ReadDiffSlim::to_string_no_offsets()`，因为它输出第四个附加信息字段；也不修改这个公共方法的既有语义。

## 3. 最小修改步骤

### 第一步：在 FSM 计数时收集明细

修改 `src/grouped_tx.rs`：

- 在 `TxAbundance` 增加 `fsm_read_info: Vec<(u32, Vec<u8>)>`，按 sample ID 有序保存非空样本的 INFO 文本，初始化为空。采用与现有 `fsm_counts` 一致的稀疏布局，避免给所有样本分配字符串。
- 在 `update_fsm_evidence_count` 中，仅当 `cli.info` 为真时创建或查找该样本的缓冲区，直接追加 `LEFT|RIGHT|STRAND` 和必要的逗号，不克隆整个 `ReadDiffSlim` 或 `PTIR`。
- 多外显子复用已有 read 遍历，避免为 INFO 再扫描一次。单外显子仅在开启 INFO 时遍历该样本的 read 切片，逐条取得实际方向；其原有批量计数逻辑保持不变。
- 更新 `TxAbundance::get_mem_size` 对新增容器和文本容量的统计。

验证：开关 INFO 前后 FSM/TSS/TES/exact 计数不变；多个 PTIR 累积到同一样本时无遗漏、无多余分隔符；关闭 INFO 时不构造 read 文本、不为明细分配堆内存。

### 第二步：通过现有临时分片传递明细

仍修改 `src/grouped_tx.rs`：

- 为 `TxAbundanceView` 增加对应字段，在现有 `encode/from_bytes` 中追加编码和解码：非空样本项数，以及每项的 sample ID、文本字节长度和文本内容。
- 沿用现有显式长度编码方式，不引入新序列化框架或辅助对象。新增长度转换、读取边界和 sample ID 检查，防止溢出或截断记录被误读。
- 未开启 INFO 时只编码零项数量，不写 read payload；接受这一小段固定临时记录开销，避免增加两套编码逻辑。
- `TmpOutputManager` 已使用可变长度记录，继续使用现有写出、排序和读回流程，不新增旁路文件或第二次索引查询。
- 同步调整现有测试中的结构体初始化，以及依赖临时记录尾部长度的字节偏移断言。

验证：明细经过 encode/decode 和多个临时分片后，仍对应正确的 transcript 和 sample；截断的新增 payload 返回错误。这里只改变本次运行的临时记录，不改变已有 profile、PTIR archive 或索引格式，不要求重建索引，也不添加旧临时文件迁移逻辑。

### 第三步：追加输出并说明格式

- `src/cmd/isoform.rs`：保留 `ISOFORM_FORMAT`，在构造 `TableOutput` 时根据 `cli.info` 追加 `:INFO`；更新 `--info` 帮助和样本格式说明。
- `src/grouped_tx.rs::write_line_directly`：已有 13 项照常输出，再按 sample ID 顺序追加 INFO 或 `.`；同时处理现有全零样本快速路径，避免它因提前 `continue` 漏掉 `:.`。不为每个输出样本重新线性搜索全部 INFO 项。
- `README.md`：将 Isoform Query 的 FORMAT 更新为代码现有的 13 字段，说明 `--info` 条件追加字段、FSM 范围、坐标约定、分隔符和空值，给出一个命令和一个样本值示例。该表内与现有检测字段直接冲突的描述同步修正，其他章节不扩展。

验证：不开启 INFO 时数据行和 FORMAT 与原行为一致；开启时为 14 项，前 13 项及其他输出列不变，全零样本、仅 EM 样本的 INFO 均为 `.`。

## 4. 最小回归检查

在 `src/grouped_tx.rs` 的现有 `#[cfg(test)]` 模块内复用 `test_cli`、`test_tx`、`test_ptir` 和现有输出/分片测试，不新增测试框架或独立测试工程。少量测试覆盖：

1. 多外显子和单外显子、正负链、重复坐标、多样本及多个 PTIR：每项内容正确，数量等于 `RC_FSM_JC`，不存在跨样本串数据；TSS/TES 不通过的 FSM read 仍在 INFO 中。
2. 开关 INFO 的最终输出比较：13/14 字段、数值列一致、无 FSM 时为 `.`，包括全零快速路径和仅 EM 的情况；关闭时明细容器为空。
3. 复用多分片测试检查明细保留和行顺序；补充新增 payload 截断的失败检查。复用现有 mono-exon 分区去重案例检查条目不重复。

实现时先补能复现当前缺失输出的检查，再完成代码。运行 `cargo test --lib grouped_tx::tests`，通过后运行 `cargo test --lib` 和 `cargo check --bins`；用 `git diff --check` 检查补丁。格式化只涉及本次改动，避免全仓库格式调整。

## 5. 范围和成本

预计只修改 `src/grouped_tx.rs`、`src/cmd/isoform.rs`、`README.md` 三个文件，测试留在原模块。不新增模块、依赖、通用 helper 或配置项，不修改 FSM 匹配、EM 算法或版本号。

开启 INFO 后，内存、临时文件和结果文件会随输出的 read 条目总数增长，read 同时支持多个 transcript 时也会重复存储。沿用现有按 transcript 数量分片的机制，不声称它能限制 read 明细的内存。实现处增加简短 `ponytail:` 注释说明这一上限；若实测高覆盖数据的内存不足，再考虑按字节预算刷分片，本次不提前设计新存储层。

## 6. 完成验证

- 首先新增跨分片输出测试，在原实现上确认 `--info` 缺失明细导致失败；实现后通过。
- `cargo test --lib grouped_tx::tests`：29 个测试通过；`cargo test --lib`：全部 39 个测试通过。
- `cargo check --bins`、`cargo build --release --bin isopedia` 和 `git diff --check` 通过。
- 合成 BAM 经 `profile → merge → index → isoform -f 0` 完整运行，使用两个样本、20 个重叠多外显子 transcript 和一个单外显子 transcript，包含正负链、重复坐标和 TSS/TES 不完全支持的 reads。
- 开启 INFO 后输出 2,000,002 条 read 明细，逐行验证数量与 `RC_FSM_JC` 一致，坐标/方向正确，原 13 个样本字段及其他固定列与关闭 INFO 时完全相同。
- Release 构建、查询线程数 2、本机单次测量：关闭/开启 INFO 分别为 0.10/0.17 秒，峰值 RSS 为 21,424/39,036 KiB，采样观察到的临时分片总量为 3,160/20,003,636 字节。合成数据重复度较高，以上不代表真实数据或持续吞吐基准。
- 现有 `merge-output` 对 isoform 表严格要求 13 字段，不接受新增 INFO 的表；本次没有扩展该命令，README 已说明。版本号保持 `1.7.3`，尚未执行发布操作。

### BRCA1 实际输入验证

- 输入：`/ssd2/projects/isoatlas/workspace/hg002/selection_BRCA1_novel.gtf`，转录本 `PB.37601.8`（23 个外显子）。
- 使用此前 HG002 查询记录中的索引 `/hdd1/isopedia_datadownload/isopedia_index_v1.1_1035`，运行当前 Release 二进制，分别使用 `-f 0` 和 `-f 0 --info`。
- 结果：1,035 个样本中 146 个有 FSM 支持，共 290 条 INFO read；`RC_FSM_JC_EXACT` 合计 290，`RC_FSM_JC_WOBBLE_ONLY` 合计 0。
- 逐样本验证明细条数等于 `RC_FSM_JC`，原 13 个样本字段及固定列在开关 INFO 前后完全一致。根据明细坐标重新计算 TSS/TES 支持也一致：TSS 274、TES 124、两端同时支持 120。585 个仅有 EM 支持的样本均输出 `INFO=.`。
- 结果、命令、日志及可重跑的 `validate.py` 保存在 `target/validation/brca1-f0-s22wiu0q/`。带 INFO 的结果为 `BRCA1.f0.info.tsv.gz`，汇总为 `summary.json`。
- 两次运行耗时分别为 8.72 秒和 0.40 秒，峰值 RSS 约 405 MiB；第二次运行可能受文件缓存预热影响，不能据此判断 INFO 更快。
