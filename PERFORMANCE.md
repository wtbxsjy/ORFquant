# ORFquant 性能优化记录（v1.3.3）

## 基准数据

采用 nf-core/riboseq 流程的测试数据（`nf-core/test-datasets` `modules` 分支，
`data/genomics/homo_sapiens/riboseq_expression`）：

| 基准 | 注释 | P-sites | 规模 |
|------|------|---------|------|
| A（真实） | GRCh38.111 chr20 GTF + chr20 FASTA | `SRX11780887_chr20.bam`（nf-core/riboseq 比对结果）经 `prepare_for_ORFquant` | 296 个基因区域，7 035 个 P-site 位点 |
| B（高深度模拟） | 同上 | chr20 上约 70% 蛋白编码基因（每基因 1–2 个异构体）的 CDS 模拟 3-nt 周期信号 + 背景 + junction reads | 195 661 个位点 / 48 万 reads |

A 反映真实数据格式，但测试 BAM 为下采样（平均每区域约 24 个 P-site），很少走到定量/注释阶段；
B 用于压测定量与注释路径。复现脚本：`inst/benchmarks/nfcore_chr20_benchmark.R`。

## 结果

参考实现 = 原 HEAD（3414a12）+ 下文 3 个最小正确性修复（原 HEAD 在用当前代码生成的注释上所有区域均失败，无法直接对比）。
同一台 4 核机器；参考与新实现成对同时运行以抵消机器负载漂移（该环境下同一程序先后两次运行可相差约 20%）。

| 基准 | 参考实现（1 核） | 新实现（1 核） | 新实现（mclapply 4 核） |
|------|-----------------|---------------|------------------------|
| A：chr20 全部 296 区域，真实 P-sites | 1636 s | 961 s（1.7×） | 249 s（6.6×） |
| B：chr20:1-3Mb，31 区域，高深度模拟 | 707 s | 402 s（1.8×） | 151 s（4.7×） |

基因组规模下的收益主要来自第 1 节的索引（chr20 子集中体现不出来），见下文微基准。

结果一致性：B 基准（chr20:1-3Mb，31 区域，42 个 ORF）上，串行 / mclapply(4 核) /
mirai（进程内 mock）三种模式输出与参考实现逐组件 `all.equal` 一致，
唯一差异是 `ORFs_tx` 中两个原本就应被删除、但在旧实现中残留为全 NA 的列
（`compatible_tx_longest`, `compatible_biotype_longest`）。

## 主要改动

### 1. 基因组规模的复杂度问题（对全基因组运行影响最大）

- `ORFquant()` 对每个区域用 `%over%` 扫描**全基因组** P-site（3 次）以及全基因组注释
  （`exons_bins`, `cds_genes`, `cds_txs`, `exons_txs`, `cds_txs_coords`），
  总代价 O(区域数 × 数据量)。
- 现在 `run_ORFquant()` 先用一次 `findOverlaps()` 建立「区域 → 元素」索引
  （`.orfquant_region_index()` / `.orfquant_region_annotation_index()`），
  每个区域只拿到自己的子集；下游仍执行原有的 strand-aware `%over%`，因此结果完全相同。
- 基因组规模微基准（5M P-site、20 万转录本/120 万外显子、2 万区域）：
  每区域子集提取 1.435 s → 0.028 s，总计约 **8 小时 → 9 分钟**。

### 2. 单区域计算热点

| 位置 | 问题 | 改动 |
|------|------|------|
| `get_orfs` | 每帧调用一次 `translate()`（每次重建 fuzzy 遗传密码表）；对每个起始密码子线性扫描终止密码子 | 三帧一次 `translate()`；`findInterval()` 二分查找下一个终止密码子 |
| `detect_translated_orfs` | 每个 ORF 单独 `translate()` 蛋白；`orfs_gr[[nam]] <- orf` 逐个增长 GRangesList（O(n²)）；ORF 结构两两 `identical()` | 每个转录本一次向量化 `translate(extractAt())`；普通 list 收集后一次合并；结构字符串键 + `outer()` |
| `detect_translated_orfs` | 三组 P-site 分别 `mapToTranscripts()` | 合并位点只映射一次 |
| `select_start` | `RleList[GRanges]` 取覆盖度；每个 ORF 生成一个 `DataFrame` 再 `rbind` | 解压一次后按坐标切片；数值矩阵一次构造 `DataFrame` |
| `calc_orf_pval` | 每个 ORF 做 `P_sites_rle[ranges(ORFs[i])]`；循环内 `ORFs$pval[i] <-`（S4 拷贝） | 向量切片；结果列在 `DataFrame` 中组装后一次赋值 |
| `from_tx_togen` | 逐个增长 GRangesList | `lapply` + 一次 `GRangesList()` |
| `select_quantify_ORFs` / `annotate_ORFs` | 对 `CompressedGRangesList` 逐元素 `[[<-` / `[[`（每次 O(总长度)） | 循环期间使用普通 list，结束时转换回 GRangesList |
| `annotate_ORFs` | 每个 ORF 对全基因组 `cds_txs_coords` 做 `as.character(seqnames())` 比较 | 先按本区域转录本预筛一次 |
| `annotate_splicing` | 每个外显子多次 `ran$spl_type <-`；循环内 `sort(c(spl_ran, ran))`（O(n²)）；逐个增长 `grliss` | 局部变量确定类型后赋值一次；收集后一次稳定排序（结果相同）；`split()` |

### 3. 并行

- `mclapply`：`mc.preschedule = TRUE` 按轮询分配任务，原来按基因组顺序分配，大基因座易集中到同一 worker。
  现按区域 P-site 数降序派发，结束后恢复原顺序（输出不变）。
- mirai（v2 与 v3）：
  - `mirai_map()` 的 worker 是定义在函数内部的闭包，序列化时会把外层帧
    （`GTF_annotation`、`genome_seq`、`for_ORFquant_data`）随**每个任务**一起发送，
    且自由变量解析到这些副本而不是 daemon 预加载的数据。chr20 小区域实测每任务 13 MB → 现 2.4 KB。
    现将 worker 的环境设为 `globalenv()`。
  - daemon 原先重新 `load(for_ORFquant_file)`：多个输入文件（向量）时出错，且跳过了主进程的合并/清洗。
    现主进程将处理后的 P-site 数据写入临时 RDS 供 daemon 读取。
  - daemon 内同样建立区域索引。
  - 删除 v3 文件末尾遗留的 `# DEBUG PATCH` 行。

## 顺带修复的正确性问题（在新数据上会导致结果为空/区域失败）

1. **`exonicParts()` 额外列**（Bioc ≥ 3.20）：`exons_bins` 多出 `tx_id/exon_id/exon_name/exon_rank`，
   `select_txs()` 中 `mcols(gene_feat)[, names(mcols(genbin))]` 报 `subscript contains invalid names`，
   **所有区域失败**。用当前代码新建的注释都会触发（旧的 rice 注释因是旧版本生成而未触发）。
   修复：`select_txs()` 只取 `tx_name/gene_id`；`prepare_annotation_files()` 只保存所需列。
2. **`compatible_with` 类型被强制为 logical**：`detect_translated_orfs` 中创建为 logical `NA`，
   新版 Bioconductor 对 `CompressedGRangesList` 逐元素赋值时会把值转换为已有列类型，
   `CharacterList` 变成 `LogicalList(NA)`，随后报 `0 elements in value to replace 1 elements`
   或 `non-character argument`（真实数据中 296 区域有 27 个因此失败）。
3. **结果为空时导出崩溃**：`ORFs_gen$type = "CDS"` 对长度 0 的 GRanges 报错。

## 回归测试

`tests/testthat/test-optimizations.R`：`get_orfs` 与旧实现逐位一致（随机序列、两种遗传密码）；
区域索引与 `x[x %over% region]` 一致；`select_start` 统计；`select_txs` 接受 `exonicParts()` 风格的 bins。
