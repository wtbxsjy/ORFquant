# ORFquant 多线程并行修复记录

## 问题概述

ORFquant 使用 `parallel::mclapply()` (Unix fork) 进行并行基因处理。当 `genome_seq` 是 `FaFile_Circ` 对象时，fork 出的子进程继承了父进程的文件描述符。子进程退出时触发 R GC → `FaFile$finalize()` → 底层 Rsamtools C++ 对象清理失败。

### 错误表现
```
Error in x$.self$finalize() : attempt to apply non-function
```
每个 fork worker 退出时均报此错误，共 ~332 次（约 20252 基因 / worker 批次数）。

## 修复尝试

### v1: DNAStringSet 预加载 (commit 00a3cdf)

**方案**: 在 `load_annotation()` 中，`getSeq(fafile)` 将 FaFile 转换为 DNAStringSet（纯 R 数据，无文件句柄）。

**结果**: ❌ 染色体名不匹配。`getSeq()` 返回的 DNAStringSet 使用长名称（如 `1 dna:chromosome chromosome:IRGSP-1.0:1:1:43270923:1 REF`），而 annotation 期望短名称（`1`, `2`, `Mt`）。

### v2: 修复 seqnames (commit 9eb061b)

**方案**: 在 `getSeq()` 前保存 `seqnames(seqinfo(fafile))`，getSeq 后恢复为短名称。

```r
orig_seqnames <- seqnames(seqinfo(genome_sequence))
genome_sequence <- getSeq(genome_sequence)
names(genome_sequence) <- orig_seqnames
```

**结果**: ❌ 332 finalizer errors。FaFile 仍在某处被引用。

### v3: load() 本地环境隔离 (commit f7f7a67)

**方案**: `load(path)` 默认加载到全局环境，造成 FaFile 全局副本。改为 `load(path, envir = _env)` 加载到本地环境，`getSeq()` 转换后通过 `<<-` 只将 DNAStringSet 赋值到全局。

```r
_env <- new.env(parent = emptyenv())
GTF_annotation <- get(load(path, envir = _env))
```

**结果**: ❌ 332 finalizer errors。`load(envir=...)` 确认为本地加载无误（无全局泄漏），但仍无法消除 FaFile 引用。

### v4: mc.cleanup=FALSE (commit 待测试)

**方案**: 设置 `parallel::mclapply(..., mc.cleanup = FALSE, mc.silent = TRUE)`。不等待子进程清理，可能抑制 finalizer 触发。

**结果**: 🔄 待测试

## 失败的方案

### SnowParam v1 — 包级别 helper 函数 (commit 072b4df)

**方案**: 使用 `BiocParallel::SnowParam(workers=n_cores, type="SOCK")` 替代 `mclapply`。独立 R 进程不继承 FaFile 文件描述符。

**结果**: ❌ `Error in .orfquant_load_worker_namespaces(): could not find function ".orfquant_load_worker_namespaces"` — SnowParam SOCK worker 是独立 R 进程，`process_gene_chunk` 闭包只捕获了数据变量 (`worker_annotation`, `genome_ref`, `worker_genome`, `process_gene`)，但包级别函数 `.orfquant_load_worker_namespaces()` / `.orfquant_open_genome_ref()` / `.orfquant_close_genome()` 的查找路径指向 ORFquant namespace，worker 未加载 ORFquant 所以找不到。且 `.orfquant_load_worker_namespaces()` 缺少 `library(ORFquant)`，即使函数可见，`FaFile_Circ` 和 `process_gene` 内调用的 ORFquant 内部函数也不可用。

### SnowParam v2 — inline 所有 helper (commit 待推送)

**方案**: 将所有 helper 函数 inline 到 `process_gene_chunk` 闭包内部，消除所有包级别函数查找；添加 `require("ORFquant")` 确保 `FaFile_Circ` 和核心 ORF 检测函数在 worker 中可用。

**结果**: 🔄 待测试

### forge_BSgenome=TRUE

**方案**: 在 PREPAREANNOTATION 中 `forge_BSgenome=TRUE` 构建真正的 BSgenome 包。

**结果**: ❌ 不可行。`annotation_name` 必须是已注册的 NCBI assembly（如 `IRGSP-1.0`），但 rice 等非模式生物不在 GenomeInfoDb 注册列表中。

### getSeq() + GTF_annotation$genome 替换

**方案**: getSeq 后再替换 `GTF_annotation$genome`。

**结果**: ❌ FaFile C++ 对象不跟随 R 层面的赋值，仍被 fork 子进程继承。

### close(FaFile)

**方案**: `close(fafile)` 关闭文件句柄。

**结果**: ❌ FaFile 不是 open 状态（lazy open）。

## 当前已应用的修复

以下修复已合并到 dev 分支并推送到远端：

1. **`disjointExons` → `exonicParts`** (line 1791, 5132): Bioc 3.20 兼容
2. **`load_annotation` FaFile 兼容** (line ~4930): 检查 `genome_package` 先再 fallback 到 `genome`；FaFile 时设置 `genome_ref` 而非预加载 DNAStringSet
3. **`load()` 本地环境** (line ~4920): 避免全局 FaFile 副本
4. **SnowParam socket worker** (line ~4545): 替代 `mclapply` fork，每个 worker 独立打开 FaFile
5. **Inline worker helpers** (line ~4566): `process_gene_chunk` 内联所有包加载和 FaFile 操作，消除 SnowParam SOCK 序列化的函数查找问题；显式 `require("ORFquant")` 确保 `FaFile_Circ` 和核心 ORF 检测函数可用

## 备选方案

- **n_cores=1**: 单线程处理，Nextflow 在 pipeline 层面并行 23 个样本
- **FaFile → TwoBitFile**: 使用 rtracklayer::TwoBitFile 替代，是否有相同问题需验证

## 文件改动位置

| 行号 | 改动 |
|------|------|
| 1791 | `exbin <- exonicParts(orfann, linked.to.single.gene.only = FALSE)` |
| 4085-4153 | `.orfquant_genome_ref` / `_open_genome_ref` / `_close_genome` / `_load_worker_namespaces` helper 函数 |
| 4212-4222 | `parallel_backend` 参数 + auto→snow 逻辑 |
| 4495 | `process_gene` 接受 `annotation` / `genome_sequence` 参数 |
| 4544-4605 | SnowParam socket 分支 (inline worker helpers + bplapply) |
| 4606-4623 | fork 分支 (mclapply, 保留为显式 legacy 选项) |
| 4920-4936 | `load_annotation`: 本地环境 + FaFile→genome_ref |
| 5132 | `nsns <- exonicParts(annotation, linked.to.single.gene.only = FALSE)` |
