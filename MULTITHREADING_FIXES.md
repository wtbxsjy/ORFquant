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

### SnowParam (BiocParallel)

**方案**: 使用 `BiocParallel::SnowParam(workers=n_cores, type="SOCK")` 替代 `mclapply`。独立 R 进程不继承 FaFile 文件描述符。

**结果**: ❌ `%over%` (GenomicRanges) 在 worker 中不可用。Worker 环境缺少必要的 S4 方法注册。

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
2. **`load_annotation` FaFile 兼容** (line 4780-4795): 不用 `is()` 检测，改为检查 `genome_package` 先再 fallback 到 `genome`
3. **`load()` 本地环境** (line 4776): 避免全局 FaFile 副本
4. **`mc.cleanup = FALSE`** (line 4472): 待测试

## 备选方案

- **SnowParam + 显式包加载**: Worker 启动时 `library(GenomicRanges)` 解决 `%over%` 不可用问题
- **n_cores=1**: 单线程处理，Nextflow 在 pipeline 层面并行 23 个样本
- **FaFile → TwoBitFile**: 使用 rtracklayer::TwoBitFile 替代，是否有相同问题需验证

## 文件改动位置

| 行号 | 改动 |
|------|------|
| 1791 | `exbin <- exonicParts(orfann, linked.to.single.gene.only = FALSE)` |
| 4170-4176 | 移除 FaFile 禁用（注释） |
| 4462-4473 | `mc.cleanup=FALSE, mc.silent=TRUE` |
| 4776 | `load(path, envir = _env)` 本地环境 |
| 4784-4792 | FaFile → DNAStringSet 预加载 + seqname 修复 |
| 5132 | `nsns <- exonicParts(annotation, linked.to.single.gene.only = FALSE)` |
