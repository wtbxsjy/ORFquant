# ORFquant Test Dataset

Minimal rice (Oryza sativa) dataset for debugging ORFquant multithreading issues.

## Files

| File | Description | Size |
|------|-------------|------|
| `mini_Rannot_fafile` | Annotation with FaFile genome (chr1 only, 50 transcripts) | ~11 MB |
| `mini_Rannot_dnastring` | Annotation with DNAStringSet genome (chr1 only, 50 transcripts) | ~11 MB |
| `mini_for_ORFquant` | P-site data (chr1 only, 10000 P-sites, 5000 unique) | ~1 MB |
| `chr1_data.rds` | Chr1 DNAStringSet sequence (43 Mb) | ~43 MB |

## Usage

```r
library(ORFquant)

# Test with FaFile version (reproduces finalizer bug with n_cores > 1)
run_ORFquant(
    for_ORFquant_file = "mini_for_ORFquant",
    annotation_file = "mini_Rannot_fafile",
    n_cores = 4,
    prefix = "test_mini",
    write_temp_files = TRUE,
    write_GTF_file = TRUE,
    write_protein_fasta = TRUE,
    interactive = FALSE
)

# Test with DNAStringSet version (should be fork-safe)
run_ORFquant(
    for_ORFquant_file = "mini_for_ORFquant",
    annotation_file = "mini_Rannot_dnastring",
    n_cores = 4,
    prefix = "test_mini_dna",
    write_temp_files = TRUE,
    write_GTF_file = TRUE,
    write_protein_fasta = TRUE,
    interactive = FALSE
)
```

## Source

- Full dataset: `~/riboseq/run/rice/process/work/`
- Chr1 size: 43,270,923 bp
- Full transcripts on chr1: 6,387
- Subset: 50 transcripts, 10,000 P-sites
