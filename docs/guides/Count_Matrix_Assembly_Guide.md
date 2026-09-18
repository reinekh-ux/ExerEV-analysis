# Assembling Count Matrices from exceRpt Output — Step-by-Step Guide

After exceRpt quantifies each of the 62 samples individually, this step
**combines** the per-sample results into count matrices — one row per small
RNA feature, one column per sample — ready for differential expression
(DESeq2). Written to match the QC / trimming / exceRpt guides.

**Input:** 62 exceRpt output folders (`exceRpt_output_mm10/<sample>/`)
**Output:** one count matrix per small RNA biotype (miRNA, tRNA, piRNA, gencode)
**Tool:** R (cluster module — no container needed for this step)

---

## Why this step is needed

exceRpt processes each sample on its own and writes a separate set of count
files per sample. For analysis you need a single table where you can compare
the same feature across all samples — a matrix:

```
                D1_1   D1_2   D1_3   ...   D6_12
mmu-miR-486a    15495  17442  91674  ...   75177
mmu-miR-92a     10080  10235  67518  ...   ...
...
```

This step walks every sample folder, pulls the count for each feature, and
stitches them together into that matrix.

> **Why not the container's `mergePipelineRuns.R`?** exceRpt ships a merge
> script, but it requires R packages (plyr, gplots, marray, Rgraphviz…) that
> failed to load when run on the cluster (the host environment overrode the
> container's R library path). The lab's solution — used for the muscle
> dataset too — is a small custom R script that runs with the cluster's own
> R + Bioconductor module, which already has `dplyr`. Simpler and robust.

---

## What counts are used

| Choice | Value | Why |
|--------|-------|-----|
| Count file per sample | `readCounts_<biotype>_sense.txt` | sense strand = real signal |
| Count column | `totalReadCount` | total (not unique) reads, matching muscle |
| Feature ID column | `ReferenceID` | the miRNA / RNA identifier |
| Missing values | set to `0` | a feature absent in a sample = 0 counts |

> These are **raw counts** — not normalized. That is correct: DESeq2 performs
> its own normalization (median-of-ratios) and requires raw integer counts as
> input. Normalizing beforehand would break its model.

---

## The biotypes assembled

exceRpt quantifies several small RNA classes. This step builds a matrix for
each of the main ones:

| Biotype | Count file | Notes |
|---------|-----------|-------|
| **miRNA** (mature) | `readCounts_miRNAmature_sense.txt` | primary analysis target |
| **tRNA** | `readCounts_tRNA_sense.txt` | tRNA fragments — abundant in EVs |
| **piRNA** | `readCounts_piRNA_sense.txt` | piwi-interacting RNAs |
| **gencode** | `readCounts_gencode_sense.txt` | other transcripts (mRNA frag, lncRNA) |

(circRNA was near-zero in this dataset and miRNA-precursor is not mature
miRNA, so both were left out — they can be added with one line each.)

---

## A key detail — the two sub-folders

Each sample's exceRpt output folder contains **two** result sub-folders:

```
exceRpt_output_mm10/D1_1/
├── ..._D1_1/                          ← main results (use this one)
├── ..._D1_1_CORE_RESULTS_v4.6.3/      ← packaged copy (ignore)
├── ..._D1_1_CORE_RESULTS_v4.6.3.tgz
├── ..._D1_1.stats
├── ..._D1_1.log
├── ..._D1_1.err
└── ..._D1_1.qcResult
```

The muscle data had only one sub-folder, so the original script (which
expected exactly one) skipped every sample. The fix: list the sub-folders and
pick the one **without** `CORE_RESULTS` in its name:

```r
main_sub <- subfolders[!grepl("CORE_RESULTS", subfolders)]
```

That reliably selects the main results folder where the count files live.

---

# PART 1 — The R assembly script

This version builds a matrix for every biotype in one run.

```bash
cat > ~/links/scratch/BoneMarrow_miRNA_EV_QC/assemble_all_biotypes.R << 'ENDOFRSCRIPT'
set.seed(910)
suppressMessages(library(dplyr))

working_dir <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
input_dir   <- "exceRpt_output_mm10"
count_column <- "totalReadCount"

# biotype -> per-sample count file -> output matrix name
biotypes <- list(
  miRNA   = list(file="readCounts_miRNAmature_sense.txt", out="BoneMarrow_miRNA_TotalRawCount_matrix.txt"),
  tRNA    = list(file="readCounts_tRNA_sense.txt",        out="BoneMarrow_tRNA_TotalRawCount_matrix.txt"),
  piRNA   = list(file="readCounts_piRNA_sense.txt",       out="BoneMarrow_piRNA_TotalRawCount_matrix.txt"),
  gencode = list(file="readCounts_gencode_sense.txt",     out="BoneMarrow_gencode_TotalRawCount_matrix.txt")
)

setwd(working_dir)
sample_folders <- list.dirs(input_dir, recursive = FALSE, full.names = FALSE)
message("Found ", length(sample_folders), " sample folders")

for (bt in names(biotypes)) {
  count_filename <- biotypes[[bt]]$file
  output_file    <- biotypes[[bt]]$out
  message("=== Building ", bt, " from ", count_filename, " ===")
  all_tabs <- list(); skipped <- c()
  for (folder in sample_folders) {
    subfolders <- list.dirs(file.path(input_dir, folder), recursive = FALSE, full.names = FALSE)
    main_sub   <- subfolders[!grepl("CORE_RESULTS", subfolders)]
    if (length(main_sub) != 1) { skipped <- c(skipped, folder); next }
    count_file <- file.path(input_dir, folder, main_sub, count_filename)
    if (!file.exists(count_file)) { skipped <- c(skipped, folder); next }
    tbl <- tryCatch(read.table(count_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char=""), error=function(e) NULL)
    if (is.null(tbl) || nrow(tbl)==0 || !all(c("ReferenceID", count_column) %in% colnames(tbl))) { skipped <- c(skipped, folder); next }
    tbl <- tbl %>% select(all_of(c("ReferenceID", count_column)))
    colnames(tbl) <- c("Gene", folder)
    all_tabs[[folder]] <- tbl
  }
  if (length(all_tabs)==0) { message("  no usable tables for ", bt); next }
  mat <- Reduce(function(x,y) full_join(x,y,by="Gene"), all_tabs)
  mat[is.na(mat)] <- 0
  write.table(mat, file=output_file, row.names=FALSE, quote=FALSE, sep="\t")
  message("  ", bt, ": ", nrow(mat), " features x ", ncol(mat)-1, " samples -> ", output_file)
  if (length(skipped)) message("  skipped: ", paste(skipped, collapse=", "))
}
message("Done.")
ENDOFRSCRIPT
```

### How the script works, step by step

1. **Settings** — the working directory, the exceRpt output folder, the count
   column (`totalReadCount`), and a list of biotypes with their count files and
   output names.
2. **Find samples** — `list.dirs(input_dir)` lists all 62 sample folders.
3. **For each biotype**, loop over every sample:
   - find the main result sub-folder (excluding `CORE_RESULTS`)
   - read that sample's count file
   - keep two columns: the feature ID (`ReferenceID`) and the count
     (`totalReadCount`), renaming the count column to the sample name
4. **Merge** — `Reduce(full_join, ...)` joins all 62 per-sample tables by
   feature into one matrix. `full_join` keeps the union of all features across
   samples; any feature missing in a sample becomes `NA`, then set to `0`.
5. **Write** — saves each biotype matrix as a tab-separated `.txt`.

> `full_join` is the right choice because not every feature appears in every
> sample — the union ensures no feature is lost, and zeros fill the gaps.

---

# PART 2 — The SLURM submission script

```bash
cat > ~/links/scratch/BoneMarrow_miRNA_EV_QC/assemble_all_biotypes.sbatch << 'ENDOFSBATCH'
#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --job-name="biotype_assembly"
#SBATCH --output=log/slurm-%x.%j.out
#SBATCH --error=log/slurm-%x.%j.err
#SBATCH --mem=50G
#SBATCH --time=3:00:00
#SBATCH --mail-user=reine.khoury@mail.mcgill.ca
#SBATCH --mail-type=ALL

cd "${SLURM_SUBMIT_DIR}"
mkdir -p log
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
Rscript assemble_all_biotypes.R
ENDOFSBATCH
```

> The `r-bundle-bioconductor/3.16` module provides `dplyr` and other R packages
> — this is why no container is needed for the assembly step.

---

# PART 3 — Run it

```bash
cd ~/links/scratch/BoneMarrow_miRNA_EV_QC
sbatch assemble_all_biotypes.sbatch
squeue -u $USER
```

It runs in a few minutes. The summary messages (how many features × samples,
and any skipped samples) go to the **`.err`** log (R's `message()` writes to
stderr), not the `.out` log:

```bash
cat ~/links/scratch/BoneMarrow_miRNA_EV_QC/log/slurm-biotype_assembly.*.err
```

Check the matrices and their dimensions:

```bash
for f in ~/links/scratch/BoneMarrow_miRNA_EV_QC/Bone*_TotalRawCount_matrix.txt; do
  echo "$(basename $f): $(($(wc -l < $f) - 1)) features x $(($(head -1 $f | wc -w) - 1)) samples"
done
```

Each should show **62 samples**.

---

# PART 4 — Results (this dataset)

| Biotype | Features | Samples |
|---------|---------:|--------:|
| miRNA | 697 | 62 |
| tRNA | (per log) | 62 |
| piRNA | (per log) | 62 |
| gencode | (per log) | 62 |

For reference, the muscle dataset detected **468 miRNAs**; bone marrow's 697
is in the same order of magnitude (higher is reasonable — bone marrow is a
hematopoietic tissue with a rich miRNA repertoire), confirming the
quantification is sound.

The matrix format:

```
Gene                                    D1_1   D1_2   ...
mmu-miR-486a-5p:MIMAT0003130:...        15495  17442  ...
mmu-miR-92a-3p:MIMAT0000539:...         10080  10235  ...
```

First column = feature ID; remaining columns = raw counts per sample.

---

# PART 5 — Organize and back up

Scratch is purged, so move the matrices into a folder and copy to project:

```bash
# folder in scratch
mkdir -p ~/links/scratch/BoneMarrow_miRNA_EV_QC/count_matrices
mv ~/links/scratch/BoneMarrow_miRNA_EV_QC/Bone*_TotalRawCount_matrix.txt \
   ~/links/scratch/BoneMarrow_miRNA_EV_QC/count_matrices/

# folder in project + copy
mkdir -p ~/links/projects/rrg-gturecki/reinekh/BoneMarrow_miRNA_EV/count_matrices
cp ~/links/scratch/BoneMarrow_miRNA_EV_QC/count_matrices/*.txt \
   ~/links/projects/rrg-gturecki/reinekh/BoneMarrow_miRNA_EV/count_matrices/

# back up the scripts too
cp ~/links/scratch/BoneMarrow_miRNA_EV_QC/assemble_*.* \
   ~/links/projects/rrg-gturecki/reinekh/BoneMarrow_miRNA_EV/scripts/
```

---

# PART 6 — Next step: DESeq2

Each raw count matrix feeds into DESeq2 in R for differential expression
(sedentary vs exercise), with **batch** as a covariate (batch 1 = 298.001
files, batch 2 = 299.001 files). miRNA is the primary analysis; tRNA/piRNA can
be analyzed the same way as secondary EV cargo.

---

# Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| Every sample "skipped: found 2 sub-folders" | exceRpt writes a main + CORE_RESULTS folder | pick the non-CORE_RESULTS one with `!grepl("CORE_RESULTS", ...)` |
| "No count tables were read" | wrong count filename or wrong sub-folder | check `readCounts_<biotype>_sense.txt` exists |
| `.sbatch` won't submit / "Unable to open file" | heredoc block didn't close (prompt showed `>`) | create R and sbatch files separately; wait for `$` between |
| Summary messages missing from `.out` | R `message()` writes to stderr | read the `.err` log instead |
| `Bone*` wildcard "No such file" | matrices already moved | point the command at `count_matrices/` |

---

*Reference for the bone marrow miRNA-EV project. Assembly uses the cluster's
R + Bioconductor (no container). Matches the muscle dataset's approach. Output
is raw counts for DESeq2.*
