# Exploratory Variance Analysis — Code Notebook (Bone Marrow miRNA-EV)

An explained record of the code used to explore sources of variation in the
miRNA count matrix **before** differential expression: exploratory PCA and
variance partition. The goal is evidence-based covariate selection — measure
which variables actually drive variation, then choose the DE model from that
evidence rather than assuming it.

All code runs in R on Rorqual via the terminal (`Rscript`), using the cluster's
R + Bioconductor.

**Load R environment (before any script):**
```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
```

**Inputs:**
- Count matrix: `count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt` (697 miRNAs × 62)
- Metadata: `metadata_BoneMarrow.txt` (62 samples; Condition, Sex, Isolation_Day, Seq_Batch, Treadmill_Batch)

---

## The workflow: explore first, then model

Rather than *assuming* which variables to correct for, the workflow first
measures which variables drive variation, then builds the model from that
evidence:

```
1. Exploratory PCA    ──► which variables drive the top principal components?
2. Variance partition ──► per-gene: how much variance does each variable explain?
3. Decide the model   ──► include the real drivers as covariates
   (DESeq2 differential expression — separate step, done after this)
```

---

## A note on the variables and how they relate

| Variable | Type | Levels | Note |
|----------|------|--------|------|
| Condition | biological | Sed, Exer, Exer3, Exer7, Exer18 | the effect of interest |
| Sex | biological | Male, Female | balanced across conditions |
| Isolation_Day | technical | 1–6 | EV isolation **and** RNA extraction done same day |
| Seq_Batch | technical | 1, 2 | two sequencing pools |
| Treadmill_Batch | technical | B, C | animal experiment run in 2 batches |

> **Isolation_Day covers extraction too:** EV isolation and RNA extraction were
> performed together on the same day, so `Isolation_Day` captures the combined
> technical variation of both wet-lab steps.

> **Isolation_Day is nested within Seq_Batch.** Isolation days 1–3 were
> sequenced in batch 1, and days 4–6 in batch 2. So each isolation day belongs
> entirely to one sequencing batch — knowing the day determines the batch.
> Isolation day is the **finer** variable; sequencing batch is a coarser
> grouping of it. This matters for modelling: adjusting for isolation day
> automatically adjusts for sequencing batch (each batch is just 3 of the days),
> so the two are not used together (that would be rank-deficient). Isolation day
> is preferred because it also captures finer day-to-day variation within each
> batch.

---

# PART 1 — Exploratory PCA (`explore_PCA_BoneMarrow.R`)

**Purpose:** Model-free look at the data. Colour a PCA by each variable and
quantify (R²) how much variance each explains on the top PCs. **No model is
fitted** (`design = ~1`), so the view is unbiased.

```r
set.seed(910)
suppressMessages({
  library(DESeq2); library(edgeR); library(ggplot2)
  library(RColorBrewer); library(pheatmap)
})

working_dir <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
setwd(working_dir)
rawcount_file <- "count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt"
metadata_file <- "metadata_BoneMarrow.txt"
sample_id_col <- "Index_Number"
group_col     <- "Condition"
reference_level <- "Sed"
inspect_vars <- c("Condition","Sex","Isolation_Day","Seq_Batch","Treadmill_Batch")
min_gene_count <- 10; min_samples_in_group <- 3; n_pcs <- 5
out_prefix <- "BoneMarrow_miRNA_EXPLORE"

# Load counts + metadata, align order, factorize
raw_counts <- read.table(rawcount_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
raw_counts <- round(raw_counts)
metadata <- read.table(metadata_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(metadata) <- metadata[[sample_id_col]]
metadata <- metadata[match(colnames(raw_counts), rownames(metadata)), ]
stopifnot(identical(colnames(raw_counts), rownames(metadata)))
for (v in inspect_vars) if (v %in% colnames(metadata)) metadata[[v]] <- factor(metadata[[v]])
metadata[[group_col]] <- relevel(factor(metadata[[group_col]]), ref=reference_level)

# Group-aware gene filter (keep miRNAs expressed in >=3 samples of any condition)
group_vec <- metadata[[group_col]]
passes <- raw_counts >= min_gene_count
max_pass <- apply(passes, 1, function(r) max(tapply(r, group_vec, sum), na.rm=TRUE))
counts_filt <- raw_counts[max_pass >= min_samples_in_group, , drop=FALSE]

# VST (model-free) then PCA on the top-500 most variable miRNAs
dds <- DESeqDataSetFromMatrix(counts_filt, metadata, design = ~ 1)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsmat <- assay(vsd)
ntop <- min(500, nrow(vsmat))
rv <- apply(vsmat, 1, var)
select <- order(rv, decreasing=TRUE)[seq_len(ntop)]
pca <- prcomp(t(vsmat[select, ]))
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# (PCA scatter coloured by each variable saved as PNG — one per variable)

# QUANTITATIVE: for each PC, regress PC scores on each variable -> R^2
pc_use <- min(n_pcs, ncol(pca$x))
r2 <- matrix(NA, length(inspect_vars), pc_use,
             dimnames=list(inspect_vars, paste0("PC", seq_len(pc_use))))
for (v in inspect_vars) for (k in seq_len(pc_use))
  r2[v,k] <- summary(lm(pca$x[,k] ~ metadata[[v]]))$r.squared
print(round(100*r2, 1))
```

**Key idea:** `lm(PC_score ~ variable)` gives R² = the fraction of that PC's
spread explained by the variable. Higher R² = bigger driver. Doing this per PC
shows which variables dominate the main axes of variation.

The script also saves a **sample-distance heatmap** (unsupervised clustering
with all variables annotated) to check for outliers.

### Result (this dataset)

```
                 PC1   PC2   PC3   PC4   PC5
Condition       10.3  15.7   1.5   5.1   3.2
Sex              2.0   3.9   5.6   0.0   0.2
Isolation_Day   35.0  24.9  14.3  32.6  15.4   <- dominant
Seq_Batch       24.0  13.9   0.6   9.3   3.8
Treadmill_Batch  0.1  11.8   0.0   3.0   0.2
```
% variance per PC: PC1 8.7, PC2 7.0, PC3 5.6, PC4 4.4, PC5 3.9

**Interpretation:**
- **Isolation_Day** is the largest variance driver (35% on PC1) → the main
  technical batch effect to correct for.
- **Condition** (biology) is real but smaller (10–16%) → partly masked by
  isolation day, so adjusting for it should help the exercise signal emerge.
- **Seq_Batch** (24% PC1) overlaps isolation day because it is nested within it
  → captured by isolation day, not used separately.
- **Sex** and **Treadmill_Batch** are small.
- **Low per-PC variance** (PC1 only 8.7%) → effects are subtle and spread across
  many PCs. Normal for EV small RNA. The sample-distance heatmap was uniform
  with no outliers → all 62 samples retained.

---

# PART 2 — Variance Partition (`variancePartition_BoneMarrow.R`)

**Purpose:** The per-gene version of Part 1. For **every miRNA**, compute how
much variance each variable explains, and show the distribution as violins. This
is the publication-standard variance decomposition (the `variancePartition`
Bioconductor package).

```r
set.seed(910)
suppressMessages({ library(variancePartition); library(edgeR); library(DESeq2) })

working_dir <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
setwd(working_dir)
rawcount_file <- "count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt"
metadata_file <- "metadata_BoneMarrow.txt"
sample_id_col <- "Index_Number"; group_col <- "Condition"
min_gene_count <- 10; min_samples_in_group <- 3
out_prefix <- "BoneMarrow_miRNA_varPart"

# Load counts + metadata, factorize, same group-aware filter as Part 1
raw_counts <- read.table(rawcount_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
raw_counts <- round(raw_counts)
metadata <- read.table(metadata_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(metadata) <- metadata[[sample_id_col]]
metadata <- metadata[match(colnames(raw_counts), rownames(metadata)), ]
stopifnot(identical(colnames(raw_counts), rownames(metadata)))
for (v in c("Condition","Sex","Isolation_Day","Seq_Batch","Treadmill_Batch"))
  if (v %in% colnames(metadata)) metadata[[v]] <- factor(metadata[[v]])
gv <- metadata[[group_col]]
passes <- raw_counts >= min_gene_count
max_pass <- apply(passes, 1, function(r) max(tapply(r, gv, sum), na.rm=TRUE))
counts_filt <- raw_counts[max_pass >= min_samples_in_group, , drop=FALSE]

# Normalize for varPart: TMM + log-CPM
dge <- DGEList(counts_filt); dge <- calcNormFactors(dge, "TMM")
logCPM <- cpm(dge, log=TRUE, prior.count=3)

# Each categorical variable as a random effect (1|var)
form <- ~ (1|Condition) + (1|Sex) + (1|Isolation_Day) + (1|Seq_Batch) + (1|Treadmill_Batch)
varPart <- fitExtractVarPartModel(logCPM, form, metadata)
vp <- sortCols(varPart)

plotVarPart(vp)                          # the violin plot
print(round(100*apply(vp,2,median),2))   # median % variance per variable
write.table(as.data.frame(vp), paste0(out_prefix,"_variance_table.txt"),
            sep="\t", quote=FALSE, col.names=NA)
```

**Key idea:** `fitExtractVarPartModel` fits a mixed model per miRNA and extracts
the % of variance attributable to each variable. `(1|var)` treats each variable
as a random effect (appropriate for categorical batch-like variables). Each dot
in the violin is one miRNA; the violin shows the distribution across all miRNAs.

### Result (this dataset)

- **Isolation_Day** — the largest identifiable component; tail to ~33% for a
  subset of miRNAs (most near 0).
- **Condition** — real biological tail to ~30% for exercise-responsive miRNAs;
  most miRNAs unaffected (expected — exercise changes specific miRNAs, not all).
- **Seq_Batch** — tail to ~27%, overlaps isolation day (nested within it).
- **Sex, Treadmill_Batch** — small.
- **Residuals** — dominant (median ~95%): most variance is gene-specific noise +
  individual animal variation. Normal for EV small RNA.

**Interpretation:** Confirms Part 1 per-gene. Isolation day is the top
identifiable technical driver; exercise has genuine signal in a responsive
subset of miRNAs; high residuals mean expect modest, specific hits (not
thousands), and covariate correction matters.

> **Caveat on correlated variables:** because some variables are related
> (isolation day nested within seq batch; treadmill batch confounded with
> condition), variance can be split between them. `variancePartition` handles
> correlated predictors better than most methods, but the correlated variables
> (Seq_Batch, Treadmill_Batch) should be read cautiously — the trustworthy
> drivers here are Isolation_Day and Condition.

---

## What this exploration establishes (for the DE model)

The two analyses agree, at the PC level and per gene:

| Variable | Verdict | Action in DE model |
|----------|---------|--------------------|
| Isolation_Day | largest technical driver | **include** as covariate |
| Condition | real biological signal (subset of miRNAs) | **tested** (reference = Sed) |
| Sex | small but present; balanced | **include** as covariate (cheap, clean) |
| Seq_Batch | nested within isolation day | **exclude** (redundant) |
| Treadmill_Batch | confounded with condition; small | **exclude** (PCA colour only) |

**Resulting model (built in the separate DESeq2 step):**
`~ Isolation_Day + Sex + Condition`, reference = Sed — chosen from this evidence,
not assumed.

---

# Outputs

| File | Content |
|------|---------|
| `*_EXPLORE_PCA_<var>.png` | PCA coloured by each variable |
| `*_EXPLORE_variance_explained_byPC.txt` | R² table (variable × PC) |
| `*_EXPLORE_sample_distance_heatmap.png` | unsupervised clustering QC |
| `*_varPart_violin.png` | variance partition violin plot |
| `*_varPart_variance_table.txt` | per-miRNA variance partition |

---

# Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| "unexpected symbol in module load" | typed a shell command inside R | `quit()` back to the `$` prompt first |
| nano paste saves an empty file | large paste dropped over SSH | use `cat > file << 'EOF'` heredoc instead |
| package not found in RStudio | JupyterHub R differs from terminal R | run in terminal, or `BiocManager::install()` |
| varPart variance split oddly | correlated / confounded variables | read correlated vars cautiously |

---

*Reference for the bone marrow miRNA-EV project. Exploratory step: measure
variance drivers (PCA R² + variancePartition) to justify covariate choice before
differential expression. Isolation day (nested-parent of sequencing batch, and
covering RNA extraction) is the dominant technical component; condition carries
real signal in a responsive subset of miRNAs.*
