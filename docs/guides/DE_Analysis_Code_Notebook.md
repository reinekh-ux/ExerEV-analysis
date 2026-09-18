# Differential Expression Prep & Analysis — Code Notebook (Bone Marrow miRNA-EV)

A complete, explained record of the code used to go from the miRNA count matrix
to differential expression: exploratory variance analysis, variance partition,
and the DESeq2 model. Matches the earlier pipeline guides. All code runs in R on
Rorqual via the terminal (`Rscript`), using the cluster's R + Bioconductor.

**Load R environment (before any script):**
```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
```

**Inputs:**
- Count matrix: `count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt` (697 miRNAs × 62)
- Metadata: `metadata_BoneMarrow.txt` (62 samples; Condition, Sex, Isolation_Day, Seq_Batch, Treadmill_Batch)

**Model decided (evidence-based):** `~ Isolation_Day + Sex + Condition`, reference = Sed

---

## Why this order: explore first, then model

Rather than *assuming* which variables to correct for, the workflow first
measures which variables actually drive variation, then builds the model from
that evidence:

```
1. Exploratory PCA  ──► which variables drive the top principal components?
2. Variance partition ──► per-gene: how much variance does each variable explain?
3. Decide model ──► include the real drivers as covariates
4. DESeq2 ──► differential expression, adjusting for those covariates
```

---

# PART 1 — Exploratory PCA (`explore_PCA_BoneMarrow.R`)

**Purpose:** Model-free look at the data. Colour a PCA by each variable and
quantify (R²) how much variance each explains on the top PCs. **No DE model is
fitted** (`design = ~1`).

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

# Load counts + metadata, align, factorize
raw_counts <- read.table(rawcount_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
raw_counts <- round(raw_counts)
metadata <- read.table(metadata_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(metadata) <- metadata[[sample_id_col]]
metadata <- metadata[match(colnames(raw_counts), rownames(metadata)), ]
stopifnot(identical(colnames(raw_counts), rownames(metadata)))
for (v in inspect_vars) if (v %in% colnames(metadata)) metadata[[v]] <- factor(metadata[[v]])
metadata[[group_col]] <- relevel(factor(metadata[[group_col]]), ref=reference_level)

# Group-aware gene filter (same threshold as DE)
group_vec <- metadata[[group_col]]
passes <- raw_counts >= min_gene_count
max_pass <- apply(passes, 1, function(r) max(tapply(r, group_vec, sum), na.rm=TRUE))
counts_filt <- raw_counts[max_pass >= min_samples_in_group, , drop=FALSE]

# VST (model-free), then PCA on top-500 variable miRNAs
dds <- DESeqDataSetFromMatrix(counts_filt, metadata, design = ~ 1)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsmat <- assay(vsd)
ntop <- min(500, nrow(vsmat))
rv <- apply(vsmat, 1, var)
select <- order(rv, decreasing=TRUE)[seq_len(ntop)]
pca <- prcomp(t(vsmat[select, ]))
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# ... (PCA plots coloured by each variable saved as PNG) ...

# QUANTITATIVE: for each PC, regress PC scores on each variable -> R^2
pc_use <- min(n_pcs, ncol(pca$x))
r2 <- matrix(NA, length(inspect_vars), pc_use,
             dimnames=list(inspect_vars, paste0("PC", seq_len(pc_use))))
for (v in inspect_vars) for (k in seq_len(pc_use))
  r2[v,k] <- summary(lm(pca$x[,k] ~ metadata[[v]]))$r.squared
print(round(100*r2, 1))
```

**Key idea:** `lm(PC_score ~ variable)` gives an R² = the fraction of that PC's
spread explained by the variable. Higher R² = bigger driver.

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

**Interpretation:** Isolation_Day is the largest variance driver (35% PC1) →
technical batch to correct for. Condition (biology) is real but smaller (10–16%)
→ partly masked by isolation day, so adjusting helps. Seq_Batch overlaps
isolation day (nested). Sex/Treadmill small. Low per-PC variance = effects are
subtle and spread across many PCs (normal for EV small RNA).

---

# PART 2 — Variance Partition (`variancePartition_BoneMarrow.R`)

**Purpose:** The per-gene version of Part 1. For **every miRNA**, compute how
much variance each variable explains, and show the distribution as violins. This
is the rigorous, publication-standard variance decomposition.

```r
set.seed(910)
suppressMessages({ library(variancePartition); library(edgeR); library(DESeq2) })

working_dir <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
setwd(working_dir)
# ... load counts + metadata, factorize, same group-aware filter as Part 1 ...

# Normalize for varPart: TMM + log-CPM
dge <- DGEList(counts_filt); dge <- calcNormFactors(dge, "TMM")
logCPM <- cpm(dge, log=TRUE, prior.count=3)

# Each categorical variable as a random effect (1|var)
form <- ~ (1|Condition) + (1|Sex) + (1|Isolation_Day) + (1|Seq_Batch) + (1|Treadmill_Batch)
varPart <- fitExtractVarPartModel(logCPM, form, metadata)
vp <- sortCols(varPart)

plotVarPart(vp)                       # the violin plot
print(round(100*apply(vp,2,median),2))# median % variance per variable
```

**Key idea:** `fitExtractVarPartModel` fits a mixed model per miRNA and extracts
the % of variance attributable to each variable. `(1|var)` treats each variable
as a random effect (appropriate for categorical batch-like variables). The
violin shows the spread across all miRNAs.

### Result (this dataset)

- **Isolation_Day** — largest identifiable component; tail to ~33% for a subset
  of miRNAs (most near 0).
- **Condition** — real biological tail to ~30% for exercise-responsive miRNAs;
  most miRNAs unaffected (expected — exercise changes specific miRNAs).
- **Seq_Batch** — tail to ~27%, overlaps isolation day (nested).
- **Sex, Treadmill_Batch** — small.
- **Residuals** — dominant (median ~95%): most variance is gene-specific noise +
  individual variation. Normal for EV small RNA.

**Interpretation:** Confirms the model per-gene. Isolation day is the top
technical driver; exercise has genuine signal in a responsive subset; high
residuals mean expect modest, specific hits (not thousands), and covariate
correction matters. Because some variables are correlated (isolation day/seq
batch nested; treadmill/condition confounded), variance can be slightly split
between them — variancePartition handles this better than most methods, but read
those correlated variables cautiously.

---

# PART 3 — DESeq2 Differential Expression (`DEG_miRNA_BoneMarrow.R`)

**Purpose:** Test for exercise-responsive miRNAs, adjusting for the covariates
identified above. Model: `~ Isolation_Day + Sex + Condition`, reference = Sed.

### Settings block

```r
working_dir   <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
rawcount_file <- "count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt"
metadata_file <- "metadata_BoneMarrow.txt"
sample_id_col <- "Index_Number"
group_col       <- "Condition"
reference_level <- "Sed"
covariate_col   <- "Isolation_Day"   # technical covariate (evidence-based)
include_sex          <- TRUE          # sex as covariate
test_sex_interaction <- FALSE         # set TRUE later for sex-specific effects
min_sample_total_count <- 3
min_gene_count       <- 10
min_samples_in_group <- 3
padj_cutoff <- 0.05; lfc_cutoff <- 1.0
out_prefix  <- "BoneMarrow_miRNA"
```

### The pipeline

```r
# 1. Load counts + metadata (round counts to integers for DESeq2)
# 2. Filter: drop failed samples (safety net) + group-aware gene filter
#    (keep miRNAs with >=10 counts in >=3 samples of at least one condition)
# 3. Normalize two ways, saved to file:
#      - edgeR TMM  -> *_TMM_normalized_counts.txt
#      - DESeq2 median-of-ratios -> *_DESeq2_normalized_counts.txt
# 4. Build design formula: covariate + Sex + Condition (+ Sex:Condition if toggled)
# 5. VST + PCA (coloured by all 5 variables) + sample-distance heatmap (QC)
# 6. dds <- DESeq(dds); pull Condition_* contrasts (each timepoint vs Sed)
# 7. For each contrast: results table + EnhancedVolcano plot; count padj<0.05
```

### The model, explained

- **Condition** — the tested variable (Exer/Exer3/Exer7/Exer18 vs Sed).
- **Sex** — biological covariate; balanced across conditions, so clean. Absorbs
  male/female variation.
- **Isolation_Day** — technical covariate; the largest variance driver (Parts
  1–2). Also captures the RNA-extraction batch (done same day as isolation).
  Subsumes the nested sequencing-batch split, so Seq_Batch is not added.
- **Excluded:** Seq_Batch (nested in isolation day → redundant); Treadmill_Batch
  (confounded with condition → unidentifiable; used only for PCA colouring).

### Reference level & comparing all groups

DESeq2 fits one model; with Sed as reference the default contrasts are each
exercise timepoint vs Sed. Any other pairwise comparison is available after
fitting, e.g.:
```r
results(dds, contrast = c("Condition","Exer7","Exer"))   # Exer7 vs Exer
```

### Later — sex-specific effects

Because sex is balanced, the interaction model is supported. Flip
`test_sex_interaction <- TRUE` to add `Sex:Condition`, which tests whether the
exercise response differs by sex. Alternatively, run the pipeline separately per
sex.

### Run it

```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
# check EnhancedVolcano is available first:
Rscript -e 'library(EnhancedVolcano); cat("ok\n")'
Rscript DEG_miRNA_BoneMarrow.R
```

### Secondary analyses (tRNA / piRNA)

Point the same script at the other biotype matrices by changing two settings:
```r
rawcount_file <- "count_matrices/BoneMarrow_tRNA_TotalRawCount_matrix.txt"
out_prefix    <- "BoneMarrow_tRNA"
```
(and likewise for piRNA). Same model, same logic — gives exercise-responsive
tRNA fragments and piRNAs alongside the miRNA results.

---

# Outputs summary

| File | Content |
|------|---------|
| `*_EXPLORE_PCA_<var>.png` | PCA coloured by each variable |
| `*_EXPLORE_variance_explained_byPC.txt` | R² table (variable × PC) |
| `*_varPart_violin.png` | variance partition violin plot |
| `*_varPart_variance_table.txt` | per-miRNA variance partition |
| `*_TMM_normalized_counts.txt` | edgeR TMM normalized counts |
| `*_DESeq2_normalized_counts.txt` | DESeq2 normalized counts |
| `*_PCA_*.png`, `*_sample_distance_heatmap.png` | DE-script QC plots |
| `*_DESeq2_Condition_*.txt` | DE results per contrast |
| `*_Volcano_Condition_*.png` | volcano per contrast |

---

# Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| "unexpected symbol in module load" | typed shell command inside R | `quit()` back to `$` prompt first |
| nano paste saves empty file | large paste dropped over SSH | use `cat > file << 'EOF'` heredoc |
| package not found in RStudio | JupyterHub R differs from terminal R | run in terminal, or `BiocManager::install()` |
| EnhancedVolcano missing | not in Bioconductor bundle | install, or run DE in RStudio |
| varPart variance split oddly | correlated/confounded variables | read correlated vars cautiously |

---

*Reference for the bone marrow miRNA-EV project. Workflow: explore variance →
confirm drivers → model `~ Isolation_Day + Sex + Condition` (Sed reference) →
DESeq2. Evidence-based covariate selection via PCA R² and variancePartition.*
