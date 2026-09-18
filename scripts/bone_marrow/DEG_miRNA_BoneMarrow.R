###############################################################################
# miRNA differential-expression analysis: filter -> normalize -> PCA -> DEG
#
# Project : Reine - Bone Marrow EV small-RNA sequencing (exercise time course)
# Adapted from the muscle DEG script (Minh's pipeline) for the bone marrow data.
#
# Pipeline:
#   1. Load raw miRNA count matrix + metadata
#   2. Filter low-count samples and lowly-expressed miRNAs (group-aware)
#   3. Normalize (TMM via edgeR AND median-of-ratios via DESeq2) -> saved
#   4. PCA / sample-similarity QC plots (coloured by all design variables)
#   5. DESeq2 differential expression vs reference (Sed)
#   6. Volcano plot + results table for every contrast
#
# Conditions: Sed, Exer, Exer3, Exer7, Exer18   (reference = Sed)
# Model     : ~ Isolation_Day + Sex + Condition
#
# How to run:
#   module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
#   Rscript DEG_miRNA_BoneMarrow.R
# (or run interactively in RStudio with the same libraries available)
###############################################################################

set.seed(910)

# ----------------------------------------------------------------------------
# Libraries
# ----------------------------------------------------------------------------
suppressMessages({
  library(edgeR)            # TMM normalization
  library(DESeq2)           # DE testing + variance-stabilizing transform
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
  library(EnhancedVolcano)  # volcano plots
})

# ============================================================================
# SETTINGS  --  edit this block, then run the whole script.
# ============================================================================
working_dir   <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
setwd(working_dir)

# Inputs
rawcount_file <- "count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt"
metadata_file <- "metadata_BoneMarrow.txt"

# Column in metadata that matches the COLUMN NAMES (sample IDs) of the matrix.
sample_id_col <- "Index_Number"

# Grouping variable to test, and the reference level all contrasts compare to.
group_col       <- "Condition"
reference_level <- "Sed"

# Technical covariate to adjust for in the DE model (set to NA to ignore).
# Isolation_Day is the recommended choice (it also tracks the sequencing batch,
# which is nested within it). Alternative: "Seq_Batch".
covariate_col <- "Isolation_Day"

# ---- Sex handling -----------------------------------------------------------
# include_sex   : add Sex as a covariate in the model (recommended = TRUE)
# test_sex_interaction : add Sex:Condition to test sex-specific exercise effects
#                        (set TRUE later for the sex-difference analysis)
include_sex           <- TRUE
test_sex_interaction  <- FALSE
# -----------------------------------------------------------------------------

# Filtering thresholds
min_sample_total_count <- 3   # drop a sample if its total counts < this (safety net)
min_gene_count       <- 10    # per-sample count threshold for the gene filter
min_samples_in_group <- 3     # how many samples in a group must pass it

# Significance cutoffs for the volcano plots
padj_cutoff <- 0.05
lfc_cutoff  <- 1.0

# Output prefix (all output files start with this)
out_prefix <- "BoneMarrow_miRNA"
# ============================================================================


# ----------------------------------------------------------------------------
# 1. Load raw counts + metadata
# ----------------------------------------------------------------------------
raw_counts <- read.table(rawcount_file, header = TRUE, row.names = 1,
                         sep = "\t", check.names = FALSE)
raw_counts <- round(raw_counts)          # DESeq2 needs integer counts
message("Raw count matrix: ", nrow(raw_counts), " miRNAs x ", ncol(raw_counts), " samples")

metadata <- read.table(metadata_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
rownames(metadata) <- metadata[[sample_id_col]]

# Sanity: every count column should have a metadata row.
missing <- setdiff(colnames(raw_counts), rownames(metadata))
if (length(missing)) stop("Samples in counts but not metadata: ", paste(missing, collapse=", "))


# ----------------------------------------------------------------------------
# 2. Filter low-count samples and lowly-expressed miRNAs
# ----------------------------------------------------------------------------
# 2a. Drop samples whose total library size is too small.
low_samples <- colnames(raw_counts)[colSums(raw_counts) < min_sample_total_count]
if (length(low_samples)) message("Dropping low-count samples: ", paste(low_samples, collapse = ", "))
counts_filt <- raw_counts[, colSums(raw_counts) >= min_sample_total_count, drop = FALSE]

# 2b. Align metadata to the (filtered) count columns and check the order.
metadata <- metadata[match(colnames(counts_filt), rownames(metadata)), ]
stopifnot(identical(colnames(counts_filt), rownames(metadata)))

# Factors. Group with chosen reference first; covariate/sex as factors.
metadata[[group_col]] <- relevel(factor(metadata[[group_col]]), ref = reference_level)
if (!is.na(covariate_col)) metadata[[covariate_col]] <- factor(metadata[[covariate_col]])
if (include_sex)           metadata[["Sex"]]         <- factor(metadata[["Sex"]])
# Treadmill_Batch kept as a factor for PCA colouring only (never in the model).
if ("Treadmill_Batch" %in% colnames(metadata))
  metadata[["Treadmill_Batch"]] <- factor(metadata[["Treadmill_Batch"]])
if ("Seq_Batch" %in% colnames(metadata))
  metadata[["Seq_Batch"]] <- factor(metadata[["Seq_Batch"]])

# 2c. Group-aware gene filter.
#     Keep a miRNA if, in at least one group, >= min_samples_in_group samples
#     have >= min_gene_count counts.
group_vec <- metadata[[group_col]]
passes_threshold <- counts_filt >= min_gene_count
max_pass_in_any_group <- apply(passes_threshold, 1, function(gene_row) {
  per_group <- tapply(gene_row, group_vec, sum)
  max(per_group, na.rm = TRUE)
})
keep <- max_pass_in_any_group >= min_samples_in_group

n_before <- nrow(counts_filt)
counts_filt <- counts_filt[keep, , drop = FALSE]
message("Gene filter: kept miRNAs with >= ", min_gene_count, " counts in >= ",
        min_samples_in_group, " samples of at least one group")
message("Dropped ", n_before - nrow(counts_filt), " lowly-expressed miRNAs")
message("After filtering: ", nrow(counts_filt), " miRNAs x ", ncol(counts_filt), " samples")


# ----------------------------------------------------------------------------
# 3. Normalization (saved for downstream / plotting use)
# ----------------------------------------------------------------------------
# 3a. edgeR TMM-normalized CPM
dge <- DGEList(counts_filt)
dge <- calcNormFactors(dge, method = "TMM")
tmm_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
write.table(tmm_counts, file = paste0(out_prefix, "_TMM_normalized_counts.txt"),
            quote = FALSE, sep = "\t", col.names = NA)

# 3b. DESeq2 median-of-ratios normalized counts
dds_tmp <- DESeqDataSetFromMatrix(countData = counts_filt,
                                  colData = metadata,
                                  design = as.formula(paste("~", group_col)))
dds_tmp <- estimateSizeFactors(dds_tmp)
write.table(counts(dds_tmp, normalized = TRUE),
            file = paste0(out_prefix, "_DESeq2_normalized_counts.txt"),
            quote = FALSE, sep = "\t", col.names = NA)


# ----------------------------------------------------------------------------
# 4. Build the design formula
# ----------------------------------------------------------------------------
# Terms: [covariate] + [Sex] + Condition  (+ Sex:Condition if requested)
terms <- c()
if (!is.na(covariate_col)) terms <- c(terms, covariate_col)
if (include_sex)           terms <- c(terms, "Sex")
terms <- c(terms, group_col)
rhs <- paste(terms, collapse = " + ")
if (test_sex_interaction && include_sex) {
  rhs <- paste0(rhs, " + Sex:", group_col)
}
design_formula <- as.formula(paste("~", rhs))
message("Design: ", deparse(design_formula))

dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                              colData = metadata,
                              design = design_formula)


# ----------------------------------------------------------------------------
# 4b. PCA + sample-similarity QC
# ----------------------------------------------------------------------------
# Variance-stabilizing transform for visualisation (not for DE testing).
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

save_pca <- function(vsd, intgroup, file) {
  if (!(intgroup %in% colnames(colData(vsd)))) return(invisible())
  p <- plotPCA(vsd, intgroup = intgroup) +
    theme_bw() +
    theme(axis.text  = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12))
  ggsave(file, p, width = 8, height = 6, dpi = 300)
  message("Saved PCA: ", file)
}

# Colour the PCA by each design variable to inspect what drives clustering.
save_pca(vsd, group_col,          paste0(out_prefix, "_PCA_Condition.png"))
save_pca(vsd, "Sex",              paste0(out_prefix, "_PCA_Sex.png"))
save_pca(vsd, "Isolation_Day",    paste0(out_prefix, "_PCA_IsolationDay.png"))
save_pca(vsd, "Seq_Batch",        paste0(out_prefix, "_PCA_SeqBatch.png"))
save_pca(vsd, "Treadmill_Batch",  paste0(out_prefix, "_PCA_TreadmillBatch.png"))

# Sample-to-sample distance heatmap.
sampleDists <- dist(t(assay(vsd)))
distMat <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
png(paste0(out_prefix, "_sample_distance_heatmap.png"), width = 9, height = 8,
    units = "in", res = 300)
pheatmap(distMat,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
message("Saved heatmap: ", out_prefix, "_sample_distance_heatmap.png")


# ----------------------------------------------------------------------------
# 5. Differential expression (DESeq2)
# ----------------------------------------------------------------------------
dds <- DESeq(dds)

# Main effect contrasts: each condition vs the reference (Sed).
# (resultsNames includes intercept + covariate levels; we keep only the
#  Condition_* contrasts so volcano/tables are the exercise effects.)
all_names <- resultsNames(dds)
message("All model coefficients: ", paste(all_names, collapse = "  |  "))
condition_contrasts <- grep(paste0("^", group_col, "_"), all_names, value = TRUE)
message("Condition contrasts: ", paste(condition_contrasts, collapse = "  |  "))


# ----------------------------------------------------------------------------
# 6. Results table + volcano plot for each contrast
# ----------------------------------------------------------------------------
make_volcano <- function(res, title, file) {
  ymax <- if (all(is.na(res$padj))) 1 else max(-log10(res$padj), na.rm = TRUE) + 1
  p <- EnhancedVolcano(res,
                       lab        = rownames(res),
                       x          = "log2FoldChange",
                       y          = "padj",
                       pCutoff    = padj_cutoff,
                       FCcutoff   = lfc_cutoff,
                       ylim       = c(0, ymax),
                       title      = title,
                       subtitle   = NULL,
                       labSize    = 3,
                       labFace    = "bold",
                       labCol     = "black",
                       pointSize  = 3,
                       ylab       = "-log10(padj)") +
    theme(axis.title = element_text(size = 16, face = "bold"),
          axis.text  = element_text(size = 14, face = "bold"))
  ggsave(file, p, width = 8, height = 6, dpi = 300)
  message("Saved volcano: ", file)
}

for (contrast in condition_contrasts) {
  res <- results(dds, name = contrast)
  res <- res[order(res$padj), ]

  tag <- gsub("[^A-Za-z0-9]+", "_", contrast)
  write.table(as.data.frame(res),
              file = paste0(out_prefix, "_DESeq2_", tag, ".txt"),
              sep = "\t", quote = FALSE, col.names = NA)

  make_volcano(res, title = contrast,
               file = paste0(out_prefix, "_Volcano_", tag, ".png"))

  n_sig <- sum(res$padj < padj_cutoff, na.rm = TRUE)
  message(contrast, ": ", n_sig, " miRNAs with padj < ", padj_cutoff)
}

message("=== DONE ===")
