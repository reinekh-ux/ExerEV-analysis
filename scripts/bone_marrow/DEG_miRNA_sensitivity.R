###############################################################################
# DESeq2 threshold sensitivity analysis — Bone Marrow miRNA-EV
#
# Runs the SAME model across several filtering / significance thresholds and
# reports how many real (baseMean-filtered) hits each yields. This is a
# robustness / sensitivity analysis, not threshold-hunting: the point is to
# show whether any signal is stable across stringency levels.
#
# All outputs are organized into a results folder with per-setting subfolders.
#
# How to run:
#   module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
#   Rscript DEG_miRNA_sensitivity.R
###############################################################################

set.seed(910)
suppressMessages({
  library(edgeR); library(DESeq2); library(dplyr)
  library(ggplot2); library(RColorBrewer)
})

# ============================================================================
# FIXED SETTINGS (same for all threshold sets)
# ============================================================================
working_dir   <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
setwd(working_dir)
rawcount_file <- "count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt"
metadata_file <- "metadata_BoneMarrow.txt"
sample_id_col <- "Index_Number"
group_col       <- "Condition"
reference_level <- "Sed"
covariate_col   <- "Isolation_Day"   # technical covariate (change to "Seq_Batch" if desired)
include_sex     <- TRUE
min_sample_total_count <- 3
min_baseMean_real      <- 10         # hits below this baseMean are low-count artifacts

# Root results folder (everything goes under here)
results_root <- "DESeq2_results"

# ============================================================================
# THRESHOLD SETS to compare  --  add/remove rows here
#   name              : label + subfolder name
#   min_gene_count    : per-sample count threshold for the gene filter
#   min_samples       : samples in a group that must pass (expression filter)
#   padj              : significance cutoff (FDR)
#   lfc               : fold-change mark for volcano (log2)
# ============================================================================
threshold_sets <- list(
  list(name="strict",         min_gene_count=10, min_samples=5, padj=0.05, lfc=1.0),
  list(name="standard",       min_gene_count=10, min_samples=3, padj=0.05, lfc=1.0),
  list(name="lenient",        min_gene_count=5,  min_samples=3, padj=0.10, lfc=0.5),
  list(name="very_lenient",   min_gene_count=5,  min_samples=2, padj=0.10, lfc=0.5)
)
# ============================================================================


# ---- Load raw data ONCE ----
raw_counts <- read.table(rawcount_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
raw_counts <- round(raw_counts)
message("Raw: ", nrow(raw_counts), " miRNAs x ", ncol(raw_counts), " samples")
metadata <- read.table(metadata_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(metadata) <- metadata[[sample_id_col]]

# Drop failed samples (safety net), align metadata
counts_all <- raw_counts[, colSums(raw_counts) >= min_sample_total_count, drop=FALSE]
metadata <- metadata[match(colnames(counts_all), rownames(metadata)), ]
stopifnot(identical(colnames(counts_all), rownames(metadata)))
metadata[[group_col]]     <- relevel(factor(metadata[[group_col]]), ref=reference_level)
metadata[[covariate_col]] <- factor(metadata[[covariate_col]])
if (include_sex) metadata[["Sex"]] <- factor(metadata[["Sex"]])

# Build design formula once
terms <- c(covariate_col, if (include_sex) "Sex", group_col)
design_formula <- as.formula(paste("~", paste(terms, collapse=" + ")))
message("Design: ", deparse(design_formula))

# ---- Create the results folder structure ----
dir.create(results_root, showWarnings=FALSE)
message("Results root: ", file.path(working_dir, results_root))

# Volcano helper
make_volcano <- function(res, title, file, padj_cut, lfc_cut) {
  df <- as.data.frame(res); df$sig <- "NS"
  df$sig[!is.na(df$padj) & df$padj < padj_cut & df$log2FoldChange >  lfc_cut] <- "Up"
  df$sig[!is.na(df$padj) & df$padj < padj_cut & df$log2FoldChange < -lfc_cut] <- "Down"
  p <- ggplot(df, aes(log2FoldChange, -log10(padj), color=sig)) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c(NS="grey70", Up="firebrick", Down="steelblue")) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed", color="grey40") +
    geom_hline(yintercept=-log10(padj_cut), linetype="dashed", color="grey40") +
    theme_bw() + labs(title=title, x="log2 fold change", y="-log10(padj)", color=NULL)
  ggsave(file, p, width=8, height=6, dpi=300)
}

# ---- Summary collector (one row per setting x contrast) ----
summary_rows <- list()

# ============================================================================
# Loop over threshold sets
# ============================================================================
for (ts in threshold_sets) {
  set_name <- ts$name
  set_dir  <- file.path(results_root, set_name)
  dir.create(set_dir, showWarnings=FALSE)
  message("\n########## THRESHOLD SET: ", set_name,
          "  (gene>=", ts$min_gene_count, " in ", ts$min_samples,
          " samples; padj<", ts$padj, "; lfc>", ts$lfc, ") ##########")

  # --- Expression filter for THIS setting ---
  group_vec <- metadata[[group_col]]
  passes <- counts_all >= ts$min_gene_count
  max_pass <- apply(passes, 1, function(r) max(tapply(r, group_vec, sum), na.rm=TRUE))
  counts_filt <- counts_all[max_pass >= ts$min_samples, , drop=FALSE]
  message("  Filtered miRNAs kept: ", nrow(counts_filt))

  # --- Save the normalized counts for this setting ---
  dge <- DGEList(counts_filt); dge <- calcNormFactors(dge, "TMM")
  write.table(cpm(dge, normalized.lib.sizes=TRUE, log=FALSE),
              file.path(set_dir, "TMM_normalized_counts.txt"),
              quote=FALSE, sep="\t", col.names=NA)

  # --- DESeq2 ---
  dds <- DESeqDataSetFromMatrix(counts_filt, metadata, design_formula)
  dds <- DESeq(dds)
  write.table(counts(dds, normalized=TRUE),
              file.path(set_dir, "DESeq2_normalized_counts.txt"),
              quote=FALSE, sep="\t", col.names=NA)

  condition_contrasts <- grep(paste0("^", group_col, "_"), resultsNames(dds), value=TRUE)

  # --- Per-contrast: results table, volcano, counts ---
  for (contrast in condition_contrasts) {
    res <- results(dds, name=contrast); res <- res[order(res$padj), ]
    tag <- gsub("[^A-Za-z0-9]+", "_", contrast)

    write.table(as.data.frame(res),
                file.path(set_dir, paste0("DESeq2_", tag, ".txt")),
                sep="\t", quote=FALSE, col.names=NA)
    make_volcano(res, paste0(set_name, ": ", contrast),
                 file.path(set_dir, paste0("Volcano_", tag, ".png")),
                 ts$padj, ts$lfc)

    n_total <- sum(res$padj < ts$padj, na.rm=TRUE)
    n_real  <- sum(res$padj < ts$padj & res$baseMean >= min_baseMean_real, na.rm=TRUE)
    n_real_lfc <- sum(res$padj < ts$padj & res$baseMean >= min_baseMean_real &
                      abs(res$log2FoldChange) > ts$lfc, na.rm=TRUE)
    message("    ", contrast, ": ", n_total, " total, ", n_real,
            " real(bM>=", min_baseMean_real, "), ", n_real_lfc, " real+lfc")

    summary_rows[[length(summary_rows)+1]] <- data.frame(
      setting=set_name, contrast=contrast,
      miRNAs_tested=nrow(counts_filt),
      sig_total=n_total, sig_real=n_real, sig_real_lfc=n_real_lfc,
      padj_cut=ts$padj, lfc_cut=ts$lfc,
      min_gene_count=ts$min_gene_count, min_samples=ts$min_samples)
  }
}

# ============================================================================
# Write the master summary table
# ============================================================================
summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, file.path(results_root, "SENSITIVITY_SUMMARY.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("\n================ SENSITIVITY SUMMARY ================")
print(summary_df[, c("setting","contrast","miRNAs_tested","sig_total","sig_real","sig_real_lfc")],
      row.names=FALSE)
message("\nColumns: sig_total = all padj<cutoff; sig_real = baseMean>=",
        min_baseMean_real, " (artifacts removed); sig_real_lfc = also |log2FC|>cutoff")
message("Full table: ", file.path(results_root, "SENSITIVITY_SUMMARY.txt"))
message("=== DONE ===")
