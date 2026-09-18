###############################################################################
# Exploratory variance analysis (PCA) — BEFORE choosing a DE model
#
# Project : Reine - Bone Marrow EV small-RNA sequencing (exercise time course)
# Purpose : Look at the data model-free. Which variables actually drive
#           variation? Colour the PCA by each design variable AND quantify how
#           much variance each explains on the top PCs. NO DE model is fitted.
#
# Use the output to DECIDE the DESeq2 design, rather than assuming it.
#
# How to run:
#   module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
#   Rscript explore_PCA_BoneMarrow.R
###############################################################################

set.seed(910)

suppressMessages({
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
})

# ============================================================================
# SETTINGS
# ============================================================================
working_dir   <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
setwd(working_dir)

rawcount_file <- "count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt"
metadata_file <- "metadata_BoneMarrow.txt"
sample_id_col <- "Index_Number"
group_col     <- "Condition"
reference_level <- "Sed"

# Variables to inspect (must be columns in the metadata).
inspect_vars <- c("Condition", "Sex", "Isolation_Day", "Seq_Batch", "Treadmill_Batch")

# Filtering (same as the DE script, so the exploration reflects the real data).
min_gene_count       <- 10
min_samples_in_group <- 3

# How many top PCs to test in the quantitative variance check.
n_pcs <- 5

out_prefix <- "BoneMarrow_miRNA_EXPLORE"
# ============================================================================


# ----------------------------------------------------------------------------
# 1. Load counts + metadata
# ----------------------------------------------------------------------------
raw_counts <- read.table(rawcount_file, header = TRUE, row.names = 1,
                         sep = "\t", check.names = FALSE)
raw_counts <- round(raw_counts)
message("Raw: ", nrow(raw_counts), " miRNAs x ", ncol(raw_counts), " samples")

metadata <- read.table(metadata_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
rownames(metadata) <- metadata[[sample_id_col]]
metadata <- metadata[match(colnames(raw_counts), rownames(metadata)), ]
stopifnot(identical(colnames(raw_counts), rownames(metadata)))

# Make inspected variables factors; put reference first for the group.
for (v in inspect_vars) if (v %in% colnames(metadata)) metadata[[v]] <- factor(metadata[[v]])
metadata[[group_col]] <- relevel(factor(metadata[[group_col]]), ref = reference_level)


# ----------------------------------------------------------------------------
# 2. Group-aware gene filter (same as DE script)
# ----------------------------------------------------------------------------
group_vec <- metadata[[group_col]]
passes <- raw_counts >= min_gene_count
max_pass <- apply(passes, 1, function(r) max(tapply(r, group_vec, sum), na.rm = TRUE))
counts_filt <- raw_counts[max_pass >= min_samples_in_group, , drop = FALSE]
message("After filter: ", nrow(counts_filt), " miRNAs")


# ----------------------------------------------------------------------------
# 3. VST normalization (model-free: design = ~1)
# ----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                              colData = metadata,
                              design = ~ 1)          # <- no model, purely exploratory
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsmat <- assay(vsd)


# ----------------------------------------------------------------------------
# 4. PCA coloured by each variable
# ----------------------------------------------------------------------------
# Compute PCA once on the top-variable genes (DESeq2's plotPCA default = 500).
ntop <- min(500, nrow(vsmat))
rv <- apply(vsmat, 1, var)
select <- order(rv, decreasing = TRUE)[seq_len(ntop)]
pca <- prcomp(t(vsmat[select, ]))
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- as.data.frame(pca$x)
pca_df <- cbind(pca_df, metadata)

save_pca_colored <- function(var, file) {
  if (!(var %in% colnames(pca_df))) return(invisible())
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[var]])) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(paste("PCA coloured by", var)) +
    theme_bw() +
    theme(axis.title = element_text(size = 14, face = "bold"),
          axis.text  = element_text(size = 12),
          plot.title = element_text(size = 15, face = "bold"))
  ggsave(file, p, width = 8, height = 6, dpi = 300)
  message("Saved: ", file)
}

for (v in inspect_vars) {
  save_pca_colored(v, paste0(out_prefix, "_PCA_", v, ".png"))
}


# ----------------------------------------------------------------------------
# 5. QUANTITATIVE check — how much variance does each variable explain per PC?
#    For each top PC, regress the PC scores on each variable (one at a time)
#    and record the R^2. Higher R^2 = that variable explains more of that PC.
# ----------------------------------------------------------------------------
message("\n==== Variance explained (R^2) by each variable, per PC ====")
pc_use <- min(n_pcs, ncol(pca$x))
r2_table <- matrix(NA, nrow = length(inspect_vars), ncol = pc_use,
                   dimnames = list(inspect_vars, paste0("PC", seq_len(pc_use))))

for (v in inspect_vars) {
  if (!(v %in% colnames(metadata))) next
  for (k in seq_len(pc_use)) {
    fit <- lm(pca$x[, k] ~ metadata[[v]])
    r2_table[v, k] <- summary(fit)$r.squared
  }
}

# Print a clean table (R^2 as %).
r2_pct <- round(100 * r2_table, 1)
print(r2_pct)

write.table(r2_pct, file = paste0(out_prefix, "_variance_explained_byPC.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
message("Saved variance table: ", out_prefix, "_variance_explained_byPC.txt")

# Also report the % variance each PC itself carries.
message("\n==== % of total variance per PC ====")
names(percentVar) <- paste0("PC", seq_along(percentVar))
print(percentVar[seq_len(pc_use)])


# ----------------------------------------------------------------------------
# 6. Sample-distance heatmap (unsupervised clustering)
# ----------------------------------------------------------------------------
sampleDists <- dist(t(vsmat))
distMat <- as.matrix(sampleDists)
# Annotate the heatmap rows with the variables so clustering is interpretable.
annot <- metadata[, intersect(inspect_vars, colnames(metadata)), drop = FALSE]
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
png(paste0(out_prefix, "_sample_distance_heatmap.png"), width = 11, height = 9,
    units = "in", res = 300)
pheatmap(distMat,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = annot,
         col = colors)
dev.off()
message("Saved heatmap: ", out_prefix, "_sample_distance_heatmap.png")

message("\n=== EXPLORATION DONE ===")
message("Interpret: whichever variable has the highest R^2 on PC1/PC2 (and lines")
message("up with the coloured PCA clustering) is a major driver of variation and")
message("should be accounted for in the DE model. Use this to confirm or revise")
message("the design before running DEG_miRNA_BoneMarrow.R.")
