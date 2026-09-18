###############################################################################
# Pre-DE QC checks — run BEFORE DESeq2 to confirm samples are clean
#
# Project : Reine - Bone Marrow EV small-RNA sequencing
# Checks  : 1. Library size / depth per sample (flag shallow samples)
#           2. Sex separation (flag possible label swaps)
#           3. Outlier detection in PCA space (flag rogue samples)
#
# How to run:
#   module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
#   Rscript preDE_QC_BoneMarrow.R
###############################################################################

set.seed(910)
suppressMessages({
  library(DESeq2); library(edgeR); library(ggplot2); library(RColorBrewer)
})

# ---- SETTINGS ----
working_dir <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
setwd(working_dir)
rawcount_file <- "count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt"
metadata_file <- "metadata_BoneMarrow.txt"
sample_id_col <- "Index_Number"
group_col     <- "Condition"
reference_level <- "Sed"
min_gene_count <- 10; min_samples_in_group <- 3
out_prefix <- "BoneMarrow_miRNA_preDEQC"
# ------------------

# Load
raw_counts <- read.table(rawcount_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
raw_counts <- round(raw_counts)
metadata <- read.table(metadata_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(metadata) <- metadata[[sample_id_col]]
metadata <- metadata[match(colnames(raw_counts), rownames(metadata)), ]
stopifnot(identical(colnames(raw_counts), rownames(metadata)))
metadata$Sex <- factor(metadata$Sex)
metadata[[group_col]] <- relevel(factor(metadata[[group_col]]), ref=reference_level)


# ============================================================================
# CHECK 1 — Library size / depth per sample
# ============================================================================
message("\n===== CHECK 1: Library size per sample =====")
libsize <- colSums(raw_counts)
libsize_sorted <- sort(libsize)
depth_df <- data.frame(Sample = names(libsize), LibSize = libsize,
                       Condition = metadata[names(libsize), group_col],
                       Sex = metadata[names(libsize), "Sex"])

# Summary stats
message("Median library size: ", format(median(libsize), big.mark=","))
message("Min: ", format(min(libsize), big.mark=","), " (", names(which.min(libsize)), ")")
message("Max: ", format(max(libsize), big.mark=","), " (", names(which.max(libsize)), ")")

# Flag samples below 25% of the median as potentially shallow
thresh <- 0.25 * median(libsize)
shallow <- names(libsize)[libsize < thresh]
if (length(shallow)) {
  message("POTENTIALLY SHALLOW (< 25% of median): ", paste(shallow, collapse=", "))
} else {
  message("No samples below 25% of median depth — all adequate.")
}

# Show the 5 lowest-depth samples
message("\n5 lowest-depth samples:")
print(head(libsize_sorted, 5))

# Barplot of library sizes, coloured by condition
depth_df$Sample <- factor(depth_df$Sample, levels = names(libsize_sorted))
p1 <- ggplot(depth_df, aes(x=Sample, y=LibSize, fill=Condition)) +
  geom_col() +
  geom_hline(yintercept = median(libsize), linetype="dashed") +
  geom_hline(yintercept = thresh, linetype="dotted", color="red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=7)) +
  labs(title="Library size per sample (dashed=median, dotted red=25% threshold)",
       y="Total miRNA counts")
ggsave(paste0(out_prefix, "_1_librarysize.png"), p1, width=12, height=6, dpi=300)
message("Saved: ", out_prefix, "_1_librarysize.png")


# ============================================================================
# Prepare normalized data for CHECKS 2 & 3 (filter + VST)
# ============================================================================
group_vec <- metadata[[group_col]]
passes <- raw_counts >= min_gene_count
max_pass <- apply(passes, 1, function(r) max(tapply(r, group_vec, sum), na.rm=TRUE))
counts_filt <- raw_counts[max_pass >= min_samples_in_group, , drop=FALSE]

dds <- DESeqDataSetFromMatrix(counts_filt, metadata, design = ~ 1)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vsmat <- assay(vsd)

# PCA
ntop <- min(500, nrow(vsmat))
rv <- apply(vsmat, 1, var)
select <- order(rv, decreasing=TRUE)[seq_len(ntop)]
pca <- prcomp(t(vsmat[select, ]))
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(pca$x, metadata)


# ============================================================================
# CHECK 2 — Sex separation (label-swap detector)
# ============================================================================
message("\n===== CHECK 2: Sex separation =====")
# Approach: how well does Sex explain each of the top PCs? And do the two sexes
# separate? A clean separation on some PC indicates labels are internally
# consistent. We report R^2 of Sex on each PC and, for the PC where Sex is
# strongest, list any sample sitting on the 'wrong' side of the sex boundary.

pc_use <- min(5, ncol(pca$x))
sex_r2 <- sapply(seq_len(pc_use), function(k)
  summary(lm(pca$x[,k] ~ metadata$Sex))$r.squared)
names(sex_r2) <- paste0("PC", seq_len(pc_use))
message("Sex R^2 per PC (%):")
print(round(100*sex_r2, 1))

best_pc <- which.max(sex_r2)
message("Sex separates best on PC", best_pc, " (R^2 = ", round(100*sex_r2[best_pc],1), "%)")

# On the best PC, compute each sex's mean; flag samples closer to the other sex's mean.
scores <- pca$x[, best_pc]
male_mean   <- mean(scores[metadata$Sex == "Male"])
female_mean <- mean(scores[metadata$Sex == "Female"])
predicted_sex <- ifelse(abs(scores - male_mean) < abs(scores - female_mean), "Male", "Female")
mismatch <- names(scores)[predicted_sex != as.character(metadata$Sex)]
if (length(mismatch)) {
  message("SAMPLES whose PC", best_pc,
          " position sits closer to the opposite sex (INSPECT for swaps):")
  print(mismatch)
  message("NOTE: some mismatches are normal if Sex is a weak PC driver; ",
          "investigate only if Sex R^2 is high and a sample is a clear outlier.")
} else {
  message("No samples sit on the wrong side of the sex boundary on PC", best_pc, ".")
}

# Plot PC (best sex PC) vs PC1, coloured by Sex
sex_plot_pc <- if (best_pc == 1) 2 else best_pc
p2 <- ggplot(pca_df, aes(x=PC1, y=.data[[paste0("PC", sex_plot_pc)]], color=Sex, label=rownames(pca_df))) +
  geom_point(size=3) +
  geom_text(size=2, vjust=-1, show.legend=FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC", sex_plot_pc, ": ", percentVar[sex_plot_pc], "%")) +
  theme_bw() + labs(title="Sex check — do the sexes separate? (labels = samples)")
ggsave(paste0(out_prefix, "_2_sexcheck.png"), p2, width=9, height=7, dpi=300)
message("Saved: ", out_prefix, "_2_sexcheck.png")


# ============================================================================
# CHECK 3 — Outlier detection in PCA space
# ============================================================================
message("\n===== CHECK 3: Outlier detection =====")
# Method: for each sample, distance from the overall centroid in the space of
# the top PCs. Flag samples > 3 SD beyond the mean distance (classic outlier).
topk <- min(5, ncol(pca$x))
coords <- pca$x[, seq_len(topk), drop=FALSE]
centroid <- colMeans(coords)
dists <- sqrt(rowSums((sweep(coords, 2, centroid))^2))
z <- (dists - mean(dists)) / sd(dists)

outlier_df <- data.frame(Sample=names(dists), Distance=round(dists,2),
                         Zscore=round(z,2))
outlier_df <- outlier_df[order(-outlier_df$Zscore), ]
message("Top 5 most distant samples from centroid (top ", topk, " PCs):")
print(head(outlier_df, 5))

flagged <- outlier_df$Sample[outlier_df$Zscore > 3]
if (length(flagged)) {
  message("OUTLIERS (> 3 SD from centroid) — INSPECT before DE: ",
          paste(flagged, collapse=", "))
} else {
  message("No samples exceed 3 SD from centroid — no clear outliers.")
}

# Labeled PCA to eyeball outliers
p3 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=.data[[group_col]], label=rownames(pca_df))) +
  geom_point(size=3) +
  geom_text(size=2, vjust=-1, show.legend=FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw() + labs(title="Outlier check — PCA with sample labels")
ggsave(paste0(out_prefix, "_3_outliercheck.png"), p3, width=9, height=7, dpi=300)
message("Saved: ", out_prefix, "_3_outliercheck.png")

# Save the numeric tables
write.table(depth_df[order(depth_df$LibSize), ],
            paste0(out_prefix, "_librarysize.txt"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(outlier_df, paste0(out_prefix, "_outlier_distances.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)


# ============================================================================
message("\n===== PRE-DE QC SUMMARY =====")
message("Samples: ", ncol(raw_counts))
message("Shallow (<25% median depth): ", if(length(shallow)) paste(shallow, collapse=", ") else "none")
message("Sex strongest on PC", best_pc, " (R^2 ", round(100*sex_r2[best_pc],1), "%); ",
        "sex mismatches: ", if(length(mismatch)) paste(mismatch, collapse=", ") else "none")
message("Outliers (>3 SD): ", if(length(flagged)) paste(flagged, collapse=", ") else "none")
message("\nReview the three PNGs and the flags above. If all clear, proceed to DESeq2.")
message("=== DONE ===")
