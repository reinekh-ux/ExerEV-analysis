set.seed(910)
suppressMessages({
  library(DESeq2); library(ggplot2); library(limma); library(reshape2)
})

# ── PATHS ──────────────────────────────────────────────────────────────────
batch1_dir    <- "/lustre09/project/6019280/reinekh/Muscle_data_From_Minh/data/Batch1"
rawcount_file <- file.path(batch1_dir, "exceRpt_rawcounts_perfect_match_mm10.csv")
metadata_file <- file.path(batch1_dir, "metadata.csv")
results_dir   <- "/lustre09/project/6019280/reinekh/Muscle_data_From_Minh/results"
fig_dir       <- file.path(results_dir, "manuscript_figures")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

# significant miRNAs to highlight
sig_mirnas <- c("miR-133b-3p","miR-199a-5p","miR-29a-3p","miR-92a-3p")
sig_timepoints <- c("Exer","Exer3","Exer7","Exer18")
contrasts  <- c("Exer","Exer3","Exer7","Exer18")

col_up   <- "#E41A1C"   # red  = upregulated
col_dn   <- "#377EB8"   # blue = downregulated
col_ns   <- "grey70"    # grey = not significant

# short label helper
short <- function(ids) sub(":.*","", sub("^mmu-","", ids))

# ── LOAD DATA for line/heatmap plots ───────────────────────────────────────
raw  <- round(read.csv(rawcount_file, row.names=1, check.names=FALSE))
meta <- read.csv(metadata_file, stringsAsFactors=FALSE)
rownames(meta) <- meta$Number_Index
raw  <- raw[, colnames(raw) %in% rownames(meta), drop=FALSE]
meta <- meta[match(colnames(raw), rownames(meta)), ]

meta$Condition <- meta$condition
meta$Condition[meta$Condition == "Ex + 3"]  <- "Exer3"
meta$Condition[meta$Condition == "Ex + 7"]  <- "Exer7"
meta$Condition[meta$Condition == "Ex + 18"] <- "Exer18"
meta$Condition[meta$Condition == "Pre Ex"]  <- "PreEx"
meta <- meta[meta$Condition != "PreEx", ]
raw  <- raw[, rownames(meta), drop=FALSE]
meta$Condition     <- factor(meta$Condition,
                             levels=c("Sed","Exer","Exer3","Exer7","Exer18"))
meta$Sex           <- factor(meta$sex)
meta$Isolation_Day <- factor(meta$animal_id)

passes <- raw >= 10
keep   <- apply(passes, 1, function(r)
             max(tapply(r, meta$Condition, sum), na.rm=TRUE)) >= 3
raw    <- raw[keep, ]

meta_d           <- meta
meta_d$Condition <- relevel(meta_d$Condition, ref="Sed")
dds <- DESeqDataSetFromMatrix(raw, meta_d, ~ Isolation_Day + Sex + Condition)
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized=TRUE)
logn <- log2(norm + 1)
design_keep <- model.matrix(~ Condition + Sex, data=meta_d)
logn_adj    <- removeBatchEffect(logn, batch=meta_d$Isolation_Day,
                                 design=design_keep)

# ── FIGURE 1: VOLCANO PLOTS (2x2 grid) ────────────────────────────────────
message("Making volcano plots...")

volcano_list <- list()

for (contrast in contrasts) {
  fn <- file.path(results_dir,
                  paste0("Muscle_miRNA_", contrast, "vsSed_allcontrasts.txt"))
  res <- read.table(fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  res$label <- short(rownames(res))

  # color coding
  res$col <- col_ns
  res$col[!is.na(res$padj) & res$padj < 0.05 &
            res$log2FoldChange > 0] <- col_up
  res$col[!is.na(res$padj) & res$padj < 0.05 &
            res$log2FoldChange < 0] <- col_dn
  res$sig <- res$col != col_ns

  # label only significant hits
  res$show_label <- res$label %in% sig_mirnas & res$sig

  # y axis
  res$neglog10p <- -log10(res$padj)
  res$neglog10p[is.na(res$neglog10p)] <- 0

  p <- ggplot(res, aes(x=log2FoldChange, y=neglog10p)) +
    geom_point(aes(color=col), size=1.8, alpha=0.7) +
    scale_color_identity() +
    geom_hline(yintercept=-log10(0.05), linetype="dashed",
               colour="black", linewidth=0.5) +
    geom_vline(xintercept=c(-1,1), linetype="dotted",
               colour="grey40", linewidth=0.4) +
    geom_text(data=subset(res, show_label),
              aes(label=label), size=3.2, vjust=-0.6, fontface="italic") +
    theme_bw(base_size=12) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=12)) +
    labs(title=paste0(contrast, " vs Sed"),
         x="log2 fold change", y="-log10(padj)") +
    xlim(c(-max(abs(res$log2FoldChange), na.rm=TRUE)-0.5,
            max(abs(res$log2FoldChange), na.rm=TRUE)+0.5))

  volcano_list[[contrast]] <- p
}

# arrange 2x2
library(gridExtra)
png(file.path(fig_dir, "Fig1_Volcano_allTimepoints.png"),
    width=12, height=10, units="in", res=300)
grid.arrange(grobs=volcano_list, ncol=2,
             top="Skeletal muscle sEV miRNA — Exercise vs Sed")
dev.off()
message("Volcano plots saved")

# ── FIGURE 2: LOG2FC TRAJECTORY LINE PLOTS ────────────────────────────────
message("Making trajectory plots...")

# collect log2FC per miRNA per timepoint from results files
traj_list <- list()
for (contrast in contrasts) {
  fn <- file.path(results_dir,
                  paste0("Muscle_miRNA_", contrast, "vsSed_allcontrasts.txt"))
  res <- read.table(fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  for (target in sig_mirnas) {
    hit <- grep(paste0("mmu-", target, ":"), rownames(res))
    if (length(hit)==0) next
    r <- res[hit[1],]
    traj_list[[length(traj_list)+1]] <- data.frame(
      miRNA    = target,
      Contrast = contrast,
      log2FC   = r$log2FoldChange,
      SE       = r$lfcSE,
      padj     = r$padj
    )
  }
}
traj_df <- do.call(rbind, traj_list)
# add Sed = 0 baseline
sed_rows <- data.frame(
  miRNA    = sig_mirnas,
  Contrast = "Sed",
  log2FC   = 0,
  SE       = 0,
  padj     = 1
)
traj_df  <- rbind(sed_rows, traj_df)
traj_df$Contrast <- factor(traj_df$Contrast,
                            levels=c("Sed","Exer","Exer3","Exer7","Exer18"))
traj_df$Direction <- ifelse(traj_df$miRNA == "miR-92a-3p", "down", "up")
traj_df$sig_point <- !is.na(traj_df$padj) & traj_df$padj < 0.05

# one panel per miRNA
traj_plots <- list()
for (target in sig_mirnas) {
  df_t <- traj_df[traj_df$miRNA == target, ]
  dir   <- unique(df_t$Direction)
  clr   <- ifelse(dir == "up", col_up, col_dn)

  p <- ggplot(df_t, aes(x=Contrast, y=log2FC, group=1)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
    geom_ribbon(aes(ymin=log2FC-SE, ymax=log2FC+SE),
                fill=clr, alpha=0.15) +
    geom_line(colour=clr, linewidth=1) +
    geom_point(aes(shape=sig_point), colour=clr, size=3, fill=clr) +
    scale_shape_manual(values=c("FALSE"=21, "TRUE"=16),
                       labels=c("ns","padj<0.05"),
                       name="") +
    theme_bw(base_size=12) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="italic", size=11),
          legend.position="bottom") +
    labs(title=target,
         x=NULL, y="log2 fold change vs Sed")
  traj_plots[[target]] <- p
}

png(file.path(fig_dir, "Fig2_Trajectory_log2FC.png"),
    width=12, height=10, units="in", res=300)
grid.arrange(grobs=traj_plots, ncol=2,
             top="Temporal log2FC trajectory — skeletal muscle sEV miRNAs")
dev.off()
message("Trajectory plots saved")

# ── FIGURE 3: HEATMAP of log2FC ───────────────────────────────────────────
message("Making heatmap...")

# build log2FC matrix
heat_data <- matrix(NA, nrow=length(sig_mirnas), ncol=length(contrasts),
                    dimnames=list(sig_mirnas, contrasts))
padj_data <- heat_data

for (contrast in contrasts) {
  fn <- file.path(results_dir,
                  paste0("Muscle_miRNA_", contrast, "vsSed_allcontrasts.txt"))
  res <- read.table(fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  for (target in sig_mirnas) {
    hit <- grep(paste0("mmu-", target, ":"), rownames(res))
    if (length(hit)==0) next
    heat_data[target, contrast] <- res[hit[1], "log2FoldChange"]
    padj_data[target, contrast] <- res[hit[1], "padj"]
  }
}

# melt for ggplot
hm <- melt(heat_data, varnames=c("miRNA","Contrast"), value.name="log2FC")
pm <- melt(padj_data, varnames=c("miRNA","Contrast"), value.name="padj")
hm$padj <- pm$padj
hm$sig_label <- ifelse(!is.na(hm$padj) & hm$padj < 0.05, "*", "")
hm$Contrast <- factor(hm$Contrast,
                       levels=c("Exer","Exer3","Exer7","Exer18"))

lim <- max(abs(heat_data), na.rm=TRUE)

ph <- ggplot(hm, aes(x=Contrast, y=miRNA, fill=log2FC)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(aes(label=sig_label), size=6, vjust=0.75, colour="black") +
  scale_fill_gradient2(low=col_dn, mid="white", high=col_up,
                       midpoint=0, limits=c(-lim, lim),
                       name="log2FC vs Sed") +
  theme_bw(base_size=12) +
  theme(axis.text.y=element_text(face="italic"),
        panel.grid=element_blank()) +
  labs(title="Skeletal muscle sEV miRNA — log2FC vs Sed",
       subtitle="* padj < 0.05",
       x=NULL, y=NULL)

ggsave(file.path(fig_dir, "Fig3_Heatmap_log2FC.png"),
       ph, width=7, height=4, dpi=300)
message("Heatmap saved")

message("\nAll figures saved to: ", fig_dir)
message("=== DONE ===")
