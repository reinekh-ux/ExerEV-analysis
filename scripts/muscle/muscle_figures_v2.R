set.seed(910)
suppressMessages({
  library(DESeq2); library(ggplot2); library(limma)
  library(reshape2); library(gridExtra); library(ggrepel)
})

results_dir <- "/lustre09/project/6019280/reinekh/Muscle_data_From_Minh/results"
batch1_dir  <- "/lustre09/project/6019280/reinekh/Muscle_data_From_Minh/data/Batch1"
fig_dir     <- file.path(results_dir, "manuscript_figures")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

sig_mirnas  <- c("miR-133b-3p","miR-199a-5p","miR-29a-3p","miR-92a-3p")
contrasts   <- c("Exer","Exer3","Exer7","Exer18")
col_up <- "#E41A1C"
col_dn <- "#377EB8"
col_ns <- "grey70"
short  <- function(ids) sub(":.*","", sub("^mmu-","", ids))

# ── FIGURE 1: VOLCANO PLOTS ────────────────────────────────────────────────
message("Making volcano plots...")
volcano_list <- list()

for (contrast in contrasts) {
  fn  <- file.path(results_dir,
                   paste0("Muscle_miRNA_", contrast, "vsSed_allcontrasts.txt"))
  res <- read.table(fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  res$label <- short(rownames(res))

  res$col <- col_ns
  res$col[!is.na(res$padj) & res$padj < 0.05 &
            res$log2FoldChange > 0 & res$baseMean >= 10] <- col_up
  res$col[!is.na(res$padj) & res$padj < 0.05 &
            res$log2FoldChange < 0 & res$baseMean >= 10] <- col_dn

  # label ALL significant hits with baseMean >= 10
  res$show_label <- res$col != col_ns
  res$neglog10p  <- ifelse(is.na(res$padj), 0, -log10(res$padj))

  xlim_val <- max(abs(res$log2FoldChange), na.rm=TRUE) + 0.5

  p <- ggplot(res, aes(x=log2FoldChange, y=neglog10p)) +
    geom_point(aes(color=col), size=1.8, alpha=0.7) +
    scale_color_identity() +
    geom_hline(yintercept=-log10(0.05), linetype="dashed",
               colour="black", linewidth=0.5) +
    geom_vline(xintercept=c(-1,1), linetype="dotted",
               colour="grey40", linewidth=0.4) +
    geom_label_repel(data=subset(res, show_label),
                     aes(label=label),
                     size=3, fontface="italic",
                     box.padding=0.4, point.padding=0.3,
                     max.overlaps=20,
                     segment.color="grey50") +
    theme_bw(base_size=12) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=12)) +
    labs(title=paste0(contrast, " vs Sed"),
         x="log2 fold change", y="-log10(padj)") +
    xlim(c(-xlim_val, xlim_val))

  volcano_list[[contrast]] <- p
}

png(file.path(fig_dir, "Fig1_Volcano_allTimepoints_v2.png"),
    width=14, height=11, units="in", res=300)
grid.arrange(grobs=volcano_list, ncol=2,
             top="Skeletal muscle sEV miRNA — Exercise vs Sed")
dev.off()
message("Volcano plots saved")

# ── FIGURE 2: TRAJECTORY PLOTS (fixed) ────────────────────────────────────
message("Making trajectory plots...")

traj_list <- list()
for (contrast in contrasts) {
  fn  <- file.path(results_dir,
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
      padj     = ifelse(is.na(r$padj), 1, r$padj)
    )
  }
}
traj_df <- do.call(rbind, traj_list)

# add Sed baseline
sed_rows <- data.frame(miRNA=sig_mirnas, Contrast="Sed",
                       log2FC=0, SE=0, padj=1)
traj_df  <- rbind(sed_rows, traj_df)
traj_df$Contrast <- factor(traj_df$Contrast,
                            levels=c("Sed","Exer","Exer3","Exer7","Exer18"))
traj_df$sig <- traj_df$padj < 0.05
traj_df$Direction <- ifelse(traj_df$miRNA == "miR-92a-3p", "down", "up")

traj_plots <- list()
for (target in sig_mirnas) {
  df_t <- traj_df[traj_df$miRNA == target, ]
  clr  <- ifelse(unique(df_t$Direction) == "up", col_up, col_dn)

  # star label for significant timepoints (not Sed)
  df_t$star <- ifelse(df_t$sig & df_t$Contrast != "Sed", "*", "")

  p <- ggplot(df_t, aes(x=Contrast, y=log2FC, group=1)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
    geom_ribbon(aes(ymin=log2FC-SE, ymax=log2FC+SE),
                fill=clr, alpha=0.15) +
    geom_line(colour=clr, linewidth=1.2) +
    geom_point(colour=clr, fill=clr, size=3.5, shape=21,
               stroke=1.2) +
    geom_text(aes(label=star, y=log2FC + sign(log2FC)*0.3),
              size=6, colour="black", vjust=0.5) +
    theme_bw(base_size=12) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="italic", size=12),
          legend.position="none") +
    labs(title=target, x=NULL, y="log2 fold change vs Sed")

  traj_plots[[target]] <- p
}

png(file.path(fig_dir, "Fig2_Trajectory_log2FC_v2.png"),
    width=12, height=10, units="in", res=300)
grid.arrange(grobs=traj_plots, ncol=2,
             top="Temporal log2FC trajectory — skeletal muscle sEV miRNAs")
dev.off()
message("Trajectory plots saved")

# ── FIGURE 3: HEATMAP ─────────────────────────────────────────────────────
message("Making heatmap...")

heat_data <- matrix(NA, nrow=length(sig_mirnas), ncol=length(contrasts),
                    dimnames=list(sig_mirnas, contrasts))
padj_data <- heat_data

for (contrast in contrasts) {
  fn  <- file.path(results_dir,
                   paste0("Muscle_miRNA_", contrast, "vsSed_allcontrasts.txt"))
  res <- read.table(fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  for (target in sig_mirnas) {
    hit <- grep(paste0("mmu-", target, ":"), rownames(res))
    if (length(hit)==0) next
    heat_data[target, contrast] <- res[hit[1], "log2FoldChange"]
    padj_data[target, contrast] <- res[hit[1], "padj"]
  }
}

hm <- melt(heat_data, varnames=c("miRNA","Contrast"), value.name="log2FC")
pm <- melt(padj_data, varnames=c("miRNA","Contrast"), value.name="padj")
hm$padj      <- pm$padj
hm$sig_label <- ifelse(!is.na(hm$padj) & hm$padj < 0.05, "*", "")
hm$Contrast  <- factor(hm$Contrast, levels=contrasts)
lim <- max(abs(heat_data), na.rm=TRUE)

ph <- ggplot(hm, aes(x=Contrast, y=miRNA, fill=log2FC)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(aes(label=sig_label), size=7, vjust=0.75, colour="black") +
  scale_fill_gradient2(low=col_dn, mid="white", high=col_up,
                       midpoint=0, limits=c(-lim, lim),
                       name="log2FC vs Sed") +
  theme_bw(base_size=12) +
  theme(axis.text.y=element_text(face="italic"),
        panel.grid=element_blank()) +
  labs(title="Skeletal muscle sEV miRNA — log2FC vs Sed",
       subtitle="* padj < 0.05", x=NULL, y=NULL)

ggsave(file.path(fig_dir, "Fig3_Heatmap_log2FC_v2.png"),
       ph, width=7, height=4, dpi=300)
message("Heatmap saved")

message("\nAll v2 figures in: ", fig_dir)
message("=== DONE ===")
