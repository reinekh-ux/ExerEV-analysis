set.seed(910)
suppressMessages({
  library(DESeq2); library(ggplot2); library(limma)
  library(reshape2); library(gridExtra); library(ggrepel)
})

results_dir <- "/lustre09/project/6019280/reinekh/Muscle_data_From_Minh/results"
batch1_dir  <- "/lustre09/project/6019280/reinekh/Muscle_data_From_Minh/data/Batch1"
fig_dir     <- file.path(results_dir, "manuscript_figures")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

contrasts  <- c("Exer","Exer3","Exer7","Exer18")
col_up <- "#E41A1C"
col_dn <- "#377EB8"
col_ns <- "grey70"
min_bm <- 10

# clean short name — handles pipe in dual annotations
short <- function(ids) {
  s <- sub(":.*","", sub("^mmu-","", ids))
  s <- sub("\\|.*","", s)
  s
}

# ── LOAD ALL RESULTS ────────────────────────────────────────────────────────
message("Loading results...")
res_list <- list()
all_sig  <- c()

for (contrast in contrasts) {
  fn  <- file.path(results_dir,
                   paste0("Muscle_miRNA_", contrast, "vsSed_allcontrasts.txt"))
  res <- read.table(fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  res$short_name <- short(rownames(res))
  res_list[[contrast]] <- res

  sig_names <- res$short_name[!is.na(res$padj) & res$padj < 0.05 &
                                res$baseMean >= min_bm]
  all_sig <- union(all_sig, sig_names)
}

message("Significant miRNAs (baseMean>=", min_bm, ", padj<0.05): ",
        length(all_sig))
message(paste(sort(all_sig), collapse=", "))

# ── FIGURE 1: VOLCANO PLOTS ─────────────────────────────────────────────────
message("\nMaking volcano plots...")
volcano_list <- list()

for (contrast in contrasts) {
  res <- res_list[[contrast]]

  # color — BOTH conditions must pass: padj < 0.05 AND baseMean >= min_bm
  res$col <- col_ns
  res$col[!is.na(res$padj) & res$padj < 0.05 &
            res$baseMean >= min_bm &
            res$log2FoldChange > 0] <- col_up
  res$col[!is.na(res$padj) & res$padj < 0.05 &
            res$baseMean >= min_bm &
            res$log2FoldChange < 0] <- col_dn

  # force grey for anything below baseMean threshold regardless
  res$col[res$baseMean < min_bm] <- col_ns

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
                     aes(label=short_name),
                     size=2.5, fontface="italic",
                     box.padding=0.35, point.padding=0.2,
                     max.overlaps=Inf,
                     min.segment.length=0,
                     segment.color="grey50",
                     label.size=0.2) +
    theme_bw(base_size=12) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=12)) +
    labs(title=paste0(contrast, " vs Sed"),
         x="log2 fold change", y="-log10(padj)") +
    xlim(c(-xlim_val, xlim_val))

  volcano_list[[contrast]] <- p
}

png(file.path(fig_dir, "Fig1_Volcano_v4.png"),
    width=14, height=11, units="in", res=300)
grid.arrange(grobs=volcano_list, ncol=2,
             top="Skeletal muscle sEV miRNA — Exercise vs Sed")
dev.off()
message("Volcano plots saved")

# ── FIGURE 2: TRAJECTORY PLOTS ──────────────────────────────────────────────
message("\nMaking trajectory plots for ", length(all_sig), " miRNAs...")

traj_list <- list()
for (contrast in contrasts) {
  res <- res_list[[contrast]]
  for (target in all_sig) {
    hit <- which(res$short_name == target)
    if (length(hit)==0) next
    r <- res[hit[1],]
    traj_list[[length(traj_list)+1]] <- data.frame(
      miRNA    = target,
      Contrast = contrast,
      log2FC   = r$log2FoldChange,
      SE       = r$lfcSE,
      padj     = ifelse(is.na(r$padj), 1, r$padj),
      baseMean = r$baseMean
    )
  }
}
traj_df <- do.call(rbind, traj_list)

sed_rows <- data.frame(miRNA=all_sig, Contrast="Sed",
                       log2FC=0, SE=0, padj=1, baseMean=NA)
traj_df  <- rbind(sed_rows, traj_df)
traj_df$Contrast <- factor(traj_df$Contrast,
                            levels=c("Sed","Exer","Exer3","Exer7","Exer18"))
traj_df$sig <- traj_df$padj < 0.05

dir_df <- traj_df[traj_df$Contrast != "Sed" & traj_df$sig, ]
dir_summary <- tapply(dir_df$log2FC, dir_df$miRNA,
                      function(x) ifelse(mean(x) > 0, "up", "down"))

traj_plots <- list()
for (target in sort(all_sig)) {
  df_t <- traj_df[traj_df$miRNA == target, ]
  dir  <- if (target %in% names(dir_summary)) dir_summary[[target]] else "up"
  clr  <- ifelse(dir == "up", col_up, col_dn)

  df_t$star        <- ifelse(df_t$sig & df_t$Contrast != "Sed", "*", "")
  star_offset      <- ifelse(dir == "up", 0.3, -0.3)

  p <- ggplot(df_t, aes(x=Contrast, y=log2FC, group=1)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
    geom_ribbon(aes(ymin=log2FC-SE, ymax=log2FC+SE),
                fill=clr, alpha=0.15) +
    geom_line(colour=clr, linewidth=1.2) +
    geom_point(colour=clr, fill=clr, size=3, shape=21, stroke=1.2) +
    geom_text(aes(label=star, y=log2FC + star_offset),
              size=5, colour="black") +
    theme_bw(base_size=10) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="italic", size=10),
          legend.position="none",
          axis.text.x=element_text(size=8)) +
    labs(title=target, x=NULL, y="log2FC vs Sed")

  traj_plots[[target]] <- p
}

n_mirnas <- length(traj_plots)
ncols    <- 4
nrows    <- ceiling(n_mirnas / ncols)
plot_h   <- max(10, nrows * 3.5)

png(file.path(fig_dir, "Fig2_Trajectory_v4.png"),
    width=16, height=plot_h, units="in", res=300)
grid.arrange(grobs=traj_plots, ncol=ncols,
             top="Temporal log2FC trajectory — skeletal muscle sEV miRNAs (baseMean>=10, padj<0.05)")
dev.off()
message("Trajectory plots saved — ", n_mirnas, " miRNAs, ",
        ncols, " cols x ", nrows, " rows")

# ── FIGURE 3: HEATMAP ───────────────────────────────────────────────────────
message("\nMaking heatmap...")

heat_data <- matrix(NA, nrow=length(all_sig), ncol=length(contrasts),
                    dimnames=list(sort(all_sig), contrasts))
padj_data <- heat_data

for (contrast in contrasts) {
  res <- res_list[[contrast]]
  for (target in all_sig) {
    hit <- which(res$short_name == target)
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

up_mirnas   <- sort(all_sig[all_sig %in%
                names(dir_summary[dir_summary=="up"])])
down_mirnas <- sort(all_sig[all_sig %in%
                names(dir_summary[dir_summary=="down"])])
hm$miRNA <- factor(hm$miRNA,
                   levels=c(rev(up_mirnas), rev(down_mirnas)))

ph <- ggplot(hm, aes(x=Contrast, y=miRNA, fill=log2FC)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(aes(label=sig_label), size=5, vjust=0.75, colour="black") +
  scale_fill_gradient2(low=col_dn, mid="white", high=col_up,
                       midpoint=0, limits=c(-lim, lim),
                       name="log2FC vs Sed") +
  theme_bw(base_size=11) +
  theme(axis.text.y=element_text(face="italic", size=9),
        panel.grid=element_blank()) +
  labs(title="Skeletal muscle sEV miRNA — log2FC vs Sed",
       subtitle="* padj < 0.05 | baseMean >= 10",
       x=NULL, y=NULL)

heatmap_h <- max(5, length(all_sig) * 0.4 + 2)
ggsave(file.path(fig_dir, "Fig3_Heatmap_v4.png"),
       ph, width=8, height=heatmap_h, dpi=300)
message("Heatmap saved")

message("\nAll v4 figures saved to: ", fig_dir)
message("=== DONE ===")
