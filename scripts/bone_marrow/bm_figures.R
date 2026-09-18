set.seed(910)
suppressMessages({
  library(DESeq2); library(ggplot2); library(limma)
  library(reshape2); library(gridExtra); library(ggrepel)
})

working_dir <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
setwd(working_dir)
results_dir <- "DESeq2_results/BM_manuscript"
fig_dir     <- file.path(results_dir, "figures")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

contrasts <- c("Exer","Exer3","Exer7","Exer18")
col_up <- "#E41A1C"
col_dn <- "#377EB8"
col_ns <- "grey70"
min_bm <- 10

short <- function(ids) {
  s <- sub(":.*","", sub("^mmu-","", ids))
  sub("\\|.*","", s)
}

# ── LOAD ALL CONTRAST RESULTS ──────────────────────────────────────────────
res_list <- list()
all_sig  <- c()

for (contrast in contrasts) {
  fn  <- file.path(results_dir,
                   paste0("BM_miRNA_", contrast, "vsSed.txt"))
  res <- read.table(fn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  res$short_name <- short(rownames(res))
  res_list[[contrast]] <- res
  sig_names <- res$short_name[!is.na(res$padj) & res$padj < 0.05 &
                                res$baseMean >= min_bm]
  all_sig <- union(all_sig, sig_names)
}

message("Significant BM miRNAs (baseMean>=", min_bm,
        ", padj<0.05): ", length(all_sig))
if (length(all_sig) > 0) message(paste(sort(all_sig), collapse=", "))

# ── FIGURE 1: P-VALUE HISTOGRAM (the key bone marrow figure) ──────────────
pooled_file <- file.path(results_dir, "BM_miRNA_PooledExercise_vs_Sed.txt")
pooled <- read.table(pooled_file, header=TRUE, sep="\t",
                     stringsAsFactors=FALSE)
pooled <- pooled[!is.na(pooled$pvalue), ]
exp_line <- nrow(pooled) * 0.05 / 20   # expected per bin

ph <- ggplot(pooled, aes(x=pvalue)) +
  geom_histogram(breaks=seq(0,1,0.05), fill="grey60",
                 colour="white", linewidth=0.3) +
  geom_hline(yintercept=exp_line, linetype="dashed",
             colour="#E41A1C", linewidth=0.8) +
  annotate("text", x=0.85, y=exp_line*1.15,
           label=paste0("expected = ", round(exp_line,1)),
           colour="#E41A1C", size=3.5) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank()) +
  labs(title="Bone marrow sEV miRNA — p-value distribution (pooled Exercise vs Sed)",
       subtitle=paste0("Observed p<0.05: ",
                       sum(pooled$pvalue < 0.05),
                       " | Expected by chance: ",
                       round(nrow(pooled)*0.05, 1)),
       x="raw p-value", y="number of miRNAs")

ggsave(file.path(fig_dir, "Fig0_pvalue_histogram.png"),
       ph, width=8, height=5, dpi=300)
message("p-value histogram saved")

# ── FIGURE 2: VOLCANO PLOTS ────────────────────────────────────────────────
message("Making volcano plots...")
volcano_list <- list()

for (contrast in contrasts) {
  res <- res_list[[contrast]]
  res$col <- col_ns
  res$col[!is.na(res$padj) & res$padj < 0.05 &
            res$baseMean >= min_bm &
            res$log2FoldChange > 0] <- col_up
  res$col[!is.na(res$padj) & res$padj < 0.05 &
            res$baseMean >= min_bm &
            res$log2FoldChange < 0] <- col_dn
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
                     max.overlaps=Inf, min.segment.length=0,
                     segment.color="grey50", label.size=0.2) +
    theme_bw(base_size=12) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=12)) +
    labs(title=paste0(contrast, " vs Sed"),
         x="log2 fold change", y="-log10(padj)") +
    xlim(c(-xlim_val, xlim_val))

  volcano_list[[contrast]] <- p
}

png(file.path(fig_dir, "Fig1_Volcano_allTimepoints.png"),
    width=14, height=11, units="in", res=300)
grid.arrange(grobs=volcano_list, ncol=2,
             top="Bone marrow sEV miRNA — Exercise vs Sed")
dev.off()
message("Volcano plots saved")

# ── FIGURE 3: TRAJECTORY + HEATMAP (only if significant hits exist) ────────
if (length(all_sig) > 0) {
  message("Making trajectory plots for ", length(all_sig), " miRNAs...")

  traj_list <- list()
  for (contrast in contrasts) {
    res <- res_list[[contrast]]
    for (target in all_sig) {
      hit <- which(res$short_name == target)
      if (length(hit)==0) next
      r <- res[hit[1],]
      traj_list[[length(traj_list)+1]] <- data.frame(
        miRNA=target, Contrast=contrast,
        log2FC=r$log2FoldChange, SE=r$lfcSE,
        padj=ifelse(is.na(r$padj),1,r$padj))
    }
  }
  traj_df <- do.call(rbind, traj_list)
  sed_rows <- data.frame(miRNA=all_sig, Contrast="Sed",
                         log2FC=0, SE=0, padj=1)
  traj_df  <- rbind(sed_rows, traj_df)
  traj_df$Contrast <- factor(traj_df$Contrast,
                              levels=c("Sed","Exer","Exer3","Exer7","Exer18"))
  traj_df$sig <- traj_df$padj < 0.05

  dir_df <- traj_df[traj_df$Contrast != "Sed" & traj_df$sig, ]
  dir_summary <- tapply(dir_df$log2FC, dir_df$miRNA,
                        function(x) ifelse(mean(x)>0,"up","down"))

  traj_plots <- list()
  for (target in sort(all_sig)) {
    df_t <- traj_df[traj_df$miRNA==target,]
    dir  <- if(target %in% names(dir_summary)) dir_summary[[target]] else "up"
    clr  <- ifelse(dir=="up", col_up, col_dn)
    df_t$star <- ifelse(df_t$sig & df_t$Contrast!="Sed","*","")
    star_offset <- ifelse(dir=="up", 0.3, -0.3)
    p <- ggplot(df_t, aes(x=Contrast, y=log2FC, group=1)) +
      geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
      geom_ribbon(aes(ymin=log2FC-SE, ymax=log2FC+SE),
                  fill=clr, alpha=0.15) +
      geom_line(colour=clr, linewidth=1.2) +
      geom_point(colour=clr, fill=clr, size=3, shape=21, stroke=1.2) +
      geom_text(aes(label=star, y=log2FC+star_offset),
                size=5, colour="black") +
      theme_bw(base_size=10) +
      theme(panel.grid.minor=element_blank(),
            plot.title=element_text(face="italic", size=10),
            legend.position="none") +
      labs(title=target, x=NULL, y="log2FC vs Sed")
    traj_plots[[target]] <- p
  }
  ncols <- min(4, length(traj_plots))
  nrows <- ceiling(length(traj_plots)/ncols)
  png(file.path(fig_dir, "Fig2_Trajectory.png"),
      width=16, height=max(5, nrows*3.5), units="in", res=300)
  grid.arrange(grobs=traj_plots, ncol=ncols,
               top="Bone marrow sEV miRNA — temporal trajectories")
  dev.off()
  message("Trajectory plots saved")

  # heatmap
  heat_data <- matrix(NA, nrow=length(all_sig), ncol=length(contrasts),
                      dimnames=list(sort(all_sig), contrasts))
  padj_data <- heat_data
  for (contrast in contrasts) {
    res <- res_list[[contrast]]
    for (target in all_sig) {
      hit <- which(res$short_name==target)
      if (length(hit)==0) next
      heat_data[target,contrast] <- res[hit[1],"log2FoldChange"]
      padj_data[target,contrast] <- res[hit[1],"padj"]
    }
  }
  hm <- melt(heat_data, varnames=c("miRNA","Contrast"), value.name="log2FC")
  pm <- melt(padj_data, varnames=c("miRNA","Contrast"), value.name="padj")
  hm$padj <- pm$padj
  hm$sig_label <- ifelse(!is.na(hm$padj) & hm$padj<0.05,"*","")
  hm$Contrast <- factor(hm$Contrast, levels=contrasts)
  lim <- max(abs(heat_data), na.rm=TRUE)
  phm <- ggplot(hm, aes(x=Contrast, y=miRNA, fill=log2FC)) +
    geom_tile(colour="white", linewidth=0.5) +
    geom_text(aes(label=sig_label), size=5, vjust=0.75) +
    scale_fill_gradient2(low=col_dn, mid="white", high=col_up,
                         midpoint=0, limits=c(-lim,lim),
                         name="log2FC vs Sed") +
    theme_bw(base_size=11) +
    theme(axis.text.y=element_text(face="italic"),
          panel.grid=element_blank()) +
    labs(title="Bone marrow sEV miRNA — log2FC vs Sed",
         subtitle="* padj < 0.05 | baseMean >= 10", x=NULL, y=NULL)
  ggsave(file.path(fig_dir, "Fig3_Heatmap.png"), phm,
         width=8, height=max(4, length(all_sig)*0.4+2), dpi=300)
  message("Heatmap saved")

} else {
  message("No significant hits — trajectory and heatmap skipped")
  message("This confirms the bone marrow null result")
}

message("\nAll figures in: ", file.path(working_dir, fig_dir))
message("=== DONE ===")
