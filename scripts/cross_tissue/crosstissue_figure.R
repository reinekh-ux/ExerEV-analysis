set.seed(910)
suppressMessages({ library(ggplot2); library(gridExtra) })

# ── DATA — from targeted analysis output ────────────────────────────────────
df <- data.frame(
  miRNA = c("miR-21a-5p","miR-128-3p","miR-652-3p","miR-191-5p",
            "miR-29a-3p","miR-5119","miR-221-3p","miR-133b-3p",
            "miR-205-5p","miR-199a-5p","miR-203-3p","miR-92a-3p",
            "miR-23b-3p","miR-223-3p"),
  log2FC_muscle = c(-0.904,-1.289,-1.289,-3.840,
                     0.965,-3.029,-3.059, 3.186,
                    -3.278, 1.182,-2.956,-1.289,
                    -2.980,-2.934),
  log2FC_BM    = c( 0.197, 0.292, 0.249,-0.298,
                    0.103, 0.436, 0.087,-0.441,
                    0.293, 0.075,-0.213,-0.047,
                    0.046,-0.025),
  pvalue_BM    = c(0.0084,0.0479,0.0778,0.1343,
                   0.1915,0.2280,0.3771,0.4195,
                   0.5102,0.5342,0.6271,0.7378,
                   0.8111,0.8678),
  padj_14      = c(0.1177,0.3353,0.3630,0.4701,
                   0.5320,0.5320,0.7341,0.7341,
                   0.7479,0.7479,0.7981,0.8608,
                   0.8678,0.8678),
  padj_muscle  = c(0.037,0.046,0.046,0.046,
                   0.046,0.046,0.046,0.013,
                   0.046,0.032,0.046,0.046,
                   0.046,0.046)
)

# direction agreement
df$direction <- ifelse(sign(df$log2FC_muscle) == sign(df$log2FC_BM),
                       "Same direction", "Opposite direction")

# significance stars for muscle
df$sig_muscle <- ifelse(df$padj_muscle < 0.05, "*", "")
# trend marker for BM
df$trend_BM <- ifelse(df$pvalue_BM < 0.05, "†", "")

# order by muscle log2FC
df$miRNA <- factor(df$miRNA,
                   levels=df$miRNA[order(df$log2FC_muscle)])

col_same <- "#9B59B6"   # purple — same direction
col_opp  <- "#E67E22"   # orange — opposite direction
col_mu_up  <- "#E41A1C" # red — muscle up
col_mu_dn  <- "#377EB8" # blue — muscle down

# ── FIGURE A: side-by-side dot plot ────────────────────────────────────────
# melt to long format
long <- rbind(
  data.frame(miRNA=df$miRNA, Tissue="Skeletal muscle",
             log2FC=df$log2FC_muscle,
             sig=df$sig_muscle,
             direction=df$direction),
  data.frame(miRNA=df$miRNA, Tissue="Bone marrow",
             log2FC=df$log2FC_BM,
             sig=df$trend_BM,
             direction=df$direction)
)
long$Tissue <- factor(long$Tissue,
                      levels=c("Skeletal muscle","Bone marrow"))

pA <- ggplot(long, aes(x=log2FC, y=miRNA, colour=Tissue, shape=Tissue)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey60") +
  geom_point(size=3.5, alpha=0.9) +
  geom_text(aes(label=sig),
            hjust=-0.3, vjust=0.5, size=5, colour="black") +
  scale_colour_manual(values=c("Skeletal muscle"="#E41A1C",
                                "Bone marrow"="#377EB8")) +
  scale_shape_manual(values=c("Skeletal muscle"=16,
                               "Bone marrow"=17)) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(),
        axis.text.y=element_text(face="italic", size=10),
        legend.position="top") +
  labs(title="Cross-tissue miRNA log2FC comparison",
       subtitle="* padj<0.05 in muscle | † raw p<0.05 in bone marrow (targeted analysis)",
       x="log2 fold change vs Sed", y=NULL,
       colour=NULL, shape=NULL)

# ── FIGURE B: paired slope plot ─────────────────────────────────────────────
pB <- ggplot(df, aes(colour=direction)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
  geom_segment(aes(x=1, xend=2,
                   y=log2FC_muscle, yend=log2FC_BM),
               linewidth=0.8, alpha=0.7) +
  geom_point(aes(x=1, y=log2FC_muscle), size=3) +
  geom_point(aes(x=2, y=log2FC_BM), size=3, shape=17) +
  geom_text(aes(x=1, y=log2FC_muscle,
                label=paste0(as.character(miRNA), sig_muscle)),
            hjust=1.15, size=2.8, fontface="italic",
            colour="grey30") +
  scale_colour_manual(values=c("Same direction"=col_same,
                                "Opposite direction"=col_opp),
                      name="") +
  scale_x_continuous(breaks=c(1,2),
                     labels=c("Skeletal\nmuscle","Bone\nmarrow"),
                     limits=c(0.5, 2.8)) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(),
        legend.position="top",
        axis.text.x=element_text(size=11, face="bold")) +
  labs(title="Exercise-induced log2FC: muscle vs bone marrow",
       subtitle="Purple = same direction | Orange = opposite direction | * sig in muscle",
       x=NULL, y="log2 fold change vs Sed")

# ── FIGURE C: heatmap both tissues ─────────────────────────────────────────
heat <- rbind(
  data.frame(miRNA=df$miRNA, Tissue="Skeletal muscle",
             log2FC=df$log2FC_muscle, direction=df$direction),
  data.frame(miRNA=df$miRNA, Tissue="Bone marrow",
             log2FC=df$log2FC_BM, direction=df$direction)
)
heat$Tissue <- factor(heat$Tissue,
                      levels=c("Skeletal muscle","Bone marrow"))
lim <- max(abs(c(df$log2FC_muscle, df$log2FC_BM)))

pC <- ggplot(heat, aes(x=Tissue, y=miRNA, fill=log2FC)) +
  geom_tile(colour="white", linewidth=0.5) +
  scale_fill_gradient2(low="#377EB8", mid="white", high="#E41A1C",
                       midpoint=0, limits=c(-lim, lim),
                       name="log2FC vs Sed") +
  theme_bw(base_size=11) +
  theme(axis.text.y=element_text(face="italic", size=9),
        axis.text.x=element_text(face="bold"),
        panel.grid=element_blank()) +
  labs(title="Cross-tissue heatmap",
       subtitle="14 muscle-significant miRNAs",
       x=NULL, y=NULL)

# ── SAVE ────────────────────────────────────────────────────────────────────
results_dir <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC/DESeq2_results/BM_manuscript"
fig_dir <- file.path(results_dir, "figures")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

ggsave(file.path(fig_dir, "CrossTissue_dotplot.png"),
       pA, width=9, height=8, dpi=300)
message("Dot plot saved")

ggsave(file.path(fig_dir, "CrossTissue_slopeplot.png"),
       pB, width=7, height=9, dpi=300)
message("Slope plot saved")

ggsave(file.path(fig_dir, "CrossTissue_heatmap.png"),
       pC, width=5, height=7, dpi=300)
message("Heatmap saved")

# summary stats
n_opp  <- sum(df$direction == "Opposite direction")
n_same <- sum(df$direction == "Same direction")
message("\n=== DIRECTION SUMMARY ===")
message("Same direction: ", n_same, "/14")
message("Opposite direction: ", n_opp, "/14")
message("Expected by chance: 7/14")
message("\nOpposite direction miRNAs:")
print(df[df$direction=="Opposite direction",
         c("miRNA","log2FC_muscle","log2FC_BM")],
      row.names=FALSE)

message("\nAll figures in: ", fig_dir)
message("=== DONE ===")
