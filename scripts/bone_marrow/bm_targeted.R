set.seed(910)
suppressMessages({ library(DESeq2) })

working_dir   <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
rawcount_file <- "count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt"
metadata_file <- "metadata_BoneMarrow.txt"
results_dir   <- "DESeq2_results/BM_manuscript"
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)
setwd(working_dir)

# the 14 significant muscle miRNAs — pre-specified hypothesis
muscle_sig <- c("miR-133b-3p","miR-191-5p","miR-199a-5p","miR-203-3p",
                "miR-205-5p","miR-21a-5p","miR-221-3p","miR-223-3p",
                "miR-23b-3p","miR-29a-3p","miR-5119","miR-652-3p",
                "miR-92a-3p","miR-128-3p")

short <- function(ids) {
  s <- sub(":.*","", sub("^mmu-","", ids))
  sub("\\|.*","", s)
}

# ── LOAD DATA ───────────────────────────────────────────────────────────────
raw  <- round(read.table(rawcount_file, header=TRUE, row.names=1,
                         sep="\t", check.names=FALSE))
meta <- read.table(metadata_file, header=TRUE, sep="\t",
                   stringsAsFactors=FALSE)
rownames(meta) <- meta$Index_Number
raw  <- raw[, colnames(raw) %in% rownames(meta), drop=FALSE]
meta <- meta[match(colnames(raw), rownames(meta)), ]
meta$Condition     <- factor(meta$Condition,
                             levels=c("Sed","Exer","Exer3","Exer7","Exer18"))
meta$Sex           <- factor(meta$Sex)
meta$Isolation_Day <- factor(meta$Isolation_Day)

# gene filter
passes <- raw >= 10
keep   <- apply(passes, 1, function(r)
             max(tapply(r, meta$Condition, sum), na.rm=TRUE)) >= 3
raw    <- raw[keep, ]

# pooled grouping
meta$ExerciseYesNo <- factor(
  ifelse(meta$Condition == "Sed", "Sed", "Exercise"),
  levels=c("Sed","Exercise"))

# ── RUN DESEQ2 ─────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(raw, meta,
                              ~ Isolation_Day + Sex + ExerciseYesNo)
dds$ExerciseYesNo <- relevel(dds$ExerciseYesNo, ref="Sed")
dds <- DESeq(dds)

res <- as.data.frame(results(dds, name="ExerciseYesNo_Exercise_vs_Sed"))
res$short_name <- short(rownames(res))

# ── EXTRACT TARGETED MIRNAS ─────────────────────────────────────────────────
targeted <- data.frame()
for (target in muscle_sig) {
  hit <- which(res$short_name == target)
  if (length(hit)==0) {
    message("NOT FOUND in BM data: ", target)
    next
  }
  r <- res[hit[1],]
  targeted <- rbind(targeted, data.frame(
    miRNA      = target,
    baseMean   = round(r$baseMean, 1),
    log2FC_BM  = round(r$log2FoldChange, 3),
    pvalue     = signif(r$pvalue, 4),
    padj_298   = signif(r$padj, 4)   # original padj across all 298
  ))
}

# ── APPLY TARGETED BH CORRECTION (across 14 tests only) ────────────────────
targeted <- targeted[order(targeted$pvalue), ]
n <- nrow(targeted)
targeted$padj_14 <- p.adjust(targeted$pvalue, method="BH")
targeted$sig_14  <- ifelse(targeted$padj_14 < 0.05, "YES", "no")
targeted$sig_298 <- ifelse(!is.na(targeted$padj_298) &
                             targeted$padj_298 < 0.05, "YES", "no")

# also add muscle direction for comparison
muscle_dir <- c(
  "miR-133b-3p"="+3.19","miR-191-5p"="-3.84","miR-199a-5p"="+1.18",
  "miR-203-3p"="-2.96","miR-205-5p"="-3.28","miR-21a-5p"="-0.90",
  "miR-221-3p"="-3.06","miR-223-3p"="-2.93","miR-23b-3p"="-2.98",
  "miR-29a-3p"="+0.97","miR-5119"="-1.29","miR-652-3p"="-1.29",
  "miR-92a-3p"="-1.29","miR-128-3p"="-0.90"
)
targeted$log2FC_muscle <- muscle_dir[targeted$miRNA]

message("\n=== TARGETED ANALYSIS: 14 muscle miRNAs in bone marrow ===")
message("Pre-specified hypothesis: based on muscle significant hits")
message("BH correction applied across 14 tests (not 298)\n")
message(sprintf("%-15s %8s %8s %8s %8s %8s %6s",
                "miRNA","baseMean","BM_log2FC","Musc_log2FC",
                "pvalue","padj_14","sig?"))
message(strrep("-", 75))
for (i in seq_len(nrow(targeted))) {
  r <- targeted[i,]
  message(sprintf("%-15s %8.1f %+8.3f %11s %8.4f %8.4f %6s",
                  r$miRNA, r$baseMean, r$log2FC_BM,
                  r$log2FC_muscle, r$pvalue,
                  r$padj_14, r$sig_14))
}

message("\nNOTE: padj_14 = BH correction across these 14 miRNAs only")
message("      padj_298 = original BH correction across all 298 BM miRNAs")
message("      These miRNAs were pre-selected based on muscle significance")
message("      Frame as: cross-tissue replication attempt, not discovery")

# save
out <- file.path(results_dir, "BM_targeted_muscle14_analysis.csv")
write.csv(targeted[, c("miRNA","baseMean","log2FC_BM","log2FC_muscle",
                        "pvalue","padj_14","padj_298","sig_14","sig_298")],
          out, row.names=FALSE)
message("\nSaved: ", file.path(working_dir, out))
message("=== DONE ===")
