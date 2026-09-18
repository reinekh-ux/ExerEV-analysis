set.seed(910)
suppressMessages({ library(DESeq2) })

working_dir   <- "/home/reinekh/links/scratch/BoneMarrow_miRNA_EV_QC"
rawcount_file <- "count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt"
metadata_file <- "metadata_BoneMarrow.txt"
results_dir   <- "DESeq2_results/BM_manuscript"
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)
setwd(working_dir)

raw  <- round(read.table(rawcount_file, header=TRUE, row.names=1,
                         sep="\t", check.names=FALSE))
meta <- read.table(metadata_file, header=TRUE, sep="\t",
                   stringsAsFactors=FALSE)
rownames(meta) <- meta$Index_Number

raw  <- raw[, colnames(raw) %in% rownames(meta), drop=FALSE]
meta <- meta[match(colnames(raw), rownames(meta)), ]
stopifnot(identical(colnames(raw), rownames(meta)))

meta$Condition     <- factor(meta$Condition,
                             levels=c("Sed","Exer","Exer3","Exer7","Exer18"))
meta$Sex           <- factor(meta$Sex)
meta$Isolation_Day <- factor(meta$Isolation_Day)

# gene filter — identical to muscle
passes <- raw >= 10
keep   <- apply(passes, 1, function(r)
             max(tapply(r, meta$Condition, sum), na.rm=TRUE)) >= 3
raw    <- raw[keep, ]
message("Features after filtering: ", nrow(raw))

meta$Condition <- relevel(meta$Condition, ref="Sed")
dds <- DESeqDataSetFromMatrix(raw, meta,
                              ~ Isolation_Day + Sex + Condition)
dds <- DESeq(dds)

contrasts <- list(
  Exer   = c("Condition","Exer","Sed"),
  Exer3  = c("Condition","Exer3","Sed"),
  Exer7  = c("Condition","Exer7","Sed"),
  Exer18 = c("Condition","Exer18","Sed")
)

targets <- c("miR-21a-5p","miR-128-3p","miR-133b-3p","miR-191-5p",
             "miR-199a-5p","miR-92a-3p","miR-29a-3p")
short <- function(ids) {
  s <- sub(":.*","", sub("^mmu-","", ids))
  sub("\\|.*","", s)
}

message("\n=== KEY miRNAs ACROSS ALL CONTRASTS ===")
message(sprintf("%-10s %-15s %8s %8s %10s",
                "Contrast","miRNA","log2FC","pvalue","padj"))

for (nm in names(contrasts)) {
  res <- as.data.frame(results(dds, contrast=contrasts[[nm]]))
  res$short_name <- short(rownames(res))
  res <- res[order(res$pvalue), ]

  # save full results
  out <- file.path(results_dir,
                   paste0("BM_miRNA_", nm, "vsSed.txt"))
  write.table(res, out, sep="\t", quote=FALSE)

  sig <- sum(!is.na(res$padj) & res$padj < 0.05 & res$baseMean >= 10)
  message("\n--- ", nm, " vs Sed | sig hits (padj<0.05, baseMean>=10): ", sig)

  for (target in targets) {
    hit <- which(res$short_name == target)
    if (length(hit)==0) next
    r <- res[hit[1],]
    message(sprintf("  %-15s %+8.3f %10.4f %10.4f",
                    target, r$log2FoldChange,
                    r$pvalue,
                    ifelse(is.na(r$padj), NA, r$padj)))
  }
}

message("\nFiles saved to: ", file.path(working_dir, results_dir))
message("=== DONE ===")
