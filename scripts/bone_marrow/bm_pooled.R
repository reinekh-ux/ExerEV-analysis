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
meta$Condition     <- factor(meta$Condition,
                             levels=c("Sed","Exer","Exer3","Exer7","Exer18"))
meta$Sex           <- factor(meta$Sex)
meta$Isolation_Day <- factor(meta$Isolation_Day)

meta$ExerciseYesNo <- factor(
  ifelse(meta$Condition == "Sed", "Sed", "Exercise"),
  levels=c("Sed","Exercise"))

message("Group sizes:")
print(table(meta$ExerciseYesNo))
message("Sex distribution:")
print(table(meta$ExerciseYesNo, meta$Sex))
message("Day distribution:")
print(table(meta$ExerciseYesNo, meta$Isolation_Day))

passes <- raw >= 10
keep   <- apply(passes, 1, function(r)
             max(tapply(r, meta$Condition, sum), na.rm=TRUE)) >= 3
raw    <- raw[keep, ]
message("Features after filtering: ", nrow(raw))

dds <- DESeqDataSetFromMatrix(raw, meta,
                              ~ Isolation_Day + Sex + ExerciseYesNo)
dds$ExerciseYesNo <- relevel(dds$ExerciseYesNo, ref="Sed")
dds <- DESeq(dds)

res <- as.data.frame(results(dds, name="ExerciseYesNo_Exercise_vs_Sed"))
res$short_name <- sub(":.*","", sub("^mmu-","",
                   sub("\\|.*","", rownames(res))))
res <- res[order(res$pvalue), ]

out <- file.path(results_dir, "BM_miRNA_PooledExercise_vs_Sed.txt")
write.table(res, out, sep="\t", quote=FALSE)

sig_total <- sum(!is.na(res$padj) & res$padj < 0.05)
sig_real  <- sum(!is.na(res$padj) & res$padj < 0.05 & res$baseMean >= 10)
message("\n=== BONE MARROW POOLED EXERCISE vs SED ===")
message("Total tested: ", nrow(res))
message("padj < 0.05: ", sig_total)
message("padj < 0.05 AND baseMean >= 10 (real hits): ", sig_real)

message("\nTop 15 by pvalue:")
print(head(res[!is.na(res$pvalue),
               c("baseMean","log2FoldChange","pvalue","padj")], 15),
      digits=3)

# p-value histogram
sig_p05 <- sum(!is.na(res$pvalue) & res$pvalue < 0.05)
exp_p05 <- nrow(res) * 0.05
message("\np-value histogram:")
message("p<0.05 observed: ", sig_p05,
        " | expected by chance: ", round(exp_p05, 1))

# key miRNAs
targets <- c("miR-21a-5p","miR-128-3p","miR-133b-3p",
             "miR-191-5p","miR-199a-5p","miR-92a-3p","miR-29a-3p")
message("\n=== KEY miRNAs in pooled test ===")
for (target in targets) {
  hit <- which(res$short_name == target)
  if (length(hit)==0) next
  r <- res[hit[1],]
  message(sprintf("%-15s | baseMean %7.1f | log2FC %+.3f | pvalue %.4f | padj %.4f",
                  target, r$baseMean, r$log2FoldChange,
                  r$pvalue, ifelse(is.na(r$padj), NA, r$padj)))
}
message("=== DONE ===")
