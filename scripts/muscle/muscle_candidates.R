set.seed(910)
suppressMessages({
  library(DESeq2); library(ggplot2); library(limma)
})

# ── PATHS ──────────────────────────────────────────────────────────────────
batch1_dir     <- "/lustre09/project/6019280/reinekh/Muscle_data_From_Minh/data/Batch1"
rawcount_file  <- file.path(batch1_dir, "exceRpt_rawcounts_perfect_match_mm10.csv")
metadata_file  <- file.path(batch1_dir, "metadata.csv")
results_dir    <- "/lustre09/project/6019280/reinekh/Muscle_data_From_Minh/results"
out_dir        <- file.path(results_dir, "candidate_boxplots")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

targets <- c("miR-133b-3p", "miR-21a-5p")

# ── LOAD DATA ──────────────────────────────────────────────────────────────
raw <- read.csv(rawcount_file, row.names=1, check.names=FALSE)
raw <- round(raw)

meta <- read.csv(metadata_file, stringsAsFactors=FALSE)
rownames(meta) <- meta$Number_Index   # RPI1, RPI2 ...

# align samples
raw  <- raw[, colnames(raw) %in% rownames(meta), drop=FALSE]
meta <- meta[match(colnames(raw), rownames(meta)), ]
stopifnot(identical(colnames(raw), rownames(meta)))

# ── CONDITION SETUP ─────────────────────────────────────────────────────────
# keep all timepoints including Pre Ex for boxplots
# standardise names to match bone marrow style
meta$Condition <- meta$condition
meta$Condition[meta$Condition == "Ex + 3"]  <- "Exer3"
meta$Condition[meta$Condition == "Ex + 7"]  <- "Exer7"
meta$Condition[meta$Condition == "Ex + 18"] <- "Exer18"
meta$Condition[meta$Condition == "Pre Ex"]  <- "PreEx"

meta <- meta[meta$Condition != "PreEx", ]
raw  <- raw[, rownames(meta), drop=FALSE]
cond_levels <- c("Sed","Exer","Exer3","Exer7","Exer18")
meta$Condition    <- factor(meta$Condition, levels=cond_levels)
meta$Sex          <- factor(meta$sex)
meta$Isolation_Day <- factor(meta$animal_id)   # 1, 2, 3 = isolation days

message("Condition counts:")
print(table(meta$Condition))
message("Isolation day counts:")
print(table(meta$Isolation_Day))
message("Condition x Day cross-tab:")
print(table(meta$Condition, meta$Isolation_Day))

# ── GENE FILTER (same as bone marrow) ──────────────────────────────────────
min_count   <- 10
min_samples <- 3
group_vec   <- meta$Condition
passes      <- raw >= min_count
keep        <- apply(passes, 1, function(r)
                 max(tapply(r, group_vec, sum), na.rm=TRUE)) >= min_samples
raw_filt    <- raw[keep, , drop=FALSE]
message("\nFeatures after filtering: ", nrow(raw_filt))

# ── DESEQ2 (model matches bone marrow + sex added for consistency) ──────────
# Use Sed as reference; PreEx included but Sed is ref
meta_deseq        <- meta
meta_deseq$Condition <- relevel(meta_deseq$Condition, ref="Sed")

dds <- DESeqDataSetFromMatrix(raw_filt, meta_deseq,
                              ~ Isolation_Day + Sex + Condition)
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized=TRUE)
logn <- log2(norm + 1)

# day-adjusted log2 (protecting Condition + Sex)
design_keep <- model.matrix(~ Condition + Sex, data=meta_deseq)
logn_adj    <- removeBatchEffect(logn,
                                 batch  = meta_deseq$Isolation_Day,
                                 design = design_keep)

# ── BOXPLOTS ───────────────────────────────────────────────────────────────
for (target in targets) {

  hit <- grep(paste0("^mmu-", target, ":"), rownames(norm))
  if (length(hit) == 0) { message("NOT FOUND: ", target); next }
  row_id <- rownames(norm)[hit[1]]
  fn     <- gsub("[^A-Za-z0-9]+", "_", target)
  message("\n=== ", target, " ===\nrow: ", row_id)

  df <- data.frame(
    Raw       = as.numeric(norm[row_id, ]),
    Adj       = as.numeric(logn_adj[row_id, ]),
    Condition = meta$Condition,
    Day       = meta$Isolation_Day,
    Sex       = meta$Sex
  )

  # summary tables
  message("Raw medians by condition:")
  print(round(tapply(df$Raw, df$Condition, median), 1))
  message("Day-adjusted log2 medians by condition:")
  print(round(tapply(df$Adj, df$Condition, median), 2))
  message("Raw medians by isolation day:")
  print(round(tapply(df$Raw, df$Day, median), 1))

  # shared theme
  base_theme <- theme_bw(base_size=13) +
    theme(panel.grid.minor=element_blank())

  # A. raw counts coloured by sex
  pA <- ggplot(df, aes(Condition, Raw)) +
    geom_boxplot(outlier.shape=NA, width=0.6, fill="grey95", colour="grey40") +
    geom_jitter(aes(colour=Sex), width=0.15, size=2.5, alpha=0.8) +
    scale_colour_manual(values=c(Female="#d55181", Male="#2a78d6")) +
    base_theme +
    labs(title=paste0(target, " — skeletal muscle sEVs"),
         subtitle="each point = one mouse; DESeq2 size-factor normalized",
         x=NULL, y="normalized count")
  ggsave(file.path(out_dir, paste0(fn, "_boxplot.png")),
         pA, width=9, height=6, dpi=300)

  # B. log10 version
  ggsave(file.path(out_dir, paste0(fn, "_boxplot_log10.png")),
         pA + scale_y_log10() + labs(y="normalized count (log10)"),
         width=9, height=6, dpi=300)

  # C. diagnostic — coloured by isolation day
  pC <- ggplot(df, aes(Condition, Raw)) +
    geom_boxplot(outlier.shape=NA, width=0.6, fill="grey95", colour="grey40") +
    geom_jitter(aes(colour=Day), width=0.15, size=2.5, alpha=0.85) +
    base_theme +
    labs(title=paste0(target, " — coloured by isolation day (diagnostic)"),
         subtitle="DIAGNOSTIC: does shift hold within days or track day?",
         x=NULL, y="normalized count", colour="Isolation\nday")
  ggsave(file.path(out_dir, paste0(fn, "_by_day_colour.png")),
         pC, width=9, height=6, dpi=300)

  # D. faceted by day
  pD <- ggplot(df, aes(Condition, Raw)) +
    geom_boxplot(outlier.shape=NA, width=0.6, fill="grey95", colour="grey40") +
    geom_jitter(width=0.12, size=1.8, alpha=0.8, colour="#2a78d6") +
    facet_wrap(~ Day, nrow=1, labeller=label_both) +
    base_theme +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(title=paste0(target, " — within each isolation day"),
         x=NULL, y="normalized count")
  ggsave(file.path(out_dir, paste0(fn, "_facet_day.png")),
         pD, width=12, height=5, dpi=300)

  # E. day-adjusted (model output, not raw data)
  pE <- ggplot(df, aes(Condition, Adj)) +
    geom_boxplot(outlier.shape=NA, width=0.6, fill="grey95", colour="grey40") +
    geom_jitter(width=0.15, size=2.5, alpha=0.8, colour="#888780") +
    base_theme +
    labs(title=paste0(target, " — log2 counts, isolation day removed"),
         subtitle="MODEL OUTPUT, not raw data (removeBatchEffect; Condition + Sex protected)",
         x=NULL, y="log2 count, day-adjusted")
  ggsave(file.path(out_dir, paste0(fn, "_dayadjusted.png")),
         pE, width=9, height=6, dpi=300)

  message("saved 5 plots for: ", target)
}

message("\nAll plots in: ", out_dir)
message("=== DONE ===")
