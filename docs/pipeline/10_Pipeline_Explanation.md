# Full Pipeline Explanation — From Raw Reads to DESeq2

**Project:** Bone-Muscle-Brain Axis | Nagy Lab | McGill / Douglas Research Centre
**Applies to:** Both skeletal muscle and bone marrow sEV miRNA analyses
**Date:** August 2026

---

## Overview

Both the muscle and bone marrow datasets went through the same pipeline:

```
Raw FASTQ → Trimming (bbduk) → Alignment + Counting (exceRpt/STAR) → Count Matrix Assembly → DESeq2
```

Both used identical exceRpt parameters. Both used totalReadCount. Both are directly comparable.

---

## Step 1 — Raw Sequencing Reads (FASTQ files)

Samples were sequenced on a NovaSeq. Each sample produces a FASTQ file containing millions of short RNA sequences (~16–50 bp), one per read. These are the raw data — just sequences with quality scores, no biological meaning attached yet.

---

## Step 2 — Trimming (bbduk)

Adapter sequences and low-quality bases are removed. Parameters used:
- Minimum read length after trimming: **16 bp**
- Reads shorter than 16 bp are discarded

This ensures only real biological sequences go into alignment.

**Output:** trimmed FASTQ files in `trimmed/`

---

## Step 3 — exceRpt Alignment (the key step)

exceRpt is a small RNA pipeline that takes each trimmed read and identifies what RNA it came from. It checks for ribosomal RNA (discards those), then aligns reads to known small RNA databases using **STAR**.

### Alignment parameters used (identical for muscle and bone marrow):

| Parameter | Value | Meaning |
|-----------|-------|---------|
| STAR_outFilterMismatchNmax | 0 | Zero mismatches allowed — perfect match alignment |
| STAR_alignEndsType | EndToEnd | Full read must align, not just part of it |
| STAR_outFilterMatchNmin | 16 | At least 16 bases must match |
| MIN_READ_LENGTH | 16 | Minimum read length 16 bp |
| ADAPTER_SEQ | none | No adapter trimming (already done by bbduk) |
| MAIN_ORGANISM_GENOME_ID | mm10 | Mouse genome reference |

> **What "perfect match" means:** The alignment is perfect — zero mismatches to the genome. This is what "perfect match" refers to in the muscle filename `exceRpt_rawcounts_perfect_match_mm10.csv`. It does NOT mean reads mapped uniquely — it means they aligned with no errors.

**Output per sample:** `readCounts_miRNAmature_sense.txt` — a table with every miRNA and its read counts, with four columns:

```
ReferenceID | uniqueReadCount | totalReadCount | multimapAdjustedReadCount | multimapAdjustedBarcodeCount
```

---

## Step 4 — Understanding uniqueReadCount vs totalReadCount

After alignment, some reads map to **one specific miRNA** (unambiguous) and some map to **multiple miRNAs** (multimappers).

### Example — miR-486a and miR-486b

These two miRNAs have nearly identical sequences. A read from either one looks almost the same. STAR will align that read to both miR-486a AND miR-486b — because it matches both perfectly with zero mismatches.

| Column | What it counts | miR-486 example |
|--------|---------------|-----------------|
| **uniqueReadCount** | Reads mapping to ONE miRNA only — no ambiguity | Read maps to both 486a and 486b → counted in **neither** |
| **totalReadCount** | All reads mapping to a miRNA, including multimappers | Read maps to both 486a and 486b → counted in **both** |
| **multimapAdjustedReadCount** | Multimappers fractionally distributed by relative abundance | 0.75 to 486a, 0.25 to 486b |

### Why bone marrow miR-486 has totalReadCount = 15,495 but uniqueReadCount = 26:

Bone marrow produces blood cells. Red blood cells and platelets are packed with miR-486. So bone marrow sEVs have enormous amounts of miR-486. Because miR-486a and miR-486b are so similar, almost all 15,495 reads map to both simultaneously — so they're all counted in totalReadCount for both, but none can be uniquely assigned, giving uniqueReadCount of only 26.

Muscle has much less miR-486 biologically, so counts are in the low hundreds — the multimapping inflation exists but is proportionally smaller.

**This is biology, not a pipeline artifact.** Both tissues used the same pipeline and counting column.

---

## Step 5 — Count Matrix Assembly

The `assemble_miRNA_count_table.R` script reads `readCounts_miRNAmature_sense.txt` from each sample folder and extracts the **totalReadCount** column, combining them into one matrix:

```
            D1_1    D1_2    D1_3   ...  D6_12
miR-21a-5p  5427    4832    6103   ...   3891
miR-92a-3p 10080    9234   11200   ...   8654
miR-486a   15495   14200   13800   ...  16200
```

**Bone marrow output:** `BoneMarrow_miRNA_TotalRawCount_matrix.txt` — 697 miRNAs × 62 samples
**Muscle output:** `exceRpt_rawcounts_perfect_match_mm10.csv` — 468 miRNAs × 36 samples

Both use totalReadCount. Both assembled the same way. Both are directly comparable.

---

## Step 6 — DESeq2 Differential Expression

DESeq2 takes raw integer counts and performs four steps:

### 6.1 Size factor estimation
Corrects for library depth differences between samples. Some samples were sequenced more deeply than others — size factors normalize for this so counts are comparable across samples.

### 6.2 Dispersion estimation
Models the biological variability of each miRNA across samples. miRNAs with high variability need a higher bar for significance than stable ones.

### 6.3 Statistical testing
Uses a negative binomial model to test whether each miRNA's counts differ significantly between conditions (e.g. Exer7 vs Sed). The model `~ Isolation_Day + Sex + Condition` accounts for isolation day and sex effects before testing the exercise effect — so those sources of variance are removed and don't contaminate the comparison.

### 6.4 FDR correction
Benjamini-Hochberg correction is applied across all tested miRNAs to control the false discovery rate. This adjusts p-values upward to account for the fact that testing hundreds of miRNAs simultaneously means some will look significant by chance.

**Output:** for each miRNA in each contrast:
- `baseMean` — average expression across all samples
- `log2FoldChange` — magnitude and direction of change (positive = up, negative = down)
- `pvalue` — raw statistical significance
- `padj` — FDR-corrected p-value

---

## Step 7 — Filtering and Interpretation

**baseMean ≥ 10** — removes miRNAs with too few counts to be reliably measured. A fold change estimate for a miRNA with 2 counts in one group and 4 in another is meaningless — any small noise creates a huge apparent effect.

**padj < 0.05** — after FDR correction, only miRNAs where the exercise effect is unlikely to be a false positive are reported.

---

## Why the Two Datasets Are Directly Comparable

| Feature | Muscle | Bone marrow | Same? |
|---------|--------|-------------|-------|
| Sequencer | NovaSeq | NovaSeq | ✓ |
| Trimming | bbduk, min 16bp | bbduk, min 16bp | ✓ |
| exceRpt version | v4.6.3 | v4.6.3 | ✓ |
| Mismatches | 0 | 0 | ✓ |
| Alignment type | EndToEnd | EndToEnd | ✓ |
| Reference genome | mm10 | mm10 | ✓ |
| Count column | totalReadCount | totalReadCount | ✓ |
| DESeq2 model | ~ Isolation_Day + Sex + Condition | ~ Isolation_Day + Sex + Condition | ✓ |
| Reference condition | Sed | Sed | ✓ |
| Gene filter | baseMean ≥ 10 | baseMean ≥ 10 | ✓ |
| FDR threshold | padj < 0.05 | padj < 0.05 | ✓ |

The count scale differs between tissues (bone marrow higher overall due to miR-486 and hematopoietic biology) but DESeq2's size factor normalization handles this — it normalizes each sample regardless of total count depth.

**The log2FoldChange values are what matter for cross-tissue comparison**, and those are computed entirely within each tissue — unaffected by the absolute count scale difference between tissues.

---

## Metadata Column Mapping

### Muscle metadata (metadata.csv):

| Metadata column | Values | Used in script as |
|----------------|--------|-------------------|
| `Number_Index` | RPI1–RPI36 | rownames(meta) |
| `sex` | Male, Female | `meta$Sex <- factor(meta$sex)` |
| `animal_id` | 1, 2, 3 | `meta$Isolation_Day <- factor(meta$animal_id)` |
| `condition` | Sed, Exer, Ex + 3, Ex + 7, Ex + 18, Pre Ex | `meta$Condition <- meta$condition` (then renamed) |

### Bone marrow metadata (metadata_BoneMarrow.txt):

| Metadata column | Values | Used in script as |
|----------------|--------|-------------------|
| `Index_Number` | D1_1–D6_12 | rownames(meta) |
| `Sex` | Male, Female | `meta$Sex <- factor(meta$Sex)` |
| `Isolation_Day` | 1–6 | `meta$Isolation_Day <- factor(meta$Isolation_Day)` |
| `Condition` | Sed, Exer, Exer3, Exer7, Exer18 | `meta$Condition <- factor(meta$Condition)` |

---

## One Sentence Summary

Raw reads → trim adapters (bbduk, min 16bp) → align to mm10 with zero mismatches (exceRpt/STAR) → count how many reads hit each miRNA (totalReadCount) → combine 62 or 36 samples into one matrix → normalize for library depth and test for differential expression with DESeq2 controlling for isolation day and sex → filter for reliable (baseMean ≥ 10) and significant (padj < 0.05) results.

---

## Methods Text (copy-paste ready)

*"Small RNA libraries were processed using exceRpt v4.6.3 with STAR alignment to the mm10 mouse genome reference. Alignment parameters enforced perfect sequence matching (STAR_outFilterMismatchNmax=0, STAR_alignEndsType=EndToEnd, STAR_outFilterMatchNmin=16, MIN_READ_LENGTH=16). Total read counts (totalReadCount) per miRNA per sample were assembled into count matrices for each tissue. Differential expression analysis was performed using DESeq2, with the model ~ Isolation_Day + Sex + Condition, sedentary animals as the reference group, and Benjamini-Hochberg FDR correction. miRNAs with baseMean < 10 were excluded from analysis. Significant differential expression was defined as padj < 0.05."*

---

*Pipeline documentation | Nagy Lab | McGill / Douglas Research Centre | August 2026*
