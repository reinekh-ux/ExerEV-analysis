# Muscle sEV miRNA Analysis — Final Scripts & Methods

**Project:** Bone-Muscle-Brain Axis | Nagy Lab | McGill / Douglas Research Centre
**Tissue:** Skeletal muscle sEVs
**Data:** Batch 1 only (Batch 2 excluded — contaminated)
**Date:** July–August 2026
**Cluster:** Alliance Canada Rorqual
**Working directory:** `~/links/projects/rrg-gturecki/reinekh/Muscle_data_From_Minh/`

---

## Overview of What We Did

Starting from Minh Nguyen's original Batch 1 analysis, we reran the muscle miRNA differential expression adding Sex as a covariate (for consistency with the bone marrow model), ran all four per-timepoint contrasts and a pooled analysis, identified significant miRNAs, and generated publication-quality figures.

**Key finding:** 14 miRNAs significantly regulated in skeletal muscle sEVs with exercise (baseMean >= 10, padj < 0.05), with distinct temporal patterns across the post-exercise recovery window (0h, 3h, 7h, 18h).

---

## Input Files

| File | Location on cluster | Description |
|------|--------------------|----|
| Count matrix | `data/Batch1/exceRpt_rawcounts_perfect_match_mm10.csv` | 36 samples x 468 miRNAs, exceRpt perfect match mm10 |
| Metadata | `data/Batch1/metadata.csv` | Sample info: RPI ID, condition, sex, isolation day |
| Minh original results | `results/DESeq2_Batch1_miRNA_Ex7vsSed_AFTERbatchcorrection.txt` | Minh Exer7 vs Sed (no Sex covariate) |

**Metadata columns:**
- Number_Index: sample ID (RPI1-RPI36)
- condition: Sed, Exer, Ex + 3, Ex + 7, Ex + 18, Pre Ex
- sex: Male / Female
- animal_id: isolation day (1, 2, 3)

---

## Key Decisions Made

| Decision | What | Why |
|----------|------|-----|
| Batch 1 only | Excluded Batch 2 | Contaminated |
| Pre Ex excluded | Removed from all analyses | Not relevant to exercise response |
| Sex added to model | ~ Isolation_Day + Sex + Condition | Matches bone marrow model exactly |
| baseMean >= 10 | Expression filter | Same as bone marrow — excludes unreliable low-count features |
| padj < 0.05 | Significance threshold | Standard FDR (Benjamini-Hochberg) |
| Sed as reference | All contrasts vs Sed | Consistent with bone marrow |
| log2FC vs Sed | Trajectory plots | Shows magnitude from baseline — most interpretable |

---

## Script 1: All Contrasts — Per-Timepoint DESeq2

**File:** `muscle_allcontrasts.R`
**Run:** `Rscript muscle_allcontrasts.R`

**What it does:**
Runs DESeq2 on all 30 samples (Pre Ex excluded) with model `~ Isolation_Day + Sex + Condition`. Extracts four separate pairwise contrasts (Exer, Exer3, Exer7, Exer18 vs Sed). Prints summary table for miR-133b-3p and miR-21a-5p. Saves full result files for each contrast.

**Key settings:**
```
design = ~ Isolation_Day + Sex + Condition
reference = Sed
gene filter = >= 10 counts in >= 3 samples of any group
```

**Output files:**
```
results/Muscle_miRNA_ExervsSed_allcontrasts.txt
results/Muscle_miRNA_Exer3vsSed_allcontrasts.txt
results/Muscle_miRNA_Exer7vsSed_allcontrasts.txt
results/Muscle_miRNA_Exer18vsSed_allcontrasts.txt
```

---

## Script 2: Pooled Exercise vs Sed

**File:** `muscle_pooled.R`
**Run:** `Rscript muscle_pooled.R`

**What it does:**
Creates a binary grouping (Sed vs Exercise — all four timepoints pooled) and runs DESeq2 with model `~ Isolation_Day + Sex + ExerciseYesNo`. Maximizes power: 6 Sed vs 24 pooled exercise mice. Reports significant hits and p-value histogram check.

**Key results:**
- 4 significant hits: miR-133b-3p (0.013), miR-199a-5p (0.032), miR-29a-3p (0.046), miR-92a-3p (0.046)
- p-value histogram: 21 observed vs 6.4 expected — real signal confirmed
- miR-21: not significant when pooled (effect diluted across timepoints)

**Output:**
```
results/Muscle_miRNA_PooledExercise_vs_Sed.txt
```

---

## Script 3: Significant Hits Summary CSV

**File:** `muscle_hits_summary.R`
**Run:** `Rscript muscle_hits_summary.R`

**What it does:**
Reads all four contrast result files plus pooled. Extracts log2FC, SE, pvalue, padj for target miRNAs. Saves clean CSV with YES/no significance flag.

**Output:**
```
results/muscle_significant_hits_summary.csv
```

---

## Script 4: Candidate Boxplots

**File:** `muscle_candidates.R`
**Run:** `Rscript muscle_candidates.R`

**What it does:**
For miR-133b-3p and miR-21a-5p, generates 8 plots each (16 total):

| Plot | Description |
|------|-------------|
| _boxplot.png | Boxplot, coloured by sex, linear scale |
| _boxplot_log10.png | Same, log10 y-axis |
| _barplot.png | Bar mean+/-SEM, coloured by sex, linear |
| _barplot_log10.png | Same, log10 y-axis |
| _by_day_colour.png | Diagnostic: coloured by isolation day |
| _facet_day.png | Diagnostic: one panel per isolation day |
| _dayadjusted.png | Day-adjusted boxplot (model output) |
| _dayadjusted_barplot.png | Day-adjusted barplot (model output) |

**Normalization:** DESeq2 size-factor normalized counts.
**Day adjustment:** limma::removeBatchEffect, protecting Condition + Sex.

**Output folder:** `results/candidate_boxplots/`

---

## Script 5: Prism Values Export

**File:** `muscle_prism_values.R`
**Run:** `Rscript muscle_prism_values.R`

**What it does:**
Extracts per-mouse values for miR-133b-3p and miR-21a-5p in four formats. Saves long-format (one row per mouse) and wide-format (conditions as columns, ready for Prism).

**Output files per miRNA:**
```
miR_133b_3p_all_values_prism.csv           (long — all 4 metrics + metadata)
miR_133b_3p_raw_count_wide_prism.csv       (wide — paste into Prism)
miR_133b_3p_normalized_count_wide_prism.csv
miR_133b_3p_log2_normalized_wide_prism.csv
miR_133b_3p_log2_day_adjusted_wide_prism.csv
```
Same 5 files for miR-21a-5p.

---

## Script 6: Manuscript Figures — FINAL (v4)

**File:** `muscle_figures_v4.R`
**Run:** `Rscript muscle_figures_v4.R`

**What it does:**
Automatically identifies all significant miRNAs (padj < 0.05, baseMean >= 10) across all four timepoints, then generates three figures.

### Figure 1 — Volcano plots (2x2 grid)

One panel per timepoint. Settings:
- X axis: log2 fold change vs Sed
- Y axis: -log10(padj)
- Dashed line: padj = 0.05 (y = 1.3)
- Dotted lines: +/-1 log2FC (~2-fold)
- Red = upregulated (padj < 0.05, baseMean >= 10)
- Blue = downregulated (padj < 0.05, baseMean >= 10)
- Grey = NS or baseMean < 10
- Grey dots ABOVE line = significant by padj but excluded (low expression)
- Labels: ggrepel, max.overlaps=Inf (all labels shown)

### Figure 2 — Trajectory plots

One panel per significant miRNA. Settings:
- X axis: timepoints (Sed through Exer18)
- Y axis: log2FC vs Sed (Sed fixed at 0)
- Line: mean log2FC
- Ribbon: +/-1 SE
- Asterisk (*): padj < 0.05 at that timepoint
- Red = predominantly up, Blue = predominantly down

### Figure 3 — Heatmap

All significant miRNAs x all timepoints. Settings:
- Color: log2FC (red = up, blue = down, white = no change)
- Asterisk: padj < 0.05
- miRNAs grouped by direction (up at bottom, down at top)

**Output files:**
```
results/manuscript_figures/Fig1_Volcano_v4.png
results/manuscript_figures/Fig2_Trajectory_v4.png
results/manuscript_figures/Fig3_Heatmap_v4.png
```

---

## Significant miRNAs — Final List (14 total)

### Downregulated with exercise (11):

| miRNA | baseMean | Peak log2FC | Significant timepoints |
|-------|---:|---:|---|
| miR-128-3p | 13.6 | -2.908 | Exer7 |
| miR-191-5p | 115.1 | -3.840 | Exer3, Exer7, Exer18 |
| miR-203-3p | 78.2 | -2.956 | Exer3, Exer7 |
| miR-205-5p | 169.9 | -3.278 | Exer3, Exer7 |
| miR-21a-5p | 136.3 | -1.202 | Exer7 |
| miR-221-3p | 50.2 | -3.059 | Exer3, Exer7 |
| miR-223-3p | 21.1 | -2.934 | Exer7 |
| miR-23b-3p | 319.7 | -2.980 | Exer7 |
| miR-5119 | 14.6 | -3.029 | Exer3, Exer7 |
| miR-652-3p | 10.1 | -3.078 | Exer7 |
| miR-92a-3p | 39.2 | -1.789 | Exer18, Pooled |

### Upregulated with exercise (3):

| miRNA | baseMean | Peak log2FC | Significant timepoints |
|-------|---:|---:|---|
| miR-133b-3p | 10.7 | +3.497 | Exer7, Pooled |
| miR-199a-5p | 221.9 | +1.499 | Exer7, Exer18, Pooled |
| miR-29a-3p | 99.4 | +1.419 | Exer3, Pooled |

---

## Model Comparison — Ours vs Minh (Exer7)

| | Minh (no Sex) | Ours (with Sex) |
|---|:---:|:---:|
| Total significant | 29 | 28 |
| In common | 24 | 24 |
| Unique | 5 | 4 |

**Minh only — sex-confounded, excluded:**
miR-124-3p, miR-15b-5p, miR-200a-3p, miR-200b-3p, miR-451a

**Ours only — revealed after Sex adjustment:**
miR-142a-5p, miR-28a-5p, miR-320-3p, miR-361-5p

**Interpretation:** The 24 hits in common are the most robust — significant regardless of whether Sex is in the model. The 5 Minh-only hits were likely driven by sex differences rather than exercise and are excluded from our analysis.

---

## Figure Legends (copy-paste ready)

**Fig 1:** Each panel shows one timepoint comparison (exercised vs sedentary). Each dot = one miRNA. X-axis = log2 fold change vs sedentary. Y-axis = -log10(adjusted p-value). Horizontal dashed line = padj 0.05 threshold (-log10(0.05) = 1.3). Vertical dotted lines = +/-1 log2FC (~2-fold change). Red = significantly upregulated (padj < 0.05, baseMean >= 10). Blue = significantly downregulated (padj < 0.05, baseMean >= 10). Grey dots above the dashed line passed padj < 0.05 but were excluded due to low expression (baseMean < 10). Differential expression assessed by DESeq2, model ~ Isolation_Day + Sex + Condition, Sed reference, Benjamini-Hochberg FDR correction. n = 6 per group, 30 samples total (Batch 1).

**Fig 2:** Each panel shows one significantly regulated miRNA. X-axis = post-exercise timepoint (Sed = sedentary baseline fixed at 0). Y-axis = log2 fold change vs Sed. Line = mean log2FC at each timepoint. Shaded ribbon = +/-1 standard error. Asterisk (*) = padj < 0.05. Red = predominantly upregulated; Blue = predominantly downregulated. Only miRNAs with padj < 0.05 and baseMean >= 10 in at least one timepoint shown.

**Fig 3:** Rows = miRNAs; columns = post-exercise timepoints vs sedentary. Color = log2FC vs Sed (red = up, blue = down, white = no change). Asterisk = padj < 0.05. Only miRNAs with padj < 0.05 and baseMean >= 10 in at least one timepoint included. miRNAs ordered by direction.

---

## How to Reproduce from Scratch

```bash
ssh reinekh@rorqual.alliancecan.ca
cd ~/links/projects/rrg-gturecki/reinekh/Muscle_data_From_Minh/
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16

Rscript muscle_allcontrasts.R      # step 1: per-timepoint DESeq2
Rscript muscle_pooled.R            # step 2: pooled analysis
Rscript muscle_hits_summary.R      # step 3: clean summary CSV
Rscript muscle_candidates.R        # step 4: boxplots
Rscript muscle_prism_values.R      # step 5: Prism export
Rscript muscle_figures_v4.R        # step 6: manuscript figures
```

---

## Remaining To-Do

- [ ] Check miR-133b in bone marrow data (cross-tissue check)
- [ ] Get miR-133b TargetScan + miRDB CSV (need to download from databases)
- [ ] GO enrichment for miR-133b with TargetScan x miRDB intersection
- [ ] Mechanistic chain notebook for miR-133b
- [ ] Common targets between miR-21 and miR-133b predicted lists
- [ ] Fill in Notebook 07 (bone marrow simple_candidates results)
- [ ] tRNA/piRNA pooled analysis (bone marrow — flagged as high value)

---

*Muscle miRNA-EV analysis | Nagy Lab | McGill / Douglas Research Centre*
*Scripts finalized August 2026*
