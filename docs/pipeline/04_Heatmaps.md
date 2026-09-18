# Heatmaps — and the circularity problem

Four heatmap variants, three of which are interpretable and one of which is a
demonstration of why the standard recipe fails on null data.

**Scripts:** `miRNA_heatmaps.R`, `miRNA_heatmaps_condordered.R`
**Output:** `DESeq2_results/heatmaps/`, `DESeq2_results/heatmaps_condordered/`

```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
Rscript miRNA_heatmaps_condordered.R
```

---

## The circularity trap

The conventional DE heatmap recipe is: take the top N genes by padj, plot their
scaled expression, cluster. On data with real signal this produces clean group
blocks.

**It also produces clean blocks on pure noise.** You select genes *for* differing
between groups, then plot them and observe that they differ between groups. The
result is guaranteed by the selection, not by biology.

Row-scaling compounds it. Converting each row to z-scores stretches whatever
range exists to fill the full colour spectrum, so a 15% difference renders as
dramatically as a 5-fold one. The colour carries no information about effect size.

Given a flat p-value distribution, a top-by-padj heatmap on this dataset would
show convincing structure representing nothing.

---

## The four variants

| # | Heatmap | Gene selection | Circular? |
|---|---------|----------------|-----------|
| 1 | Sample-distance | none — all 298 miRNAs | No |
| 2 | Top-variable | highest variance across all samples | No |
| 3 | Top-LRT | lowest LRT padj | **Yes** |
| 4 | Top-variable, condition-ordered | highest variance | No |

### 1. Sample-distance

62 × 62 grid, every sample against every other. Each cell is overall
dissimilarity computed across **all** 298 miRNAs. No gene selection, so no
opportunity to cherry-pick.

*Looking for:* dark blocks along the diagonal aligning with the condition
annotation bar, indicating exercised samples resemble each other.

*Limitation:* a real effect in 3 miRNAs out of 298 would be invisible — the
other 295 dominate the distance.

### 2. Top-variable genes

30 miRNAs with the highest variance across all 62 samples. **Variance is
computed ignoring group membership** — condition plays no role in selection.

*Looking for:* if exercise is among the things making miRNAs vary, condition
should emerge in the column pattern unprompted. If the pattern tracks
Isolation_Day instead, day-to-day technical variation is the driver — consistent
with what variancePartition showed.

*Why it's fair:* genes chosen blind to the hypothesis. Structure appearing here
would be a real observation.

### 3. Top-LRT genes — circular, labelled as such

Structurally identical to #2, but the 30 miRNAs are those with the lowest LRT
padj — i.e. selected *because* they looked condition-associated.

Filename and title both carry `CIRCULAR`. It exists to sit beside #2: same
samples, same ordering, same scaling, differing only in the gene-selection step.
The visual difference between them demonstrates the mechanism more effectively
than any verbal explanation.

### 4. Condition-ordered top-variable

Columns ordered Sed → Exer → Exer3 → Exer7 → Exer18 with gaps between blocks,
clustering disabled on columns. Genes still selected by variance alone.

```r
ord <- order(metadata$Condition, metadata$Isolation_Day)
gaps <- cumsum(table(meta_ord$Condition)); gaps <- gaps[-length(gaps)]
pheatmap(mat, scale="row", cluster_cols=FALSE, gaps_col=gaps, ...)
```

This is the presentation figure. Grouped by condition so the comparison is
visible, genes selected fairly so what it shows is real.

**Caption:** *Top 30 most variable miRNAs (selected independently of condition),
samples grouped by condition. No condition-associated structure is apparent.*

---

## Implementation notes

- **VST with `blind=TRUE`** — the transformation must not use the design, or it
  leaks condition information into the visualisation.
- **Row labels** — exceRpt IDs are long (`mmu-miR-21a-5p:MIMAT0000530:Mus:...`);
  shortened with `sub(":.*", "", sub("^mmu-", "", ids))`.
- **Annotation bars** — Condition and Isolation_Day both shown, so day-driven
  structure is distinguishable from condition-driven structure.

---

## Which to use where

| Purpose | Figure |
|---------|--------|
| Presentation slide | #4 (condition-ordered top-variable) |
| Whole-profile check | #1 (sample-distance) |
| Demonstrating circularity | #2 beside #3 |
| Anything requiring evidence | never #3 alone |

---

*Part of the bone marrow miRNA-EV analysis series.*
