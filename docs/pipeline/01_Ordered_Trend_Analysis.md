# Ordered Trend Analysis — miRNA (Bone Marrow EV)

Tests whether miRNA expression follows a **progression** across the ordered
sequence Sed → Exer (0h) → Exer3 → Exer7 → Exer18, rather than differing in
any single pairwise comparison.

**Script:** `DEG_miRNA_OrderedTrend.R`
**Output:** `DESeq2_results/OrderedTrend/`

```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
Rscript DEG_miRNA_OrderedTrend.R
```

---

## Why this test

Pairwise contrasts ask "is Exer3 different from Sed?" one pair at a time. A
miRNA that shifts consistently *across* the recovery window — never enough to
win any single contrast — would be missed. The trend test looks for that
pattern directly.

Two components are tested:

| Component | Detects |
|-----------|---------|
| **Linear** | steady monotonic progression across the ordered sequence |
| **Quadratic** | a hump or U-shape (e.g. peak mid-recovery, return by 18h) |

---

## Implementation gotcha

DESeq2 **rejects ordered factors** in the design formula:

```
Error in DESeqDataSet(se, design = design, ignoreRank) :
  the design formula contains an ordered factor. The internal steps
  do not work on ordered factors as a formula.
```

Workaround: build the polynomial contrasts manually as numeric columns and pass
those instead.

```r
ordered_levels <- c("Sed","Exer","Exer3","Exer7","Exer18")
metadata$CondPos <- match(metadata$Condition, ordered_levels)   # 1..5
polys <- poly(metadata$CondPos, degree=2)
metadata$trend_linear    <- polys[,1]
metadata$trend_quadratic <- polys[,2]

design <- ~ Isolation_Day + trend_linear + trend_quadratic
```

`poly()` generates orthogonal polynomial contrasts — mathematically identical to
what an ordered factor would produce, but as plain numeric predictors DESeq2
accepts.

**Caveat on spacing:** positions are 1–5, i.e. *equally spaced*. The real
timepoints are not (baseline, 0h, 3h, 7h, 18h). So "linear trend" means linear
across ordered positions, not across true hours. Fine for detecting monotonic
patterns; not a true time-scaled regression.

---

## Results

Filtered: 298 miRNAs × 62 samples.
Design: `~ Isolation_Day + trend_linear + trend_quadratic`

| Component | sig total | sig real (baseMean ≥ 10) |
|-----------|:---:|:---:|
| Linear | 0 | 0 |
| Quadratic | 1 | 0 |

**Linear — top candidates (all padj 0.64):**

| miRNA | baseMean | log2FC | pvalue |
|-------|---:|---:|---:|
| miR-21a-5p | 7124 | +0.45 | 0.028 |
| miR-148a-3p | 5013 | +0.55 | 0.034 |
| miR-223-3p | 1355 | −0.95 | 0.026 |
| miR-26b-5p | 827 | −0.61 | 0.037 |
| let-7g-5p | 630 | −0.94 | 0.007 |
| miR-423-3p | 676 | +1.15 | 0.013 |

All well-expressed, none significant.

**Quadratic — the single hit is an artifact:**

miR-134-5p, baseMean **1.79**, log2FC −19.4, padj 0.016. Near-zero expression
with a physically impossible fold change — the low-count artifact signature.
Zero hits survive the baseMean ≥ 10 filter.

---

## Interpretation

No miRNA shows a significant progression across the recovery sequence, in either
a monotonic or a hump-shaped form. Combined with the null pairwise and pooled
results, this rules out the "consistent small shift across timepoints" scenario
that the trend test exists to catch.

Note miR-21a appears at +0.45 linear and −0.60 quadratic — inconsistent
direction across components, which is itself a sign of noise rather than a
coherent trajectory.

---

*Part of the bone marrow miRNA-EV analysis series.*
