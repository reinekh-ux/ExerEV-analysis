# Simple Candidate Ranking — abundance filter + effect size

A stripped-down ranking: **baseMean ≥ 50, sorted by largest median shift from
Sed.** No power calculation, no consistency filter, no p-values.

**Script:** `simple_candidates.R`
**Output:** `DESeq2_results/simple_candidates.txt`

```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
Rscript simple_candidates.R
```

---

## What it computes

For each miRNA passing the gene filter:

1. DESeq2 size-factor normalised counts
2. Median per condition
3. The exercise group whose median is furthest from Sed
4. That shift as both % change and log2FC

```r
med <- t(apply(norm, 1, function(x) tapply(x, meta$Condition, median)))
sed <- med[, "Sed"]
ex  <- med[, c("Exer","Exer3","Exer7","Exer18")]
idx  <- apply(abs(ex - sed), 1, which.max)
peak <- ex[cbind(seq_len(nrow(ex)), idx)]
pct  <- ifelse(sed > 0, (peak/sed - 1) * 100, NA)
```

Then filtered to `baseMean >= 50` and sorted by `abs(change_pct)`.

**Output columns:** `miRNA`, `baseMean`, `Sed_median`, `peak_group`,
`peak_median`, `change_pct`, `log2FC`.

The `peak_group` column matters — it shows *which* timepoint drove the shift,
which distinguishes a coherent response from a single aberrant group.

---

## The NA bug, fixed

The earlier `qpcr_candidate_screen.R` returned `NA` for `peak_change_pct` on
miR-21a-5p and miR-25-3p. The `ifelse(sed > 0, ...)` guard resolves it — both now
receive proper values.

---

## Why the abundance filter replaces the p-value

Ranking purely by fold change puts near-zero miRNAs on top: a feature with 0
counts in most mice and 3 in two of them produces an enormous, meaningless fold
change. This is what generated the earlier padj artifacts (miR-1843a at −24.8,
baseMean 0.98).

The p-value was partly encoding *measurement reliability*. Removing it entirely
means removing the thing that separates a real 24% shift from a noise-driven
1000% one.

`baseMean ≥ 50` does that job more transparently — the rule is "only miRNAs I can
measure," not "only miRNAs that passed a test." Defensible, and not circular in
the way padj-ranking is.

---

## Interpretation caution

With a flat p-value distribution, the top of an abundance-filtered fold-change
list is **also** largely noise — just larger noise. The ranking is dominated by
whichever well-expressed miRNAs had the widest scatter.

A bigger apparent effect in a null dataset indicates more noise, not more signal.
This is the same trap as let-7e in different clothing: rank by padj → let-7e;
rank by fold change → something else. Both rankings are drawn from a distribution
the histogram shows to be flat.

Cross-reference against the **detection floor** (notebook 06) before treating any
entry as a candidate. A 58% shift on a miRNA whose MDE is 150% is inside the
noise band.

An alternative sort worth running: order the same filtered list by **baseMean**
rather than by change. That surfaces the miRNAs where measurement is most
reliable and noise is tightest.

---

## Results

*[Pending — paste the printed table here after running.]*

| miRNA | baseMean | Sed_median | peak_group | peak_median | change_pct | log2FC |
|-------|---:|---:|---|---:|---:|---:|
| | | | | | | |

**miR-21a-5p (corrected value):**

*[Pending]*

---

*Part of the bone marrow miRNA-EV analysis series.*
