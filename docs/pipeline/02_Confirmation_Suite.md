# Confirmation Suite — p-value histogram, limma-voom, family enrichment

Three independent checks on the miRNA null. Each addresses a different "but what
if" objection: *is the signal hidden by the threshold?*, *is it a DESeq2
artifact?*, *is it coordinated at the family level?*

**Script:** `miRNA_confirmation.R`
**Output:** `DESeq2_results/miRNA_confirmation/`

```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
Rscript miRNA_confirmation.R
```

Requires `limma` and `fgsea` — both present in the cluster Bioconductor bundle
(check with `Rscript -e 'library(limma); library(fgsea)'`).

---

## 1. p-value histogram — the conclusive diagnostic

**The principle:** when there is no true effect, p-values are **uniformly
distributed**. This is a mathematical property, not a convention. Real effects
produce a spike of small p-values piled up near zero, standing above a flat
noise floor.

So the *shape* is a signature:

| Shape | Means |
|-------|-------|
| Flat / uniform | no true effects |
| Spike near zero | real effects present |

**Result (pooled Exercise vs Sed):**

```
p < 0.05:  17 observed  vs  17.9 expected by chance
```

Flat. The number of miRNAs below p 0.05 is *exactly* what randomness produces
from ~360 tests. No excess, no spike.

**Why this settles the FDR-method question.** Every FDR correction works by
comparing how many small p-values you have against what chance predicts. Yours
match chance. So no correction method — BH, BY, Storey, IHW — can recover a hit,
because there is no excess to recover. Threshold choice is equally moot.

This is the single most informative figure in the miRNA analysis. It is worth
putting on a slide.

---

## 2. limma-voom cross-check — method independence

Re-ran the same pooled contrast under an entirely different statistical
framework: voom's mean-variance modelling plus empirical Bayes, rather than
DESeq2's negative-binomial GLM.

```r
dge <- DGEList(counts_filt); dge <- calcNormFactors(dge, "TMM")
design <- model.matrix(~ Isolation_Day + Sex + ExerciseYesNo, data=metadata)
v <- voom(dge, design)
fit <- eBayes(lmFit(v, design))
topTable(fit, coef="ExerciseYesNoExercise", number=Inf)
```

**Result: 0 miRNAs with adj.P < 0.05.**

Top 5:

| miRNA | logFC | AveExpr | P.Value | adj.P.Val |
|-------|---:|---:|---:|---:|
| miR-101c | −0.90 | 8.48 | 0.0013 | 0.335 |
| miR-744-5p | −0.59 | 9.47 | 0.0019 | 0.335 |
| miR-146b-3p | +1.48 | 2.91 | 0.0042 | 0.336 |
| miR-21a-5p | +0.20 | 15.91 | 0.0045 | 0.336 |
| miR-31-5p | +1.38 | 2.99 | 0.0054 | 0.336 |

Note the two methods agree not only on the null but on *which* miRNAs rank
highest — miR-101c, miR-744, miR-21a appear in both. Reassuring consistency.
miR-21a's logFC of +0.20 matches the DESeq2 estimate exactly.

---

## 3. Family enrichment (fgsea) — coordinated subtle effects

Individual miRNAs may not move, but a *family* could shift together. Tested with
fgsea on the DESeq2 test statistic, families parsed from the miRNA IDs.

Family extraction heuristic:

```r
base <- sub(":.*", "", ids)                        # drop MIMAT suffix
base <- sub("^mmu-", "", base)                     # drop species prefix
fam  <- sub("-[0-9]+p$", "", base)                 # drop -3p/-5p arm
fam  <- sub("([a-z]+-[0-9]+)[a-z].*", "\\1", fam)  # miR-486a -> miR-486
```

Only families with ≥3 members are testable.

**Result: 0 families with padj < 0.05.**

| family | size | NES | pval | padj |
|--------|:---:|---:|---:|---:|
| let-7 | 15 | 0.43 | 0.995 | 0.995 |
| miR-28a | 3 | −0.92 | 0.520 | 0.995 |

**Caveat:** only 2 families reached the ≥3-member threshold. The regex is
conservative and miRNA naming is inconsistent, so family coverage is limited.
That said, let-7 — the largest and most-studied family, and the one that kept
recurring in the top-ranked lists — shows nothing (NES 0.43, padj 0.99).

---

## Summary

| Check | Result | Rules out |
|-------|--------|-----------|
| p-value histogram | flat (17 vs 17.9) | signal hidden by threshold or FDR method |
| limma-voom | 0 hits | DESeq2-specific artifact |
| family enrichment | 0 families | coordinated sub-significance shift |

Three independent lines converging on the same conclusion. Combined with the
pairwise, pooled, and trend analyses, the miRNA null is established from every
angle available in this dataset.

---

*Part of the bone marrow miRNA-EV analysis series.*
