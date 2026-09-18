# Likelihood Ratio Test (LRT) — Condition across all five levels

An **omnibus** test: does Condition explain a significant amount of variation
across all five groups jointly, rather than in any one pairwise comparison?

**Script:** `DEG_miRNA_LRT.R`
**Output:** `DESeq2_results/LRT/`
**Reference:** HBC training, [DGE LRT lesson](https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html)

```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
Rscript DEG_miRNA_LRT.R
```

---

## How it works

Compares two nested models and asks whether dropping Condition significantly
worsens the fit:

```r
dds <- DESeqDataSetFromMatrix(counts, metadata,
                              design = ~ Isolation_Day + Sex + Condition)
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ Isolation_Day + Sex)
res <- results(dds_lrt)
```

| Model | Formula |
|-------|---------|
| Full | `~ Isolation_Day + Sex + Condition` |
| Reduced | `~ Isolation_Day + Sex` |

The difference in deviance between them becomes the test statistic. One p-value
per miRNA, covering all five levels at once.

**Note on the reduced model:** HBC's example uses `reduced = ~ 1` because their
full model contains only one factor. Here the covariates are retained in the
reduced model so the test isolates **Condition specifically**, rather than
asking "does anything at all matter."

---

## Two interpretation rules

**1. log2FoldChange is not the tested quantity.** DESeq2 still reports a fold
change under the LRT, but it defaults to the last contrast and has no
relationship to the hypothesis being tested. Judge hits on **padj and baseMean
only**. Never apply a fold-change filter to LRT results.

**2. HBC recommend a stricter FDR threshold.** Their reasoning, quoted in
substance: the LRT typically returns many more genes than pairwise tests, and
because you cannot filter by fold change, you lose your usual second layer of
stringency. They suggest padj < 0.001 as compensation.

That advice addresses a flood of hits. This dataset produced two, so the premise
does not apply here — the 0.001 line was reported for completeness, not because
it is the operative criterion.

---

## Results

Filtered: 298 miRNAs × 62 samples.
Group sizes: Sed=10, Exer=14, Exer3=12, Exer7=14, Exer18=12.

```
raw p < 0.05:  17 observed  vs  14.9 expected by chance
padj < 0.05:    2 total  |  1 real (baseMean >= 10)
padj < 0.001:   0 total  |  0 real
```

**Top hits:**

| miRNA | baseMean | log2FC | pvalue | padj |
|-------|---:|---:|---:|---:|
| let-7e-5p | 166.8 | −0.03 | 5.0e-06 | **0.0015** |
| miR-7b-5p | 9.4 | +5.81 | 7.6e-05 | 0.011 |
| miR-223-3p | 1354.5 | −0.17 | 0.0022 | 0.221 |
| miR-19b-3p | 70.6 | −0.25 | 0.0034 | 0.253 |
| miR-21a-5p | 7124.5 | +0.27 | 0.0062 | 0.308 |
| miR-125b-5p | 1847.1 | +0.12 | 0.0053 | 0.308 |

---

## Interpretation — why let-7e is not a finding

let-7e-5p is the only hit in the entire miRNA analysis that is both significant
at padj < 0.05 and well-expressed (baseMean 167, not a low-count artifact). It
warrants a careful look, and it does not survive one.

**The p-value distribution is flat.** 17 observed below 0.05 against 14.9
expected — an excess of ~2 out of 298 tests. Genuine structure that pairwise
tests missed would produce a substantial pile-up of small p-values, not two.

**It fails the recommended threshold.** padj 0.0015 > 0.001.

**Its log2FC is −0.03.** Essentially zero. The LRT ranked it on the *pattern
across groups*, not on any coherent shift.

**Its group medians zigzag:** Sed 176 → Exer 119 → Exer3 198 → Exer7 149 →
Exer18 161. Down, up past baseline, down, up. Not a trajectory — the groups
differ from one another without differing in any direction.

**Isolation day explains a comparable range.** Day medians span 119–189; condition
medians span 119–198. Nearly identical.

**Literature caution.** let-7e is a documented exercise-responsive miRNA (HERITAGE
Family Study: decreased in human serum after 20 weeks of endurance training,
FDR q < 0.05), and let-7 family members are established regulators of osteogenesis
in bone marrow MSCs. This is a coherent story — which is precisely the danger.
With 298 miRNAs tested, whichever ranked top would have *some* literature; let-7
is a large, heavily-studied family. Post-hoc plausibility is nearly free and is
not evidence. Note also that the HERITAGE direction (down with training) does not
match these data, where Sed sits mid-range.

**Write-up sentence:**

> let-7e-5p ranked top under the likelihood ratio test (padj = 0.0015) and is a
> reported exercise-responsive miRNA, but did not survive the stricter threshold
> recommended for LRT, and the uniform p-value distribution (17 observed vs 14.9
> expected below 0.05) indicates no genuine differential expression.

miR-7b-5p (baseMean 9.4, log2FC +5.8) is a borderline low-count artifact.

---

*Part of the bone marrow miRNA-EV analysis series.*
