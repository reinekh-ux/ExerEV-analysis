# Candidate Boxplots — miR-21a-5p and let-7e-5p

Per-mouse expression plots for the two miRNAs that recurred across analyses.
Shows what the data actually look like, rather than arguing about p-values.

**Scripts:** `candidate_boxplots.R`, `candidate_boxplots_day.R`
**Output:** `DESeq2_results/candidate_boxplots/`

```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
Rscript candidate_boxplots_day.R
```

---

## Why boxplots

A p-value argument is abstract. Sixty-two dots showing near-total overlap is not.
These plots let a reader judge effect size directly — and they are honest in both
directions: if something were visible, it would show.

Each point is one mouse; values are DESeq2 size-factor normalised counts (the
same normalisation the DE tests used). Both linear and log10 versions are
produced — miR-21 is abundant with a wide spread, so log scale often reads better.

---

## Condition × Isolation_Day cross-tab

Checked first, because it determines whether day-adjustment is meaningful:

```
        1 2 3 4 5 6
Sed     2 2 2 2 2 0
Exer    2 2 2 2 2 4
Exer3   2 2 2 2 2 2
Exer7   2 2 2 2 2 4
Exer18  2 2 2 2 2 2
```

Close to 2 per cell throughout — condition and isolation day are **not
confounded**. Good design, and it means the day-adjustment is doing real work
rather than fighting an entangled structure.

---

## miR-21a-5p

**Raw medians by condition:**

| Sed | Exer | Exer3 | Exer7 | Exer18 |
|---:|---:|---:|---:|---:|
| 6375 | 6548 | 7681 | 7525 | 7173 |

**Day-adjusted log2 medians:**

| Sed | Exer | Exer3 | Exer7 | Exer18 |
|---:|---:|---:|---:|---:|
| 12.60 | 12.70 | 12.91 | 12.86 | 12.81 |

**Medians by isolation day:**

| 1 | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|
| 6450 | 7023 | 7159 | 7728 | 6892 | 7130 |

**Reading it.** All four exercise groups sit above Sed, peaking around +18–20% at
Exer3/Exer7. Direction survives day adjustment: Sed→Exer3 is 0.31 log2, ~24%.

But the day-median range (6450–7728, span ~1280) is comparable to the
condition-median range (6375–7681, span ~1300). **Isolation day moves this miRNA
as much as condition does.** Day 4 alone exceeds every condition median.

Group SDs are 677–1097 against a ~1150 median gap, so distributions overlap
heavily — Sed's upper points sit above several Exer3 points. That overlap is what
padj 0.80 encodes.

**Statistical context:** pooled Exercise vs Sed gave log2FC +0.20, raw p 0.008,
padj 0.80. The pooled estimate is diluted because Exer (6548) barely moves while
Exer3/Exer7 do — averaging all four drags it down.

---

## let-7e-5p

**Raw medians by condition:**

| Sed | Exer | Exer3 | Exer7 | Exer18 |
|---:|---:|---:|---:|---:|
| 176.3 | 118.7 | 198.0 | 149.4 | 160.8 |

**Day-adjusted log2 medians:**

| Sed | Exer | Exer3 | Exer7 | Exer18 |
|---:|---:|---:|---:|---:|
| 7.43 | 6.90 | 7.75 | 7.14 | 7.28 |

**Medians by isolation day:**

| 1 | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|
| 144.7 | 146.3 | 119.1 | 163.8 | 188.8 | 181.4 |

**Reading it.** A zigzag, not a trend: down, up past baseline, down, up. Sed is
neither highest nor lowest. Day adjustment does not resolve it — the adjusted
values zigzag identically.

Day-median range (119–189, span 70) is essentially equal to the condition-median
range (119–198, span 79).

Pooled Sed vs exercise: 176.3 vs 160.3. Essentially nothing.

**This is why let-7e ranked top under LRT** — the LRT rewards *any* pattern of
group differences, including non-monotonic ones. Its pooled log2FC was −0.03
while its LRT padj was 0.0015. The groups differ from one another without
differing in any direction.

---

## On "correcting for" isolation day

Two distinct operations, only one of which is a diagnostic:

**Showing** day structure (colour or facet by day) reveals whether a shift holds
within days or tracks day. This is informative.

**Removing** it (`limma::removeBatchEffect`, protecting Condition + Sex) produces
a cleaner plot, but plots model output rather than data — and the model already
reported padj 0.80. A residual plot cannot overturn the test that generated it.

```r
design_keep <- model.matrix(~ Condition + Sex, data=metadata)
logn_adj <- removeBatchEffect(logn, batch=metadata$Isolation_Day, design=design_keep)
```

Both are produced; the day-adjusted version is labelled `MODEL OUTPUT, not raw
data` in its subtitle.

---

## Conclusion

**let-7e:** drop it. The zigzag is not an exercise response, day explains as much
as condition, and the direction contradicts the literature prior (HERITAGE
predicted down with training).

**miR-21a:** direction is consistent across four groups and survives day
adjustment at ~24%. That is a real pattern in the medians. Whether it is a
detectable effect is a separate question — see the detection floor notebook,
where miR-21's MDE is ~15%, meaning the screen *should* have caught a 24% effect
and did not.

---

*Part of the bone marrow miRNA-EV analysis series.*
