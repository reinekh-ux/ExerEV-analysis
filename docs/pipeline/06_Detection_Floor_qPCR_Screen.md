# Detection Floor & qPCR Candidate Screen

Establishes **what effect size this study could detect**, then uses that to
decide which miRNAs are worth independent validation.

**Script:** `qpcr_candidate_screen.R`
**Output:** `DESeq2_results/qpcr_candidates/`

```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
Rscript qpcr_candidate_screen.R
```

---

## Why this rather than another ranking

Any test run on these 298 miRNAs to rank candidates draws from the same flat
p-value distribution. A new test gives a new ranking, equally noisy, with the
false comfort of a fresh number attached — the trap let-7e fell into, since it
*was* top-ranked by a legitimate test.

Ranking by fold change fails differently: the largest log2FCs in this dataset
belong to near-zero miRNAs (miR-1843a at −24.8, miR-134-5p at −19.4). Fold change
ignores measurement reliability.

The detection floor sidesteps both. It is a property of the **design**, not of
the results, so it cannot be contaminated by the rankings it is used to filter.

---

## The three filters

Starting from 298 miRNAs:

**1. Abundance** — baseMean ≥ 50. qPCR needs enough material to measure
reliably. Removes the low-count artifacts by construction.

**2. Direction consistency** — all four exercise group medians on the same side
of Sed (`up` or `down`). Mixed patterns dropped — this is why let-7e does not
appear despite topping the LRT.

*Important limitation:* four independent groups landing the same way occurs by
chance 2 × (1/2)⁴ = **1 in 8**. On 298 miRNAs that is ~37 by luck. This filter
narrows; it does not select.

**3. Detection floor** — minimum detectable effect at 80% power, computed per
miRNA from its actual DESeq2 dispersion and the group sizes (10 Sed vs 52
exercise).

```r
se_log2 <- sqrt((1/(n_sed*bm) + 1/(n_exr*bm) + disp*(1/n_sed + 1/n_exr))) / log(2)
mde <- (qnorm(1-alpha/2) + qnorm(power)) * se_log2
mde_pct <- (2^mde - 1) * 100
```

*Caveat:* normal approximation to negative-binomial power. An order-of-magnitude
estimate, not a number to quote precisely.

---

## Results — the floor is high

Across miRNAs with baseMean ≥ 50, the minimum detectable effect is typically
**20–50%**, and worse for mid-abundance features:

| miRNA | baseMean | MDE |
|-------|---:|---:|
| miR-9-5p | 118 | 150% |
| miR-455-3p | 58 | 83% |
| miR-101c | 50 | 73% |
| miR-133a-3p | 151 | 74% |
| miR-744-5p | 89 | 45% |
| miR-150-5p | 300 | 41% |
| miR-21a-5p | 7125 | **14.6%** |
| let-7i-5p | 3346 | 14.3% |
| let-7b-5p | 3106 | 14.7% |
| miR-29a-3p | 2162 | 15.2% |

**44 of 46 shortlisted miRNAs have `below_floor = TRUE`** — the observed shift is
smaller than what the design could detect. That is a finding about the study, not
about biology.

`padj_Exer3` across the entire shortlist runs 0.87–1.0.

---

## What this reframes

The miRNA result becomes a **precise** claim rather than a vague one:

> Not "we found nothing" — but "we were powered to detect effects above
> approximately 15% in the most abundant miRNAs and 40–80% in mid-abundance
> ones, and found none above that threshold."

The first invites *did you look hard enough*. The second answers it.

It also explains the shape of every candidate list produced in this project: the
largest apparent changes consistently land on the **least abundant** miRNAs,
because that is where the noise is largest. The same low-count effect that
generated the padj artifacts, operating at baseMean 50–150 instead of baseMean 1.

---

## What it says about miR-21a

miR-21's MDE is **14.6%** — one of the lowest in the dataset, because it is by far
the most abundant miRNA (baseMean 7125). The observed peak shift is ~24%
(day-adjusted, Sed→Exer3).

So miR-21 is **not** below its detection floor. The screen had the power to catch
a 24% effect in this miRNA and did not return it as significant.

This weakens rather than strengthens the case for miR-21, and it corrects an
earlier working assumption that the floor would sit around 40%. Where power
existed, nothing moved; where things appear to move, power was absent.

---

## Known bug

`peak_change_pct` returned `NA` for miR-21a-5p and miR-25-3p — a division guard
issue on the most abundant rows. Fixed in `simple_candidates.R` with a
`sed > 0` condition. The MDE values are unaffected.

---

## How to read `below_floor`

`TRUE` means *"this study could not have seen it either way"* — not *"it is
probably real."* Most of the shortlist is noise pointing consistently by chance
(~37 expected). Being on the list is weak evidence.

What separated miR-21 was an **independent prior** — the exercise/bone-marrow
literature would have named it before the data were examined. A data-driven filter
combined with an independent prior holds up; either alone does not.

---

*Part of the bone marrow miRNA-EV analysis series.*
