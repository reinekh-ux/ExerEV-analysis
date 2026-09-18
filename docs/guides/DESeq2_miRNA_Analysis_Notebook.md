# DESeq2 Differential Expression — miRNA (Bone Marrow EV) — Analysis Notebook

Complete record of the differential expression analysis testing whether exercise
alters bone marrow small-EV miRNA cargo. Every model, threshold, and framing was
tested; the result is a robust null. Runs in R on Rorqual via `Rscript`.

**Load environment:**
```bash
module load StdEnv/2020 gcc/9.3.0 r/4.2.1 r-bundle-bioconductor/3.16
```

**Inputs:** `count_matrices/BoneMarrow_miRNA_TotalRawCount_matrix.txt` (697 miRNAs × 62),
`metadata_BoneMarrow.txt` (62 samples).
**Conditions:** Sed, Exer (0h), Exer3 (3h), Exer7 (7h), Exer18 (18h) — hours post
final bout after 3 weeks treadmill training. Reference = Sed.
**Model (evidence-based):** `~ Isolation_Day + Sex + Condition`.

---

## Pipeline (shared across all analyses)

```
Raw counts (697 miRNAs)
  → drop failed samples (total < 3; none dropped)
  → group-aware gene filter (>=10 counts in >=3 samples of any group) → 298 kept
  → DESeq2 (median-of-ratios normalization, internal)
  → fit model, extract contrasts
  → judge hits by: padj < 0.05 AND baseMean >= 10 (removes low-count artifacts)
```

**Key interpretation rule — baseMean:** a miRNA's average normalized abundance.
Low-baseMean miRNAs (<~5) produce huge fake fold changes from a few stray reads
(e.g. log2FC −24 at baseMean 1). These are **low-count artifacts** — statistically
"significant" but biologically meaningless. Every result is judged on *real* hits
(baseMean ≥ 10), not the raw significant count.

---

## Analyses run

### 1. Per-timepoint (each exercise timepoint vs Sed)
Model `~ Isolation_Day + Sex + Condition`. Four contrasts.
**Result:** 0 real hits. The 1–2 "significant" hits were baseMean <2 artifacts.

### 2. Covariate robustness (3 models)
Same contrasts under three technical-covariate choices:

| Model | Real hits |
|-------|:---:|
| `~ Isolation_Day + Sex + Condition` | 0 |
| `~ Seq_Batch + Sex + Condition` | 0 |
| `~ Condition` (muscle-matched) | 0 |

Result is independent of covariate choice.

### 3. Threshold sensitivity (4 stringency levels)
One run, four filter/significance settings, organized into subfolders:

| Setting | filter | padj | miRNAs tested | sig_total | **sig_real** | sig_real_lfc |
|---------|--------|------|:---:|:---:|:---:|:---:|
| strict | 10 in 5 | 0.05 | 233 | 0 | 0 | 0 |
| standard | 10 in 3 | 0.05 | 298 | 2 | 0 | 0 |
| lenient | 5 in 3 | 0.10 | 354 | 21 | 1 | 0 |
| very_lenient | 5 in 2 | 0.10 | 423 | 63 | 4 | 1 |

`sig_total` grows only because looser filters admit more low-count artifacts;
`sig_real` (baseMean-filtered) stays ~0. Hits appear only at maximum permissiveness
= noise crossing a relaxed threshold, not signal.

### 4. Pooled Exercise vs Sedentary (maximum power)
Collapsed all exercise into one group: **52 Exercise vs 10 Sed**, `~ Isolation_Day
+ Sex + ExerciseYesNo`.
**Result:** 0 real hits. Most abundant miRNA (miR-21a-5p, baseMean 7124) changed by
only log2FC +0.20 (~15%). The best-measured miRNA barely moves.

### 5. Ordered trend (Sed → 0 → 3 → 7 → 18h)
Polynomial trend contrasts (built via `poly()` on ordered positions, since DESeq2
rejects ordered factors in the formula), `~ Isolation_Day + trend_linear +
trend_quadratic`.
- **Linear trend** (steady progression): 0 real hits
- **Quadratic trend** (hump/U-shape): 0 real hits (1 total = a baseMean 1.8 artifact)

No dose-response, no non-linear recovery dynamics.

---

## Summary of every test

| Analysis | Hypothesis tested | Real hits |
|----------|-------------------|:---:|
| Per-timepoint (×4) | each timepoint differs from Sed | 0 |
| 3 covariate models | robustness to technical covariate | 0 |
| Sensitivity (×4 thresholds) | robustness to stringency | 0 |
| Pooled Exercise vs Sed | any exercise effect (max power) | 0 |
| Linear ordered trend | steady progression over recovery | 0 |
| Quadratic ordered trend | hump / U-shape over recovery | 0 |

**Conclusion:** No differential expression of bone marrow EV miRNA cargo with
exercise — no pairwise difference, no dose-response, no linear or non-linear trend
— robust across covariate specifications and stringency thresholds.

---

## The consistent (non-significant) whisper

Several abundant miRNAs recur near the top across analyses, mostly trending down,
never surviving correction: **miR-21a-5p, let-7 family (let-7e/f/g), miR-148a,
miR-223, miR-101**. Worth one honest sentence in the write-up ("a consistent but
non-significant downward trend in several abundant miRNAs"), not more. miR-21 is
notable — prior literature links it to exercise and bone marrow osteogenesis — but
here it does not reach significance.

---

## Technical notes / gotchas encountered

| Issue | Resolution |
|-------|-----------|
| EnhancedVolcano not in cluster Bioconductor bundle | swapped for base ggplot2 volcano |
| Volcano showed more "significant" dots than sig_real | volcano lacked the baseMean filter; added baseMean-aware colouring + point sizing + x-axis cap |
| DESeq2 rejects ordered factors in design formula | built polynomial trend columns via `poly()`, passed as numeric predictors |
| Low-count artifacts inflating hit counts | judge hits on baseMean ≥ 10, not raw padj alone |
| nano paste saved empty files / heredoc mangled over SSH | use `cat > file << 'EOF'` carefully; verify with `tail -3` |

---

## Output organization

```
DESeq2_results/
├── SENSITIVITY_SUMMARY.txt
├── strict/  standard/  lenient/  very_lenient/     (sensitivity)
├── ExerciseVsSed/                                   (pooled)
├── OrderedTrend/                                    (linear + quadratic)
└── each subfolder: results tables, volcano, normalized counts
```

---

## Biological context (for discussion)

Prior work shows exercise *training* can remodel bone marrow EV miRNA cargo, but the
clearest effects appear in **metabolically stressed** models (high-fat diet,
hypertension, irradiation) where exercise counteracts a perturbation. In healthy,
lean, chronically-trained mice, this well-powered analysis (n=62, far larger than
typical n=3–6 EV cargo screens) finds bone marrow EV miRNA cargo is stable. This
defines a boundary condition: exercise does not substantially alter the miRNA
compartment of bone marrow EVs in the absence of metabolic stress.

---

## Next steps (miRNA settled)

- **tRNA & piRNA** — does any small RNA respond where miRNA doesn't? (matrices ready)
- **Method cross-check** (limma-voom / edgeR) — confirm the null is method-independent
- **Muscle comparison** — cross-tissue contrast
- **Literature cross-check** — are prior studies' reported miRNAs among our candidates?

---

*Reference for the bone marrow miRNA-EV project. miRNA differential expression
tested exhaustively (6 analytical framings); robust null. Judged on baseMean-filtered
hits to exclude low-count artifacts.*
