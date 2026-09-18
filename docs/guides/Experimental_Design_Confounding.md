# Experimental Design & Confounding Analysis — Bone Marrow miRNA-EV

Documentation of the study design, sample metadata, and the confounding
analysis that determines how the DESeq2 model is built. This bridges
quantification (count matrix) and differential expression (DESeq2).

**Study:** July 2025 Treadmill Experiment (B + C), Bone Marrow sEVs, mouse (mm10)
**Samples:** 62 bone marrow small-EV small RNA libraries

---

## The five variables per sample

| Variable | Type | Levels | Role |
|----------|------|--------|------|
| **Condition** | biological | Sed, Exer, Exer3, Exer7, Exer18 | the effect we test |
| **Sex** | biological | Male, Female | covariate (or interaction, later) |
| **Isolation_Day** | technical | 1–6 (day EVs were isolated) | candidate batch covariate |
| **Seq_Batch** | technical | 1, 2 (two sequencing pools) | candidate batch covariate |
| **Treadmill_Batch** | technical | B, C (animal experiment run in 2 batches) | recorded; **cannot** be modelled (see below) |

> `D#_#` = isolation day _ sample number. Conditions renamed for R safety:
> `Exer+3 → Exer3`, etc. (the `+` causes issues in some R contexts).

> **Seq_Batch boundary:** batch 1 = D1_1 through **D4_1**; batch 2 = **D4_2**
> through D6_12 (each pool also included a water control, not in the count
> matrix). The split falls inside isolation day 4.

---

## Group sizes

| Condition | n |
|-----------|---|
| Sed | 10 |
| Exer | 14 |
| Exer3 | 12 |
| Exer7 | 14 |
| Exer18 | 12 |

**Sex:** 31 Male / 31 Female (perfectly balanced)
**Isolation days:** 10 each on days 1–5, 12 on day 6
**Sequencing batches:** 31 / 31

---

## Confounding analysis

Before modelling: **can the exercise effect be separated from the technical
variables?** The cross-tabs below check this.

### Condition × Sex — balanced ✓

| Condition | Male | Female |
|-----------|------|--------|
| Sed | 5 | 5 |
| Exer | 7 | 7 |
| Exer3 | 6 | 6 |
| Exer7 | 7 | 7 |
| Exer18 | 6 | 6 |

Sex is **not** confounded with condition → clean covariate, and a
sex×condition interaction (sex-specific exercise effects) is feasible.

### Condition × Isolation Day — balanced ✓

| Condition | D1 | D2 | D3 | D4 | D5 | D6 |
|-----------|----|----|----|----|----|----|
| Sed | 2 | 2 | 2 | 2 | 2 | 0 |
| Exer | 2 | 2 | 2 | 2 | 2 | 4 |
| Exer3 | 2 | 2 | 2 | 2 | 2 | 2 |
| Exer7 | 2 | 2 | 2 | 2 | 2 | 4 |
| Exer18 | 2 | 2 | 2 | 2 | 2 | 2 |

Condition is **not** confounded with isolation day → exercise effect can be
separated from isolation-day variation.

### Condition × Seq Batch — spread ✓

| Condition | Batch 1 | Batch 2 |
|-----------|---------|---------|
| Sed | 7 | 3 |
| Exer | 6 | 8 |
| Exer3 | 6 | 6 |
| Exer7 | 6 | 8 |
| Exer18 | 6 | 6 |

Every condition appears in both sequencing batches → not confounded.

### Seq Batch × Isolation Day — (almost) NESTED ⚠️

| | D1 | D2 | D3 | D4 | D5 | D6 |
|--|----|----|----|----|----|----|
| Batch 1 | 10 | 10 | 10 | 1 | 0 | 0 |
| Batch 2 | 0 | 0 | 0 | 9 | 10 | 12 |

Days 1–3 are batch 1, days 5–6 are batch 2, day 4 is split (1 vs 9). They
carry nearly the same information → **use one technical covariate, not both.**
Isolation day (finer, 6 levels) is preferred and tracks the batch split.

### Condition × Treadmill Batch — CONFOUNDED ✗ (cannot model)

| Condition | B | C |
|-----------|---|---|
| Sed | 10 | 0 |
| Exer | 10 | 4 |
| Exer3 | 0 | 12 |
| Exer7 | 10 | 4 |
| Exer18 | 0 | 12 |

**Sed is 100% batch B; Exer3 and Exer18 are 100% batch C.** Treadmill batch is
partly collinear with condition — for these groups, batch is perfectly
determined by condition. Therefore **Treadmill_Batch cannot be a model
covariate** (DESeq2 cannot separate the Exer3 effect from the batch-C effect).

It is still **recorded** in the metadata and used for **PCA colouring** — to
visually check whether B/C drives clustering — but never enters the design
formula.

> **Interpretive caveat for the write-up:** because Sed/Exer3/Exer18 are each
> entirely one treadmill batch, any treadmill-batch technical effect for those
> conditions is inseparable from the biological effect. Worth one sentence in
> the methods/limitations. Common in staged animal experiments; does not
> invalidate the analysis.

---

## The model

The design is well-balanced: condition is not confounded with sex, isolation
day, or sequencing batch. Constraints: (1) batch and isolation day are nested —
pick one; (2) treadmill batch is confounded with condition — exclude from model.

**Main model:**
```
~ Isolation_Day + Sex + Condition
```
- `Isolation_Day` (6 levels) — technical batch covariate (subsumes seq batch)
- `Sex` — balanced covariate; exercise effects averaged across sexes
- `Condition` — tested variable, **reference = Sed**

**Alternative** if isolation day is over-parameterized (uses 5 df): substitute
`Seq_Batch` (1 df). Decide from the PCA (does clustering follow isolation day
or just the two batches?).

**Later — sex-specific effects** (sex is balanced, so supported):
```
~ Isolation_Day + Sex + Condition + Sex:Condition
```
The `Sex:Condition` term tests whether the exercise response differs by sex.
Alternatively, run the pipeline separately per sex.

**Reference level & "comparing all":** DESeq2 fits one model; contrasts are
pulled from it. With Sed as reference the default contrasts are
Exer/Exer3/Exer7/Exer18 vs Sed. Any other pairwise comparison (e.g. Exer7 vs
Exer) is still available via
`results(dds, contrast = c("Condition","Exer7","Exer"))`.

---

## Filtering (for the DESeq2 step)

1. **Drop failed samples** — total counts below a tiny threshold (safety net;
   none dropped here).
2. **Group-aware gene filter** — keep a miRNA only if ≥10 counts in ≥3 samples
   of at least one condition. Removes sparse noise, protects condition-specific
   miRNAs. Smallest group (Sed, n=10) > 3, so the threshold is safe.

---

## Full sample metadata

`metadata_BoneMarrow.txt` — `Index_Number` matches the count-matrix columns.

| Sample | Condition | Sex | Isolation_Day | Seq_Batch | Treadmill_Batch |
|--------|-----------|-----|---------------|-----------|-----------------|
| D1_1 | Sed | Male | 1 | 1 | B |
| D1_2 | Sed | Female | 1 | 1 | B |
| D1_3 | Exer | Male | 1 | 1 | B |
| D1_4 | Exer | Female | 1 | 1 | C |
| D1_5 | Exer3 | Male | 1 | 1 | C |
| D1_6 | Exer3 | Female | 1 | 1 | C |
| D1_7 | Exer7 | Male | 1 | 1 | B |
| D1_8 | Exer7 | Female | 1 | 1 | C |
| D1_9 | Exer18 | Male | 1 | 1 | C |
| D1_10 | Exer18 | Female | 1 | 1 | C |
| D2_1 | Sed | Male | 2 | 1 | B |
| D2_2 | Sed | Female | 2 | 1 | B |
| D2_3 | Exer | Male | 2 | 1 | C |
| D2_4 | Exer | Female | 2 | 1 | B |
| D2_5 | Exer3 | Male | 2 | 1 | C |
| D2_6 | Exer3 | Female | 2 | 1 | C |
| D2_7 | Exer7 | Male | 2 | 1 | C |
| D2_8 | Exer7 | Female | 2 | 1 | B |
| D2_9 | Exer18 | Male | 2 | 1 | C |
| D2_10 | Exer18 | Female | 2 | 1 | C |
| D3_1 | Sed | Male | 3 | 1 | B |
| D3_2 | Sed | Female | 3 | 1 | B |
| D3_3 | Exer | Male | 3 | 1 | B |
| D3_4 | Exer | Female | 3 | 1 | C |
| D3_5 | Exer3 | Male | 3 | 1 | C |
| D3_6 | Exer3 | Female | 3 | 1 | C |
| D3_7 | Exer7 | Male | 3 | 1 | B |
| D3_8 | Exer7 | Female | 3 | 1 | C |
| D3_9 | Exer18 | Male | 3 | 1 | C |
| D3_10 | Exer18 | Female | 3 | 1 | C |
| D4_1 | Sed | Male | 4 | 1 | B |
| D4_2 | Sed | Female | 4 | 2 | B |
| D4_3 | Exer | Male | 4 | 2 | C |
| D4_4 | Exer | Female | 4 | 2 | B |
| D4_5 | Exer3 | Male | 4 | 2 | C |
| D4_6 | Exer3 | Female | 4 | 2 | C |
| D4_7 | Exer7 | Male | 4 | 2 | C |
| D4_8 | Exer7 | Female | 4 | 2 | B |
| D4_9 | Exer18 | Male | 4 | 2 | C |
| D4_10 | Exer18 | Female | 4 | 2 | C |
| D5_1 | Sed | Male | 5 | 2 | B |
| D5_2 | Sed | Female | 5 | 2 | B |
| D5_3 | Exer | Male | 5 | 2 | B |
| D5_4 | Exer | Female | 5 | 2 | B |
| D5_5 | Exer3 | Male | 5 | 2 | C |
| D5_6 | Exer3 | Female | 5 | 2 | C |
| D5_7 | Exer7 | Male | 5 | 2 | B |
| D5_8 | Exer7 | Female | 5 | 2 | B |
| D5_9 | Exer18 | Male | 5 | 2 | C |
| D5_10 | Exer18 | Female | 5 | 2 | C |
| D6_1 | Exer | Male | 6 | 2 | B |
| D6_2 | Exer | Female | 6 | 2 | B |
| D6_3 | Exer | Male | 6 | 2 | B |
| D6_4 | Exer | Female | 6 | 2 | B |
| D6_5 | Exer3 | Male | 6 | 2 | C |
| D6_6 | Exer3 | Female | 6 | 2 | C |
| D6_7 | Exer7 | Male | 6 | 2 | B |
| D6_8 | Exer7 | Female | 6 | 2 | B |
| D6_9 | Exer7 | Male | 6 | 2 | B |
| D6_10 | Exer7 | Female | 6 | 2 | B |
| D6_11 | Exer18 | Male | 6 | 2 | C |
| D6_12 | Exer18 | Female | 6 | 2 | C |

---

## Summary of design decisions

| Variable | Decision |
|----------|----------|
| Condition | tested; reference = Sed |
| Sex | covariate (balanced); interaction available later |
| Isolation_Day | technical covariate in main model |
| Seq_Batch | recorded; nested with isolation day → not used alongside it |
| Treadmill_Batch (B/C) | recorded; confounded with condition → PCA only, never in model |
| Model | `~ Isolation_Day + Sex + Condition` |

---

*Reference for the bone marrow miRNA-EV project. Balanced for condition vs sex
and vs isolation day; sequencing batch nested within isolation day; treadmill
batch confounded with condition (visualize only).*
