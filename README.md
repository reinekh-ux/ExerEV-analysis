# ExerEV — Exercise-Induced Small Extracellular Vesicle miRNA Analysis

Bioinformatics analysis code and methods documentation for the ExerEV project (Nagy Lab, McGill University / Douglas Research Centre), investigating exercise-induced small extracellular vesicles (sEVs) as candidate mediators of exercise's antidepressant effects.

Small RNA-seq data (mm10, exceRpt pipeline) from a C57BL/6 treadmill mouse model was used to profile sEV miRNA cargo across four post-exercise timepoints (0h/3h/7h/18h vs. sedentary) in two tissues: **skeletal muscle** and **bone marrow**.

> This repo contains analysis code and methods write-ups only — no raw sequencing data, count matrices, or figures are included (see `.gitignore`).

## Repo structure

```
docs/
  pipeline/    Numbered, ordered walkthrough of the analysis (01–10)
  guides/      Standalone methods guides & notebooks (QC, DESeq2, quantification, etc.)
scripts/
  bone_marrow/   DE analysis, QC, PCA, targeted/pooled/all-contrasts models, figures
  muscle/        Candidate screening and figure generation (versioned iterations)
  cross_tissue/  Cross-tissue comparison figure (muscle vs. bone marrow log2FC)
  pipeline/      exceRpt quantification cluster job script
```

## Suggested reading order

1. `docs/guides/exceRpt_Quantification_Guide.md` — raw read → miRNA quantification (exceRpt, mm10)
2. `docs/guides/Count_Matrix_Assembly_Guide.md` — assembling count matrices from exceRpt output
3. `docs/guides/Exploratory_Variance_Analysis_Notebook.md` + `Experimental_Design_Confounding.md` — QC and design considerations
4. `docs/guides/DESeq2_miRNA_Analysis_Notebook.md` / `DE_Analysis_Code_Notebook.md` — core DE modeling
5. `docs/pipeline/01` → `10` — ordered walkthrough from trend analysis through the final muscle scripts and pipeline explanation
6. `docs/guides/BoneMarrow_miRNA_EV_QC_Workflow.md` and `docs/guides/mir21_mechanistic_notebook.md` — tissue-specific / mechanistic deep-dives

## Pipeline & tools

- **Quantification:** exceRpt (STAR alignment, `outFilterMismatchNmax=0`, EndToEnd, mm10)
- **Differential expression:** DESeq2 (negative binomial), `~ Isolation_Day + Sex + Condition`
- **Clustering:** k-means (Euclidean, Hartigan-Wong, nstart=50, iter.max=100, seed=910), k=5, amplitude-normalized by max |log2FC|
- **Enrichment:** multiMiR → clusterProfiler (BP/CC/MF)
- **Environment:** R 4.4.2 / Bioconductor 3.20 (local); Alliance Canada Rorqual cluster (R 4.2.1 / Bioconductor 3.16) for HPC jobs

## Notes on the scripts

- Several scripts (`bm_*`, `muscle_figures_*`) exist in multiple numbered versions reflecting iterative development; the highest version number in each family is the most current.
- Cluster scripts contain hard-coded Alliance Canada (Rorqual) filesystem paths under `/lustre09/project/6019280/reinekh/...` — update these before running elsewhere.
- Scripts assume the R/Bioconductor package set noted above (DESeq2, clusterProfiler, org.Mm.eg.db, enrichplot, multiMiR, limma, edgeR, fgsea, ggplot2).
