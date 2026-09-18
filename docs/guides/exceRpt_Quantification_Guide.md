# miRNA Quantification with exceRpt — Step-by-Step Guide

Quantifying trimmed bone marrow miRNA-EV libraries on Rorqual using the
**exceRpt** small RNA-seq pipeline (containerized with Apptainer), then
merging the per-sample outputs into one count matrix. Written to match the
QC/trimming guides — beginner-friendly, with explanations.

**Input:** 62 trimmed FASTQ files (`*_bbduk_min16.fastq.gz`)
**Output:** per-sample small RNA counts → one combined miRNA × 62 matrix
**Reference:** mouse genome **mm10**
**Tool:** exceRpt v4.6.3, run from a pre-built container (`exceRpt.sif`)

---

## What exceRpt does

exceRpt (extra-cellular RNA processing toolkit) is purpose-built for
**EV / exRNA small RNA-seq**. For each trimmed read it:

1. Filters out low-quality, homopolymer, contaminant (UniVec), and rRNA reads
2. Aligns the rest to the mm10 genome and small RNA databases (miRBase for
   miRNAs, plus tRNA, piRNA, gencode, circRNA)
3. Counts how many reads map to each small RNA → per-sample count tables

It runs from an **Apptainer container** — a single sealed `.sif` file
bundling exceRpt and all its dependencies. No installation needed (this is
why it avoids the environment problems that plagued the MultiQC setup).

---

## The big-picture flow

```
62 trimmed FASTQ
       │
       ▼
[launcher script]  ──► submits 62 separate SLURM jobs (one per sample)
       │
       ▼
[exceRpt container, per sample]  ──► filter → align to mm10 → count small RNAs
       │
       ▼
62 output folders (one per sample, each with a .stats file + count tables)
       │
       ▼
[merge job: mergePipelineRuns.R]  ──► one combined miRNA × 62-sample matrix
       │
       ▼
   DESeq2 (sedentary vs exercise)
```

---

## Key paths

| What | Path |
|------|------|
| Trimmed input | `~/links/scratch/BoneMarrow_miRNA_EV_QC/trimmed` |
| Per-sample output | `~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_output_mm10` |
| Merged output | `~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_merged` |
| Container + merge scripts | `~/links/projects/rrg-gturecki/shared/Software_and_Installation/exceRpt` |
| mm10 reference DB | `~/links/projects/rrg-gturecki/shared/Software_and_Installation/exceRpt_mm10/mm10` |
| Compute account | `def-cnagy` |

---

## The scripts (how the work is structured)

Unlike the trimming job (one job looping over all samples), exceRpt uses
**separate scripts**, because each sample is heavy (35 GB RAM, up to a day)
and running them as separate jobs lets them run **in parallel**:

1. **`exceRpt_job_mm10.sh`** — the SLURM job (`#SBATCH` header). Processes
   ONE sample. Takes two arguments: the FASTQ filename and the sample name.
2. **`launch_exceRpt_all.sh`** — a plain loop (no `#SBATCH`) that calls
   `sbatch` 62 times, once per sample. Run with `bash`, not `sbatch`.
3. **`merge_exceRpt.sh`** — a SLURM job that runs the merge script inside the
   container after all 62 are done.

---

# PART 1 — One-time check

Confirm the shared container and reference exist and are readable:

```bash
ls ~/links/projects/rrg-gturecki/shared/Software_and_Installation/exceRpt/
ls ~/links/projects/rrg-gturecki/shared/Software_and_Installation/exceRpt_mm10/
```

Look for `exceRpt.sif` and the merge scripts in the first, and an `mm10`
folder in the second.

---

# PART 2 — Create the per-sample job script

```bash
cat > ~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_job_mm10.sh << 'EOF'
#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=35G
#SBATCH -o exceRpt_%x_%j.log
#SBATCH --mail-user=reine.khoury@mail.mcgill.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020
module load apptainer/1.1.8

echo "Input file: $1"
echo "Sample name: $2"

SHARED=~/links/projects/rrg-gturecki/shared/Software_and_Installation
INPUT_DIR=~/links/scratch/BoneMarrow_miRNA_EV_QC/trimmed
OUTPUT_DIR=~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_output_mm10

apptainer run \
  -B $INPUT_DIR:/exceRptInput \
  -B $OUTPUT_DIR/$2:/exceRptOutput \
  -B $SHARED/exceRpt_mm10/mm10:/exceRpt_DB/mm10 \
  $SHARED/exceRpt/exceRpt.sif \
  INPUT_FILE_PATH=/exceRptInput/$1 \
  SAMPLE_NAME=$2 \
  ADAPTER_SEQ=none \
  N_THREADS=8 \
  REMOVE_LARGE_INTERMEDIATE_FILES=TRUE \
  STAR_outFilterMatchNmin=16 \
  STAR_outFilterMismatchNmax=0 \
  STAR_alignEndsType=EndToEnd \
  MAIN_ORGANISM_GENOME_ID=mm10 \
  MIN_READ_LENGTH=16
EOF
```

### What each setting means

| Setting | Meaning |
|---------|---------|
| `--time=1-00:00:00` | up to 1 day per job (exceRpt is slow) |
| `--mem=35G` | exceRpt holds the genome index in memory |
| `apptainer run` | runs the sealed container |
| `-B X:Y` | "bind" — makes folder X on the cluster visible inside the container as Y |
| `ADAPTER_SEQ=none` | DON'T re-trim — reads were already trimmed with BBDuk |
| `STAR_outFilterMismatchNmax=0` | zero mismatches (matches muscle dataset) |
| `STAR_outFilterMatchNmin=16` | require ≥16 bp aligned |
| `STAR_alignEndsType=EndToEnd` | align whole read, no soft-clipping |
| `MAIN_ORGANISM_GENOME_ID=mm10` | align to mouse mm10 |
| `MIN_READ_LENGTH=16` | ignore reads under 16 nt |

> All alignment parameters match the muscle dataset exactly, so the two
> tissues are directly comparable.

---

# PART 3 — Test ONE sample first (always)

exceRpt is the heaviest step — test one before launching all 62.

```bash
# make the output folder for the test sample
mkdir -p ~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_output_mm10/D1_1

# submit just D1_1
cd ~/links/scratch/BoneMarrow_miRNA_EV_QC
sbatch exceRpt_job_mm10.sh NX.VH01519_298.001.RPI1.D1_1_R1_bbduk_min16.fastq.gz D1_1

# check status
squeue -u $USER
```

Wait for it to finish (≈13 min for this dataset), then check the stats
(see PART 6 for how to read it):

```bash
cat ~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_output_mm10/D1_1/*.stats
```

Only proceed to all 62 if the test produced a sensible stats file with
miRNAs detected.

---

# PART 4 — Create the launcher and submit all 62

```bash
cat > ~/links/scratch/BoneMarrow_miRNA_EV_QC/launch_exceRpt_all.sh << 'EOF'
#!/bin/bash
# Launches one exceRpt SLURM job per trimmed sample (62 total)

INPUT_DIR=~/links/scratch/BoneMarrow_miRNA_EV_QC/trimmed
OUTPUT_DIR=~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_output_mm10
JOB=~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_job_mm10.sh

mkdir -p $OUTPUT_DIR

for fq in $INPUT_DIR/*_bbduk_min16.fastq.gz; do
    filename=$(basename "$fq")
    sample_id=$(echo "$filename" | grep -oP 'D[0-9]+_[0-9]+')
    echo "Submitting $sample_id  ($filename)"
    mkdir -p $OUTPUT_DIR/$sample_id
    sbatch $JOB "$filename" "$sample_id"
done
EOF
```

Then run it (with **bash**, not sbatch — it submits the jobs):

```bash
bash ~/links/scratch/BoneMarrow_miRNA_EV_QC/launch_exceRpt_all.sh
```

> **Important — sample naming:** the launcher extracts the `D#_#` sample ID
> with `grep -oP 'D[0-9]+_[0-9]+'`. The muscle launcher used
> `${filename%%_*}`, which would NOT work on these filenames (it would name
> every sample `NX`). The regex approach gives each sample its correct
> unique name (D1_1, D6_12, etc.).

---

# PART 5 — Monitor the jobs

```bash
# how many exceRpt jobs still in the queue
squeue -u $USER | grep exceRpt | wc -l
```

- A number → still running/queued
- `0` → all finished

You can disconnect; jobs run on their own. Email notifications fire on
start/end of each.

When the queue is empty, confirm all 62 succeeded:

```bash
ls ~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_output_mm10/*/*.stats | wc -l
```

Should print `62`. If fewer, rerun the missing samples (resubmit just those
with the job script).

Quick miRNA check across all samples:

```bash
grep -H "^miRNA_sense" ~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_output_mm10/*/*.stats
```

---

# PART 6 — Reading a .stats file

The stats file tracks reads through sequential stages. Example (D1_1):

| Stage | Meaning |
|-------|---------|
| `input` | reads that went in |
| `successfully_clipped` | NA (adapter trimming skipped — already done) |
| `failed_quality_filter` | reads dropped for low quality |
| `failed_homopolymer_filter` | reads dropped as homopolymer artifacts (AAAA…) |
| `calibrator` | NA (no spike-in used) |
| `UniVec_contaminants` | reads matching common lab vector/primer contamination |
| `rRNA` | ribosomal RNA reads (removed; ~9% is normal) |
| `reads_used_for_alignment` | reads surviving filters, sent to alignment |
| `genome` | reads mapping to mm10 genome generally |
| **`miRNA_sense`** | **reads mapping to mature miRNAs — the key signal** |
| `miRNA_antisense` | wrong-strand (noise; should be ~0) |
| `miRNAprecursor_sense` | reads on pre-miRNA hairpins (not mature) |
| `tRNA_sense` | tRNA fragments (often abundant in EVs) |
| `piRNA_sense` | piwi-interacting RNAs (small) |
| `gencode_sense` | other annotated transcripts (mRNA frag, lncRNA…) |
| `circularRNA_sense` | circular RNAs (usually ~0) |
| `not_mapped_to_genome_or_libs` | reads matching nothing (often high in EVs) |

`_sense` = same orientation as the annotated RNA (real signal);
`_antisense` = opposite (noise).

> **Note on EV data:** a high `not_mapped` fraction and a relatively low
> miRNA fraction are common in extracellular vesicle small RNA, partly
> because EV RNA is fragmented and partly because of the strict zero-mismatch
> setting. What matters is consistency across samples and comparability to
> the muscle dataset — not any single number.

---

# PART 7 — Merge into one count matrix

Once all 62 exceRpt jobs are done, combine the per-sample outputs into single
matrices (miRNA × samples, plus tRNA, piRNA, etc.) using
`mergePipelineRuns.R` from the shared folder, run inside the container.

The merge script is invoked as:
`Rscript mergePipelineRuns.R <data path> <output path>` — it scans all sample
folders under the data path and writes the combined matrices.

**First confirm all 62 finished:**
```bash
squeue -u $USER | grep exceRpt | wc -l                                            # should be 0
ls ~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_output_mm10/*/*.stats | wc -l   # should be 62
```

**Create the merge job script:**
```bash
cat > ~/links/scratch/BoneMarrow_miRNA_EV_QC/merge_exceRpt.sh << 'EOF'
#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o merge_exceRpt_%j.log

module load StdEnv/2020
module load apptainer/1.1.8

SHARED=~/links/projects/rrg-gturecki/shared/Software_and_Installation
DATA_DIR=~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_output_mm10
MERGE_DIR=~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_merged

mkdir -p $MERGE_DIR

apptainer exec \
  -B $DATA_DIR:/exceRptOutput \
  -B $MERGE_DIR:/mergedOutput \
  -B $SHARED/exceRpt:/scripts \
  $SHARED/exceRpt/exceRpt.sif \
  Rscript /scripts/mergePipelineRuns.R /exceRptOutput /mergedOutput
EOF
```

**Submit it:**
```bash
cd ~/links/scratch/BoneMarrow_miRNA_EV_QC
sbatch merge_exceRpt.sh
squeue -u $USER
```

### How the merge works

| Piece | Meaning |
|-------|---------|
| `apptainer exec` | runs a specific command (`Rscript`) inside the container, instead of exceRpt's default pipeline |
| `-B $DATA_DIR:/exceRptOutput` | makes the 62 sample folders visible inside the container |
| `-B $MERGE_DIR:/mergedOutput` | where the combined matrices are written |
| `-B $SHARED/exceRpt:/scripts` | makes the merge `.R` scripts visible inside the container |
| `Rscript .../mergePipelineRuns.R /exceRptOutput /mergedOutput` | runs the merge: data path, then output path |

The merge is much lighter than alignment and finishes quickly.

### Key output (in `exceRpt_merged/`)

- **`exceRpt_miRNA_ReadCounts.txt`** — the **miRNA × 62-sample count matrix**
  (feeds DESeq2; compare distinct-miRNA count to muscle's **468**)
- equivalent matrices for tRNA, piRNA, gencode, etc.
- a combined QC summary across all samples

**Check the matrix after merging:**
```bash
# number of distinct miRNAs detected (rows minus header)
wc -l ~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_merged/exceRpt_miRNA_ReadCounts.txt

# confirm sample columns + peek at the top
head ~/links/scratch/BoneMarrow_miRNA_EV_QC/exceRpt_merged/exceRpt_miRNA_ReadCounts.txt
```

---

# PART 8 — Next step: DESeq2

The `exceRpt_miRNA_ReadCounts.txt` matrix is the input for differential
expression in R with DESeq2 (sedentary vs exercise), including batch as a
covariate (samples split across two sequencing pools: batch 1 = the 298.001
files, batch 2 = the 299.001 files). Set up separately. Compare the detected
miRNA count and downstream results against the muscle dataset for the
cross-tissue analysis.

---

# Quick reference — running it again

```bash
# 1. SSH into a login node
ssh reinekh@rorqual.alliancecan.ca

# 2. Test one sample
cd ~/links/scratch/BoneMarrow_miRNA_EV_QC
mkdir -p exceRpt_output_mm10/D1_1
sbatch exceRpt_job_mm10.sh <one_file>.fastq.gz D1_1

# 3. If good, launch all 62
bash launch_exceRpt_all.sh

# 4. Monitor
squeue -u $USER | grep exceRpt | wc -l   # 0 = done

# 5. Confirm 62 stats files
ls exceRpt_output_mm10/*/*.stats | wc -l

# 6. Merge into one matrix
sbatch merge_exceRpt.sh
```

---

# Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| All samples named the same / overwrite | launcher used `${filename%%_*}` | use `grep -oP 'D[0-9]+_[0-9]+'` |
| Job fails instantly | wrong `--account` or missing output folder | check account; `mkdir` the sample folder |
| Container not found | wrong shared path | verify `exceRpt.sif` path |
| `<jobid>` syntax error in tail | typed the `< >` brackets | use just the number, no brackets |
| Fewer than 62 stats files | some jobs failed | rerun the missing samples individually |
| Merge can't find functions script | scripts folder not bound | ensure `-B $SHARED/exceRpt:/scripts` is present |

---

*Reference for the bone marrow miRNA-EV project. Tool: exceRpt v4.6.3 via
Apptainer. Parameters match the muscle dataset for cross-tissue comparison.*
