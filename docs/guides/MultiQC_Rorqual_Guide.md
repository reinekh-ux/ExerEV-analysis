# Running FastQC + MultiQC on Rorqual — Step-by-Step Guide

A beginner-friendly walkthrough for QC of sequencing libraries on the
Alliance (Compute Canada) Rorqual cluster. Written for someone with no
prior command-line experience.

---

## What this pipeline does

1. **FastQC** runs on each FASTQ file separately and produces one quality
   report per sample.
2. **MultiQC** scans all those FastQC reports and combines them into a
   single HTML page so you can compare every sample at once.

**Input:** your `.fastq.gz` files
**Output:** one combined `multiqc_report.html`

---

## The single most important lesson

> **Do the setup on a LOGIN NODE, not in the JupyterHub terminal.**

The JupyterHub terminal sets a hidden `PYTHONPATH` variable that breaks
Python virtual environments and makes MultiQC impossible to install. A
login node has a clean environment where the standard install just works.

You reach a login node by SSH-ing in from your own computer's terminal.

---

## Key paths (edit these to match your project)

| What | Path |
|------|------|
| Raw data | `~/links/projects/rrg-gturecki/raw_files/20260617_Reine_miRNAEv` |
| Working folder (scratch) | `~/links/scratch/BoneMarrow_miRNA_EV_QC` |
| Compute account | `def-cnagy` |

> Note: `~/links/scratch` points to your real scratch
> (`/lustre10/scratch/reinekh`). Always use `~/links/scratch`, not
> `~/scratch`.

---

# PART 1 — One-time setup (install MultiQC)

You only do this once. After it's done, the environment lives in your
scratch and you just reuse it.

### Step 1.1 — Connect to a login node

Open **Terminal** on your Mac (Cmd+Space, type "Terminal", Enter), then:

```bash
ssh reinekh@rorqual.alliancecan.ca
```

Enter your Alliance password (you won't see it as you type) and complete
multifactor authentication if prompted. You'll know you're connected when
the prompt looks like `[reinekh@rorqual2 ~]$`.

### Step 1.2 — Clean the environment

```bash
unset PYTHONPATH
export PYTHONNOUSERSITE=1
```

This removes the leftover variable that would otherwise break the install.

### Step 1.3 — Load the software modules

```bash
module load python/3.11 gcc arrow
```

- `python/3.11` — the Python interpreter
- `gcc` + `arrow` — provide `pyarrow`, a dependency MultiQC needs

> The order matters and `arrow` must be loaded BEFORE creating the
> virtual environment.

### Step 1.4 — Create a virtual environment

A "virtual environment" is a private, self-contained folder for Python
packages so they don't conflict with anything else.

```bash
virtualenv --no-download ~/links/scratch/multiqc_env
```

`--no-download` tells it to use Alliance's pre-built packages instead of
the internet.

### Step 1.5 — Activate it

```bash
source ~/links/scratch/multiqc_env/bin/activate
```

Your prompt now starts with `(multiqc_env)`. This means the environment
is active.

### Step 1.6 — Confirm pyarrow is visible (important check)

```bash
python -c "import pyarrow; print(pyarrow.__version__)"
```

This should print a version number like `24.0.0`. If it errors with
"No module named pyarrow", stop — the modules in Step 1.3 didn't load
correctly. Re-run Step 1.3, then re-activate (Step 1.5).

### Step 1.7 — Install MultiQC

```bash
pip install --no-index --upgrade pip
pip install --no-index multiqc
```

`--no-index` means "use Alliance's tested packages, not the internet."

### Step 1.8 — Confirm it works

```bash
multiqc --version
```

Should print `multiqc, version 1.33` (or similar). Setup is done.

---

# PART 2 — Running the QC (do this each time)

### Step 2.1 — Make your working folders

```bash
mkdir -p ~/links/scratch/BoneMarrow_miRNA_EV_QC/fastqc
mkdir -p ~/links/scratch/BoneMarrow_miRNA_EV_QC/multiqc
```

### Step 2.2 — Test on ONE sample first (always do this)

Running one file first confirms the pipeline works before committing
hours to all of them.

```bash
module load fastqc
fastqc ~/links/projects/rrg-gturecki/raw_files/20260617_Reine_miRNAEv/NX.VH01519_298.001.RPI1.D1_1_R1.fastq.gz \
  -o ~/links/scratch/BoneMarrow_miRNA_EV_QC/fastqc -t 2
```

When you see `Analysis complete for ...`, it worked.

### Step 2.3 — Find your compute account

```bash
sshare -U -u $USER | head
```

Look for an account like `def-cnagy`. Use the part before `_cpu`.

### Step 2.4 — Confirm your file count

```bash
ls ~/links/projects/rrg-gturecki/raw_files/20260617_Reine_miRNAEv/*.fastq.gz | wc -l
```

Should print the number of samples you expect (e.g. `62`).

### Step 2.5 — Create the batch job script

A "batch job" is a recipe you hand to the cluster's scheduler so it runs
on a powerful compute node while you do other things (even with your
laptop closed).

Paste this whole block to write the script:

```bash
cat > ~/links/scratch/BoneMarrow_miRNA_EV_QC/run_qc.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=miRNA_QC
#SBATCH --account=def-cnagy
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=%x_%j.out

module load fastqc python/3.11 gcc arrow

RAW=~/links/projects/rrg-gturecki/raw_files/20260617_Reine_miRNAEv
QC=~/links/scratch/BoneMarrow_miRNA_EV_QC
FASTQC_OUT=$QC/fastqc
MULTIQC_OUT=$QC/multiqc

mkdir -p $FASTQC_OUT $MULTIQC_OUT

# Run FastQC on all samples (8 threads in parallel)
fastqc $RAW/*.fastq.gz -o $FASTQC_OUT -t 8

# Activate MultiQC and aggregate all reports into one
source $QC/multiqc_env/bin/activate
multiqc $FASTQC_OUT -o $MULTIQC_OUT --title "BoneMarrow miRNA-EV QC" --force
EOF
```

**What each `#SBATCH` line means:**

| Line | Meaning |
|------|---------|
| `--job-name=miRNA_QC` | a label to spot your job in the queue |
| `--account=def-cnagy` | charges compute to your lab's allocation |
| `--time=03:00:00` | reserve up to 3 hours |
| `--cpus-per-task=8` | request 8 CPU cores |
| `--mem=16G` | request 16 GB memory |
| `--output=%x_%j.out` | save a log file named after the job |

### Step 2.6 — Check the script saved correctly

```bash
cat ~/links/scratch/BoneMarrow_miRNA_EV_QC/run_qc.sh
```

It should print the script back, starting with `#!/bin/bash`.

### Step 2.7 — Submit the job

```bash
cd ~/links/scratch/BoneMarrow_miRNA_EV_QC
sbatch run_qc.sh
```

It prints `Submitted batch job XXXXXXXX`. That number is your job ID.

### Step 2.8 — Check job status

```bash
squeue -u $USER
```

- `PD` = pending (waiting in queue)
- `R` = running
- gone from the list = finished

You can now close your laptop or disconnect — the job keeps running.

### Step 2.9 — Watch the live log (optional)

```bash
tail -f ~/links/scratch/BoneMarrow_miRNA_EV_QC/miRNA_QC_<jobID>.out
```

Replace `<jobID>` with your number. Press `Ctrl+C` to stop watching
(this does NOT stop the job).

---

# PART 3 — After the job finishes

### Step 3.1 — Find your report

```
~/links/scratch/BoneMarrow_miRNA_EV_QC/multiqc/BoneMarrow-miRNA-EV-QC_multiqc_report.html
```

### Step 3.2 — View it

Open it through the **JupyterHub Files panel** (scratch is visible there):
navigate to `links / scratch / BoneMarrow_miRNA_EV_QC / multiqc` and click
the HTML file.

### Step 3.3 — Back it up (important!)

Scratch is **purged after ~60 days and is not backed up.** Copy the final
report somewhere permanent:

```bash
cp ~/links/scratch/BoneMarrow_miRNA_EV_QC/multiqc/*.html \
   ~/links/projects/rrg-gturecki/
```

(Adjust the destination to wherever you keep results in your project space.)

---

# Quick reference — running it again from scratch

Once setup (Part 1) is done once, a fresh run is just:

```bash
# 1. SSH into a login node
ssh reinekh@rorqual.alliancecan.ca

# 2. Submit the job (script already exists)
cd ~/links/scratch/BoneMarrow_miRNA_EV_QC
sbatch run_qc.sh

# 3. Check status
squeue -u $USER
```

---

# Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| `No module named multiqc` after install | You're in the JupyterHub terminal; PYTHONPATH broke the install | Do it from an SSH login node instead |
| `No module named pyarrow` | Arrow module not loaded, or env activated before loading it | `module load python/3.11 gcc arrow`, then re-activate the env |
| pip tries the "dummy pyarrow" wheel and errors | PYTHONPATH unset so pyarrow invisible | Load the arrow module so pyarrow is visible before installing |
| `i/o timeout` reaching Docker / the internet | Compute nodes have no internet | Only download things on a login node |
| Job fails instantly | Wrong `--account` name | Check with `sshare -U -u $USER` |
| Can't find scratch | `~/scratch` is wrong | Use `~/links/scratch` |

---

*Generated as a personal reference. Adjust paths, account, and sample
names to match your own project.*
