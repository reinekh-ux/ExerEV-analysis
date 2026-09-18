# Quality Control and Adapter Trimming of Bone Marrow miRNA-EV Libraries

**Dataset:** 62 bone marrow small extracellular vesicle (sEV) small RNA
libraries, single-end, 101 bp, sequenced in two pools (batches).
**Library prep:** Galas Lab 4N small RNA protocol (4 random nucleotides at
each ligation junction to reduce ligation bias).
**Processing environment:** Alliance (Compute Canada) Rorqual cluster.

This document describes the three-stage workflow: (1) initial quality
control of the raw reads, (2) adapter and quality trimming, and (3) a
second round of quality control on the trimmed reads to confirm the
trimming worked.

---

## Overview of the workflow

```
Raw FASTQ (62 samples)
        │
        ▼
[1] FastQC + MultiQC  ──►  assess raw read quality, confirm library structure
        │
        ▼
[2] BBDuk trimming    ──►  remove 5' 4N, 3' adapter, barcode, junk tail; length filter
        │
        ▼
[3] FastQC + MultiQC  ──►  confirm clean ~21–22 nt miRNA length distribution
        │
        ▼
   Clean reads ready for miRNA quantification
```

The same pipeline (identical parameters) was used to process a parallel
muscle miRNA-EV dataset, so the two tissues can be compared without
processing-related confounds.

---

# Stage 1 — Initial quality control (raw reads)

## Purpose

Before any processing, the raw reads are assessed to (a) confirm the
sequencing succeeded, (b) understand the library structure, and (c)
identify what needs to be trimmed.

## Method

FastQC (v0.12.1) was run on each of the 62 raw FASTQ files, producing one
quality report per sample. MultiQC (v1.33) was then used to aggregate all
62 reports into a single interactive summary so that all samples could be
compared on the same axes.

## What the raw QC showed

The raw reports displayed the expected signature of a 4N small RNA
library:

- **Read depth:** all 62 samples sequenced successfully, ranging roughly
  10–36 million reads — adequate depth for small RNA quantification.
- **Per-base quality:** high quality (Phred ~38–40) for the first ~60 bp,
  then a sharp decline toward the read end. This is expected: the miRNA
  insert is short (~21–22 nt), so the latter part of each 101 bp read is
  adapter and low-quality read-through.
- **Adapter content:** all 62 samples flagged for adapter content, with
  the Illumina small RNA 3' adapter (TGGAATTCTCGGGTGCCAAGG) rising from
  ~15 bp onward. This is the hallmark of a working small RNA library —
  the read passes through the short insert and into the adapter early.
- **Duplication / overrepresented sequences:** high duplication and
  overrepresented sequences, dominated by adapter and the most abundant
  miRNAs. Expected for small RNA, where a small number of molecules make
  up most of the reads.
- **Per-base N content:** essentially 0% across the read, indicating clean
  base calls.

## Interpretation

Every flag FastQC raised is the normal, expected behaviour of a small RNA
library — not a quality problem. The raw QC confirmed the libraries were
successful and that the next required step was adapter and quality
trimming.

## Read structure identified

Inspection of the raw reads confirmed the following structure (single-end,
read from the 5' side):

```
[5' 4N] [====miRNA insert====] [3' 4N] [3' adapter] [index] [low-quality tail]
```

- **5' 4N:** 4 random bases from the 5' adapter (Galas 4N design)
- **miRNA insert:** the biological molecule of interest (~18–24 nt)
- **3' 4N:** 4 random bases from the 3' adapter
- **3' adapter:** TGGAATTCTCGGGTGCCAAGG (Illumina small RNA 3' adapter)
- **index:** the RPI barcode (e.g. ATCACG for RPI1); present in the read
  but downstream of the adapter
- **tail:** low-quality read-through to fill the 101 bp read length

Only the miRNA insert is biologically meaningful; everything else is
protocol scaffolding to be removed.

---

# Stage 2 — Adapter and quality trimming

## Purpose

Remove all non-biological sequence (the two 4N blocks, the 3' adapter, the
barcode, and the low-quality tail) so that only the miRNA insert remains.
Reads that contain no real insert (adapter dimers) or are too short to be
miRNAs are discarded.

## Tool

BBDuk and reformat.sh from BBMap (v38.86), matching the muscle dataset
processing. The adapter to trim was supplied as a FASTA file:

```
>3_prime_adapter
NNNNTGGAATTCTCGGGTGCCAAGG
```

The leading `NNNN` instructs BBDuk to account for the 3' 4N random bases
together with the fixed adapter.

## The three trimming steps

**Step 1 — Hard cap at 50 bp (`forcetrimright=50`).**
Each read is truncated to its first 50 bp. The miRNA insert plus 4N plus
adapter all fall well within the first 50 bp, so no biological signal is
lost; this removes the worst of the low-quality tail up front and matches
the muscle dataset's read-length handling.

**Step 2 — Adapter trimming (BBDuk, `ktrim=r`).**
- `forcetrimleft=4` removes the 4 random bases at the 5' end (the 5' 4N).
- `ref=...NNNNTGGAATTCTCGGGTGCCAAGG`, `ktrim=r` finds the 3' adapter and
  removes it and everything to its right — which, because the adapter sits
  between the insert and the barcode, also removes the 3' 4N, the barcode,
  and the junk tail in a single operation.
- Supporting parameters: `k=8` (8-base seed for adapter matching),
  `minoverlap=6` (trim even with short adapter overlap), `mininsert=16`
  (require at least 16 bp of insert), `rcomp=f` (do not search the reverse
  complement), `copyundefined` (expand the N bases when matching).

**Step 3 — Length filter (`reformat.sh minlength=16`).**
Reads shorter than 16 bp after trimming are discarded. Mature miRNAs are
~18–24 nt, so anything below 16 nt is not a miRNA (typically adapter
dimers with no insert).

## After both ends are trimmed

```
Before:  [5'4N][====miRNA====][3'4N][adapter][index][tail]   (101 → 50 bp)
                  │
                  ▼  forcetrimleft=4 (5' side) + ktrim=r adapter (3' side)
After:   [====miRNA====]                                      (~18–24 nt)
```

## Single-sample test result (sample D1_1, RPI1)

The pipeline was validated on one sample before running all 62:

| Stage | Reads | Note |
|-------|-------|------|
| Raw input | 16,321,496 | 101 bp reads |
| After 50 bp cap | 16,321,496 | all retained, now 50 bp |
| After adapter trim | 11,843,329 | adapter found in 99.69% of reads; 27.4% of reads removed as dimers |
| After length filter (≥16 bp) | 7,443,950 | final clean reads (~46% of raw) |

Inspection of individual trimmed reads confirmed they became short,
insert-only sequences in the small RNA range (e.g. 16–19 nt), with the
adapter, 4N, barcode, and tail all removed. A retention of ~46% is normal
for 4N small RNA libraries, where adapter dimers and too-short fragments
are expected to be discarded.

> **Note on read lengths after trimming:** trimmed reads vary in length
> (e.g. 16, 19, 22 nt) because different miRNAs are naturally different
> lengths. The uniform 101 bp of the raw reads was an artifact of fixed
> sequencing length; trimming reveals each molecule's true size. A spread
> in the ~18–24 nt range is the expected, correct result.

## Batch processing

All 62 samples were processed with the identical three-step pipeline,
submitted as a single SLURM batch job looping over every FASTQ file. The
final output per sample is `<sample>_bbduk_min16.fastq.gz`.

---

# Stage 3 — Post-trimming quality control

## Purpose

Confirm that trimming worked across all 62 samples — that the adapter is
gone and the reads now show the characteristic miRNA length distribution.

## Method

FastQC and MultiQC are run again, this time on the 62 trimmed
(`_bbduk_min16.fastq.gz`) files, and the aggregated report is compared
against the raw QC report.

## What a successful post-trim report should show

- **Read-length distribution:** a sharp peak around **21–22 nt** (the
  dominant mature-miRNA size), with a spread from ~18–25 nt. This is the
  single most important confirmation — the appearance of this peak means
  the adapter was removed and intact miRNAs remain.
- **Adapter content:** the adapter flag should now be cleared (or
  drastically reduced) — the adapter sequence should no longer be present.
- **Per-base quality:** improved overall, since the low-quality tail has
  been removed.
- **Read counts:** reduced relative to raw (adapter dimers and short
  fragments removed), as expected.

## Interpretation

If the trimmed report shows the ~21–22 nt peak and a cleared adapter flag,
the libraries are clean and ready for the next stage — miRNA
identification and quantification (e.g. alignment/annotation against a
miRNA reference such as miRBase).

---

# Summary

| Stage | Tool | Input | Output | Confirms |
|-------|------|-------|--------|----------|
| 1. Raw QC | FastQC + MultiQC | 62 raw FASTQ | 1 aggregated report | sequencing succeeded; library structure; trimming needed |
| 2. Trimming | BBDuk + reformat.sh (BBMap 38.86) | 62 raw FASTQ | 62 trimmed FASTQ | adapter/4N/barcode/tail removed; inserts isolated |
| 3. Post-trim QC | FastQC + MultiQC | 62 trimmed FASTQ | 1 aggregated report | clean ~21–22 nt miRNA distribution; adapter cleared |

The identical pipeline (BBMap 38.86, Galas 4N adapter, `forcetrimleft=4`,
`forcetrimright=50`, `ktrim=r`, `mininsert=16`, `minlength=16`) was used
for the parallel muscle miRNA-EV dataset, enabling direct cross-tissue
comparison of the exercise vs. sedentary response without processing
confounds. The only difference is that the muscle data were sequenced
paired-end (PE100) but analyzed on R1 only (`skipr2=t`), while the bone
marrow data are single-end; the biologically meaningful trimming is
identical.

---

*Prepared as a workflow/methods reference. Tool versions: FastQC v0.12.1,
MultiQC v1.33, BBMap v38.86. Adjust sample-specific details as needed.*
