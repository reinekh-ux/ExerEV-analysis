#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=35G
#SBATCH -o run_exceRpt_jobs_mm10_%j.log
#SBATCH --mail-user=minh.nguyen.comtl@ssss.gouv.qc.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020
module load apptainer/1.1.8

echo $1
echo $2

apptainer run -B /scratch/minhng/Reine_miRNA/bbduk_output/bbduk_step3_trimmed_min16:/exceRptInput \
		-B exceRpt_output_mm10/$2:/exceRptOutput \
		-B /project/rrg-gturecki/Software_Installations_and_Databases/exceRpt_mm10//mm10:/exceRpt_DB/mm10 \
		/project/rrg-gturecki/Software_Installations_and_Databases/exceRpt/exceRpt.sif \
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
