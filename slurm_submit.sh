#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=rnaseq.%j.log
##SBATCH --partition=your_partition
##SBATCH --account=your_account
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=you@example.com
#
# Example SLURM job for rnaseq_pipeline.sh. Edit the paths below, then:
#   sbatch slurm_submit.sh

set -euo pipefail

# ------------------------------- edit me -------------------------------------
CONDA_ENV=rnaseq                                   # env name or full prefix path
PIPELINE=/path/to/RNAseq_pipeline/rnaseq_pipeline.sh
INPUT=/path/to/fastq_dir                           # flat or one folder per sample
OUTDIR=/path/to/results
STAR_INDEX=/path/to/STAR_index
GTF=""                                             # empty = use the GTF recorded in the STAR index
STRANDEDNESS=0                                     # 0 unstranded, 1 stranded, 2 reverse
# -----------------------------------------------------------------------------

# If conda is provided by an environment module on your cluster, load it first, e.g.:
# module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

bash "$PIPELINE" \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    --star-index "$STAR_INDEX" \
    ${GTF:+--gtf "$GTF"} \
    --strandedness "$STRANDEDNESS" \
    --threads "${SLURM_CPUS_PER_TASK:-16}"
