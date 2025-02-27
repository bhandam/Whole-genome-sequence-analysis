#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --output=slurm-seqserro2-%j.output
#SBATCH --cpus-per-task=4

# Arguments
R1=$1
outdir=$2

#load software
module load python/3.6-conda5.2
source activate Seqsero2-env

## Infer the name of the R2 file
R2=${R1/_R1_/_R2_}
outdir_full="$outdir"/$(basename "$R1" _R1_001.fastq.gz)

echo "## starting script seqsero2.sh"
date
echo
echo "## R1 file is:      $R1"
echo "## R2 file is:      $R2"
echo "## Outdir is:       $outdir_full"

mkdir -p "$outdir_full"

SeqSero2_package.py \
    -t 2 \
    -i "$R1" "$R2" \
    -d "$outdir_full" \
    -p $SLURM_CPUS_PER_TASK

echo "## Done with script"
date