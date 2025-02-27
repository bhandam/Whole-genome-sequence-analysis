#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --output=slurm-roary-%j.output
#SBATCH --cpus-per-task=16

#load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/roary-3.13
echo "## starting script roary.sh"
date
echo

## Bash strict settings
set -euo pipefail

## Arguments
outdir=$1
indir=$2

## Report
echo "Outdir: $outdir"
echo "Input: $indir"

mkdir -p "$outdir"

## Run roary
roary -f "$outdir" \
    -e -n -p "$SLURM_CPUS_PER_TASK" \
    -v "$indir"/*gff

## Report
echo "## Done with script"
