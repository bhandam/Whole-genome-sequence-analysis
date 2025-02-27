#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --output=slurm-snippy-%j.output
#SBATCH --cpus-per-task=16

module load python/3.6-conda5.2
source activate /users/PAS0471/menuka/.conda/envs/snippy

## Arguments
cpus=$1
outdir=$2
ref=$3
ctgs=$4

mkdir -p $outdir

snippy --cpus 16 --outdir "$outdir" --ref "$ref" --ctgs "$ctgs"