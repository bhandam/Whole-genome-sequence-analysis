#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-coverage-%j.out

## Software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/bbmap-38.96

## Bash strict settings
set -euo pipefail

## Arguments
input=$1
outfile=$2
#sample_ID=$(basename "$input")


## Report
echo "Running script pileup.sh"
date
echo "R1:                 $input"
echo "outfile:                 $outfile"
echo

# determine the coverage
pileup.sh in="$input" out="$outfile"

## Report
echo "Done with script"
date
