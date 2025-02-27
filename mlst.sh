#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --output=slurm-mlst-%j.output
#SBATCH --cpus-per-task=16

#load software
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/mlst
echo "## starting script mlst.sh"
date
echo

## Bash strict settings
set -euo pipefail

## Arguments
indir=$1
outfile=$2


## Report
echo "Input: $indir"
echo "Outfile: $outfile"


## Run mlst
mlst "$indir" --legacy --scheme senterica_achtman_2 --csv "$outfile"

## Report
echo "## Done with script"
