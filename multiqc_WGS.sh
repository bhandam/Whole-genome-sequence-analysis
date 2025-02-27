#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --output=slurm-multiqc-%j.out

## Bash strict settings, stop if there is error, if varibales cannot be found
set -euo pipefail

##load our conda environement
module load python
source activate multiqc-env

## process aruguments
indir=$1
outdir=$2

echo "The input directory is $indir"
echo "The output directory file is $outdir"
##create output directory, if needed
mkdir -p "$outdir"
#-p allow to create multiple level at once
multiqc --interactive --force "$indir" -o "$outdir"

# --force: overwrite old file

echo "Done with script"
date
