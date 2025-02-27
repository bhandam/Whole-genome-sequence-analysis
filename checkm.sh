#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --output=slurm-checkm-%j.output
#SBATCH --cpus-per-task=16

#load software
module load python
source activate /fs/project/PAS0471/jelmer/conda/checkm-1.2.0

echo "## starting script checkm.sh"
date
echo

## Bash strict settings
set -euo pipefail

## Arguments

input=$1
outdir=$2

# If you do not have the output directory then make the directory
# that is passed as an argumnet for the code.
mkdir -p "$outdir"

checkm lineage_wf "$input" "$outdir"