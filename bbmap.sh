#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --output=slurm-bbmap-%j.out

## Software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/bbmap-38.96

## Bash strict settings
set -euo pipefail

## Arguments
R1=$1
R2=$2
assembly=$3
outfile=$4

## Make output dir
outdir=$(dirname $outfile)
mkdir -p "$outdir"

## Report
echo "Running script bbmap.sh"
date
echo "R1:                 $R1"
echo "R2:                 $R2"
echo "assembly:           $assembly"
echo "output file:        $outfile"
echo

# To index our assembly for that same sample
##echo "Now indexing the assembly..."
##bbmap.sh ref="$assembly"

# Map the FASTQs to the assembly
echo "Now mapping..."
bbmap.sh in="$R1" in2="$R2" ref="$assembly" t="$SLURM_CPUS_PER_TASK" out="$outfile" nodisk

# determine the coverage
pileup.sh in= out=stats.txt

## Report
echo "Done with script"
date
