#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --output=slurm-quast-%j.output
#SBATCH --cpus-per-task=16

#load software
module load python
source activate /fs/project/PAS0471/jelmer/conda/quast-5.0.2

echo "## starting script quast.sh"
date
echo

## Bash strict settings
set -euo pipefail

## Arguments
outdir=$1


# If you do not have the output directory then make the directory
# that is passed as an argumnet for the code.
mkdir -p $outdir



quast.py\
    -o "$outdir" \
    -t 16 \
    --conserved-genes-finding \
    results_29_new_seq/spades/GRS-S-1/contigs.fasta

#-g ref.gff \
#-r ref.fna \

echo "## Done with script"
date