#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --output=slurm-fastqc-%j.out

## Bash strict settings, stop if there is error, if varibales cannot be found
set -euo pipefail

echo "This script will run FastQC"
date

## process the command line argument
fastq_file=$1
outdir=$2

echo "The FASTQ file is $fastq_file"
echo "The output directory file is $outdir"

# If you do not have the output directory then make the directory
# that is passed as an argumnet for the code.
mkdir -p $outdir

#load software module
module load fastqc

#Run fastQC
fastqc --outdir=$outdir $fastq_file

echo "Done with script"
date

#tail -f: to see the end of file
