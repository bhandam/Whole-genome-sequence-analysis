#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --output=slurm-abricate-%j.output
#SBATCH --cpus-per-task=16

#load software
source activate abricate

echo "## starting script abricate.sh"
date
echo

## Bash strict settings
set -euo pipefail

## Arguments
input=$1
outfile=$2
db=$3

## Create output dir if needed
outdir=$(dirname "$outfile")
mkdir -p "$outdir"

## Run abricate
abricate --db "$db" "$input" > "$outfile"

## Report
echo "## Done with script"
