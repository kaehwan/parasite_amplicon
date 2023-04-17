#!/usr/bin/env bash

# please run dereplicate_fasta.py to remove replicate sequences from input fasta

SECONDS=0

eval "$(conda shell.bash hook)"
conda activate clustalw

input=$1

echo "MSA-ing ${input}..."

clustalw \
	-infile=${input} \
	-type=dna \
	-pwgapopen=15.00 \
	-pwgapext=6.66 \
	-gapopen=15.00 \
	-gapext=6.66 \
	-maxdiv=30 \
	-endgaps \
	-transweight=0.50 \
	-output=clustal \
	-outorder=aligned \
	-iteration=alignment

echo "DONE MSA-ing"

DURATIONS=$SECONDS
echo "$(($DURATIONS / 60)) minutes and $(($DURATIONS % 60)) seconds elapsed."
