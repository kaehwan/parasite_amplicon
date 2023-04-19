#!/usr/bin/env bash

# please run dereplicate_fasta.py to remove replicate sequences from input fasta

SECONDS=0

eval "$(conda shell.bash hook)"
conda activate muscle

input=$1

echo "MSA-ing ${input}..."

muscle \
	-super5 ${input} \
	-output ${input%.fasta}.aligned.fasta \
	-nt \
	-threads 50

echo "DONE MSA-ing"

DURATIONS=$SECONDS
echo "$(($DURATIONS / 60)) minutes and $(($DURATIONS % 60)) seconds elapsed."
