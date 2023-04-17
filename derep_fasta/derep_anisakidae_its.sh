#!/usr/bin/env bash

clustalw \
-infile=derep_anisakidae_coi5p.fasta \
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
