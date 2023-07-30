#!/bin/bash

# This script takes ONT minion reads and a fasta reference, indexes fasta reference and aligns the reads to the reference supplied 
# $1 = fastq reads, $2 = fasta reference

sample_id=$(basename $1 .fastq)

minimap2 -d $2 h37rv.mmi
minimap2 -a h37rv.mmi $1 -R @RG\\tID:$sample_id\\tSM:$sample_id -o "${sample_id}_sam"