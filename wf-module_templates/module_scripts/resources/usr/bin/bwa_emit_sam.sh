#!/bin/bash
# raw_fastq1=$1; raw_fastq2=$2; fastafile=$3; sample_ID =$4


# Index fasta file
bwa index $3

# Align  reads to reference
bwa mem -R "@RG\tID:run_$4\tPL:ILLUMINA\tLB:paired_end\tSM:$4" $3 $1 $2 -o "${4}_sam"