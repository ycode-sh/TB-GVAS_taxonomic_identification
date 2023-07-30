#!/bin/bash

# This script is relevant to deconntaminating a bam file. It relies on thre fact that reads that are too distant from the reference should 
# not map and therefore are flagged as unmapped. The scrip that uses the -F option to include only mapped reads in the fastq file(s) that are
#subsequently generated.    $1=bamfile, $2=platform

filename=$(basename $1 _sam.bam)

if [[ $2 == "ont" ]]; then
    samtools fastq -0 "${filename}.fastq" -F 4 $1

else
    samtools fastq -1 "${filename}.fastq_1" -2 "${filename}.fastq_2" -F 4 $1

fi