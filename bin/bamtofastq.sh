#!/bin/bash

# This script is relevant to deconntaminating a bam file. It relies on thre fact that reads that are too distant from the reference should 
# not map and therefore are flagged as unmapped. The scrip that uses the -F option to include only mapped reads in the fastq file(s) that are
#subsequently generated.    $1=bamfile, $2=platform

filename=$(basename $1 _sam.bam)

if [[ $2 == "minion_ont_reads" || $2 == "se_illumina_reads" ]]; then
    samtools fastq -0 "${filename}.fastq" -F 4 $1

elif [[ $2 == "pe_illumina_reads" ]]; then
    samtools fastq -1 "${filename}_1.fastq" -2 "${filename}_2.fastq" -F 4 $1

fi