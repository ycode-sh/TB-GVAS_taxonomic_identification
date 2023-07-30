#!/bin/env/bash

# This scripts contains code to check the quality of reads, trim adapters and low quality reads
# $1= fastq1; $2 = fastq2; $3 = pair_id; $4 = adapter sequence;  

TrimmomaticPE -trimlog trimlogfile $1 $2 \
    "${3}_1P.fastq" "${3}_1U.fastq" "${3}_2P.fastq" "${3}_2U.fastq" \
    ILLUMINACLIP:$4:2:30:10:1:True SLIDINGWINDOW:30:20 MAXINFO:170:0.5 MINLEN:150