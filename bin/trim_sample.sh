#!/bin/env/bash

# This scripts contains code to check the quality of reads, trim adapters and low quality reads
# $1= fastq1; $2 = fastq2; $3 = pair_id; $4 = adapter sequence;  $5 = input_read_type

# $1 = seq, $2 = adapter sequence for se reads, $3 = input_read_type

if [[ $3 == "se_illumina_reads" ]]; then
    file_name=$(basename $1 .fastq)
    TrimmomaticSE  -trimlog trimLogFile $1 "${file_name}_P.fastq" ILLUMINACLIP:$2:2:30:10:1:True SLIDINGWINDOW:30:20 MAXINFO:170:0.5 MINLEN:150

elif [[ $5 == "pe_illumina_reads" || $5 == "pe_min" || $5 == "pe_only" ]]; then
    TrimmomaticPE -trimlog trimlogfile $1 $2 \
        "${3}_1P.fastq" "${3}_1U.fastq" "${3}_2P.fastq" "${3}_2U.fastq" \
        ILLUMINACLIP:$4:2:30:10:1:True SLIDINGWINDOW:30:20 MAXINFO:170:0.5 MINLEN:150
fi
