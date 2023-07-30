#!/bin/env/bash

# This scripts takes a pair of reads (if PE reads) - or single end reads  - 
# and checks for taxonomic contamination against a database

# $1 = read1, $2 = read2/db if ont or se, $3 = pair_id/read-type if ont or se, $4 = kraken_database if illumina to be supplied by user on the terminal
if [[ $3 == "se" ]]; then
    sample_id=$(basename $1 .fastq)
    kraken2 -db $2 --classified-out classified_$sample_id#.fq --unclassified-out unclassified_$sample_id#.fq \
    --output "${sample_id}_kraken_output" --report "${sample_id}_kraken_report" ${1} --use-names
else
    kraken2 --db  $4  --classified-out classified_$3#.fq --unclassified-out unclassified_$3#.fq \
    --output ${3}_kraken_output --report ${3}_kraken_report --paired ${1} ${2} --use-names

