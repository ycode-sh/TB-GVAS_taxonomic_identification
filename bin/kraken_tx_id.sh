#!/bin/env/bash

# This scripts takes a pair of reads (if PE reads) - or single end reads  - 
# and checks for taxonomic contamination against a database

# $1 = read1, $2 = read2/db if ont or se, $3 = pair_id/read-type if ont or se, $4 = kraken_database if illumina to be supplied by user on the terminal, $5 = input read type if illumina
if [[ $3 == "se_illumina_reads" || $3 == "minion_ont_reads" ]]; then
    if [[ $3 == "se_illumina_reads" ]]: then
        sample_id=$(basename $1 _P.fastq)
    elif [[ $3 == "minion_ont_reads" ]]; then
        sample_id=$(basename $1 .fastq)
    fi
    kraken2 --db $2 --classified-out classified_$sample_id.fastq --unclassified-out unclassified_$sample_id.fastq \
    --output "${sample_id}_kraken_output" --report "${sample_id}_kraken_report" ${1} --use-names
else
    kraken2 --db  $4  --classified-out classified_$3#.fastq --unclassified-out unclassified_$3#.fastq \
    --output ${3}_kraken_output --report ${3}_kraken_report --paired ${1} ${2} --use-names
fi

#classified_*.fastq