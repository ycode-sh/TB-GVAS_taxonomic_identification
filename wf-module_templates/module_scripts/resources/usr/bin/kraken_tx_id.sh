#!/bin/env/bash

# This scripts takes a pair of filtered or unfiltered reads (if PE reads) - or single end reads - 
# and checks for taxonomic contamination against a database

# $1 = read1, $2 = read2, $3 = pair_id, $4 = kraken_database to be supplied by user on the terminal

kraken2 --db  $4  --classified-out classified_$3#.fq --unclassified-out unclassified_$3#.fq \
    --output ${3}_kraken_output --report ${3}_kraken_report --paired ${1} ${2} --use-names