#!/bin/env/bash

# $1= fastq1, $2= fastq2, $3= pair_id

spades.py --careful -1 $1 -2 $2 -o ${3}_spades_result
cd ${3}_spades_result
mv contigs.fasta ${3}_contig.fasta

mv ${3}_contig.fasta ../
