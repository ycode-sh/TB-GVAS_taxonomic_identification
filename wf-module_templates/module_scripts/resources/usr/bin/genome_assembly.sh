#!/bin/bash

# $1= fastq1, $2= fastq2, $3= pair_id , $4 = input_read_type for pe 
# $1= fastq, $2 = input_read_type for se

if [[  $2 == "se_illumina_reads" ]]; then
    filename=$(basename $1 P.fastq)
    spades.py --careful -s $1 -o ${filename}_spades_result
    cd ${filename}_spades_result
    mv contigs.fasta ${filename}_contig.fasta
elif [[ $2 == "minion_ont_reads" ]]; then # uses flye assmbler to generate a draft assembly. The draft assembly is passed to medaka which then polishes it to generate a consensus fasta file.
    filename=$(basename $1 .fastq)
    flye --nano-raw $1 --out-dir nano_output
    cd nano_output
    mv assembly.fasta ../

    medaka_consensus -i $1 -d assembly.fasta  -o medaka_output 
    cd medaka_output
    mv consensus.fasta ${filename}_contig.fasta
    mv ${filename}_contigs.fasta ../

elif [[ $4 == "pe_illumina_reads" ]]; then
    mv ${filename}_contig.fasta ../
    spades.py --careful -1 $1 -2 $2 -o ${3}_spades_result
    cd ${3}_spades_result
    mv contigs.fasta ${3}_contig.fasta

    mv ${3}_contig.fasta ../
fi



