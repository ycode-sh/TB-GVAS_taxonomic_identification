#!/bin/bash
# fastq1=$1; fastq2=$2; fastafile=$3; sample_ID =$4; $5= input_read_type for pe
# $1 = fastq_read, $2 = fastafile, $3 = input_read_type  for se 


if [[ $3 == "se_illumina_reads" ]]; then
    filename=$(basename $1 .fastq)
    bwa index $2
    bwa mem -R "@RG\tID:run_$filename\tPL:ILLUMINA\tLB:paired_end\tSM:$filename" $2 $1 -o "${filename}_sam"

elif [[ $3 == "minion_ont_reads" ]]; then # indexes fasta reference and aligns the reads to the reference supplied using minimap
    sample_id=$(basename $1 .fastq)
    minimap2 -d $2 h37rv.mmi
    minimap2 -a h37rv.mmi $1 -R @RG\\tID:$sample_id\\tSM:$sample_id -o "${sample_id}_sam"
elif [[ $5 == "pe_illumina_reads" ]]; then
    bwa index $3
    bwa mem -R "@RG\tID:run_$4\tPL:ILLUMINA\tLB:paired_end\tSM:$4" $3 $1 $2 -o "${4}_sam"
fi
