#!/bin/bash

# Marks duplicates in  alignement file. It takes aligned sam file as input and converts to marked bam files. Then uses the genotype likelihood results from bcfpileup to 
# call variants for each genomic position. Indexed fasta file must be available in execution context. 
# The file may be compressed). Also, the --read-groups FILE option may be included. 

# Parameters
# $1 = aligned_sam, $2 = fasta file, $3 = min-BQ, $4 = min-MQ, $5 = ploidy, $6 = tandem-qual, $7 = indel-bias
sample_id=$(basename $1 _sam)   
 
samtools view -b $1 | samtools sort -n -o namesort_${sample_id}.bam
samtools fixmate -m namesort_${sample_id}.bam fixmate_${sample_id}.bam
samtools sort -o positionsort_${sample_id}.bam fixmate_${sample_id}.bam
samtools markdup positionsort_${sample_id}.bam ${sample_id}_marked.bam

bcftools mpileup --full-BAQ \
    --fasta-ref $2 --min-BQ $3 --min-MQ $4 \
   --skip-any-set UNMAP,SECONDARY,QCFAIL,DUP  --tandem-qual $6 --indel-bias $7 \
   --output-type u ${sample_id}_marked.bam | bcftools call --ploidy $5 \
   --multiallelic-caller --variants-only --output ${sample_id}_bt.vcf






