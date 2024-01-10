#!/bin/bash

# Takes a coordinate-sorted bam file as input, marks duplicates and calls variants from marked barm files. The marked bam (bam file passed to Haplotype caller), fasta file,
# must be indexed with samtools index and faidx respectively. Reference dictionary for Haplotype caller must be generated and supplied in real time

# Parameters
# $1 = Coord_sort bam, $2 = fasta file, $3 = min-BQ, $4 = ploidy, $5 = gvcf_choice ("gvcf")

sample_id=$(basename $1 _sam.bam)

# Mark duplicates
#gatk MarkDuplicates \
#    --INPUT "$1" \
#    --METRICS_FILE marked_dup_metrics.txt \
#    --OUTPUT "${sample_id}_marked.bam"

# Create fasta and bam indexes
samtools faidx $2
#samtools index "${sample_id}_marked.bam"
samtools index $1


# Create reference dictionary
gatk CreateSequenceDictionary -R $2 -O h37rv.dict

# Call variants
if [[ $5 == "gvcf" || $5 == "cgvcf" ]]; then
    gatk HaplotypeCaller \
        --reference $2 \
        --input  $1 \
        --min-base-quality-score $3 \
        --sample-ploidy $4 \
        --output "${sample_id}_gatk.vcf" \
        --emit-ref-confidence GVCF
else
    gatk HaplotypeCaller \
        --reference $2 \
        --input "${sample_id}_marked.bam" \
        --min-base-quality-score $3 \
        --sample-ploidy $4 \
        --output "${sample_id}_gatk.vcf" 
fi
