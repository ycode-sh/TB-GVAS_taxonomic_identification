#!/bin/bash

# Takes flattened gvcf files and combine into s single joint_genotyped.vcf 
# Parameters 
# $@ = flattened gvcfs

# Initialize an empty array
vcf_array=""
# Iterate through the command line options and grow the initialized array
for file in "$@"; do
    if [ -s "$file" ]; then
        if [[ $file =~ [[:alnum:]]*\.fas?t?a?$ ]]; then # Match and extract fastafile
            ref=$file
        else
            vcf_array="$vcf_array --variant $file"
        fi
    fi
done

# Index reference file, and run CombineGVCFs to combine the vcf files in vcf_array
samtools index $ref

gatk CombineGVCFs \
   --reference $ref \
   $vcf_array \
   --output combined_gvcf.vcf

# Create sequence dictionary and genotype combined_gvcf

if [[ $ref =~ [[:alnum:]]*\.fasta$ ]]; then
    ref_base=$(basename $ref .fasta)
elif [[ $ref =~ [[:alnum:]]*\.fa$ ]]; then
    ref_base=$(basename $ref .fa)
fi


gatk CreateSequenceDictionary -R $ref -O ${ref_base}.dict

gatk --java-options "-Xmx4g" GenotypeGVCFs \
    --reference $ref \
    --variant combined_gvcf.vcf \
    --output gatk_joint_genotyped.vcf

