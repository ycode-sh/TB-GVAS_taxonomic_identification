#!/bin/bash


# Filter only called snps. Don't recalculate INFO fields and output alleles present in the ALT field

# $1 = referenece_fasta, $2 = vcf_file

filename=$(basename $2 _subtract_repeat.vcf) 

bcftools view --exclude-types indels,mnps,bnd,other --no-update --trim-alt-alleles --exclude 'GT="0/0"' --exclude-uncalled $2 --output "${filename}_snps_temp.vcf"

bcftools view --include "FORMAT/DP >= 10" "${filename}_snps_temp.vcf" --output "${filename}_snps.vcf"

rm "${filename}_snps_temp.vcf"

bgzip "${filename}_snps.vcf"

bcftools index "${filename}_snps.vcf.gz"

bcftools consensus --fasta-ref $1 \
    --mark-del D --mark-ins lc --mark-snv lc "${filename}_snps.vcf.gz" --output ${filename}.vcf




# Read the input FASTA file
while IFS= read -r line
do
    if [[ $line =~ ^\> ]]; then
        # Replace the header line with the new header
        echo ">$filename" >> "${filename}_consensus.vcf"
    else
        # Write the sequence lines as they are
        echo "$line" >> "${filename}_consensus.vcf"
    fi
done < "${filename}.vcf"

echo "Header renamed successfully. Output file: ${filename}_consensus.vcf"




