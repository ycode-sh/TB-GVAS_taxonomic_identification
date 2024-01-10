#1/bin/bash
# $1 = vcf_file, $2 = amr_genomic_interval, $3 = lineage_snps.bed, $4 = assay_type:  clin_assay
# $1 = vcf_file, $2 = variable_region_interval, $3 = assay_type:   epi_assay

if [[ $3 == "epi_assays" || $3 == "all_assays" ]]; then
    filename=$(basename "$1" .vcf)
    bcftools view --header-only $1 > "${filename}_subtract_repeat.vcf"
    bedtools subtract -a "${1}" -b $2>> "${filename}_subtract_repeat.vcf"
elif [[ $4 == "clin_assays" || $4 == "all_assays" ]]; then
    filename=$(basename "$1" _ann.vcf)
    bcftools view --header-only "$1" > "${filename}_intersect.vcf"
    bedtools intersect -a "$1" -b $2  $3 -wa -wb -loj -names amr_regions lineage_snps >> "${filename}_intersect.vcf"
fi
