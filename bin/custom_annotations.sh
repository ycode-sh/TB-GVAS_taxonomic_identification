#!/bin/bash

# This script performs a sequential custom annotation on vcf files using annotion files. 
# It annotates the variants in a vcf file that intersects with lineage, amr_region, gene_annotation, and variable_regions information
# $1= gzipped and tabix-indexed gene_annotation_file, $2 = gene_annotation_header, $3 = lineage_annotation_file, $4 = lineage_annotation_header, 
# $5 = amr_region_file, $6 = amr_region_header, $7 = variable_region_annotation_file, $8 = variable_region_annotation_header $9 = vcf_file

#sample_2_bt_ann.vcf

file_name=$(basename $9 _ann.vcf)
bcftools annotate --annotations  ${3} --header-lines $4  --columns chrom,from,to,lineage,exp_allele,lsp,spoligotype,lin_ref  --merge-logic lineage:append --merge-logic exp_allele:append --merge-logic lsp:append --merge-logic spoligotype:append --merge-logic lin_ref:append $9 > ${file_name}_la.vcf
bcftools annotate --annotations  ${5} --header-lines $6  --columns CHROM,from,to,antibiotics_gene,related_antibiotics ${file_name}_la.vcf > ${file_name}_ar.vcf
bcftools annotate --annotations  ${7} --header-lines $8  --columns CHROM,from,to,vr_locus_tag,comment ${file_name}_ar.vcf > ${file_name}_vr.vcf
bcftools annotate --annotations  ${1} --header-lines $2  --columns CHROM,from,to,gene_ann_region,gene_annotation --merge-logic gene_annotation:append --merge-logic gene_ann_region:append ${file_name}_vr.vcf > ${file_name}_c_ann.vcf


rm ${file_name}_la.vcf ${file_name}_ar.vcf ${file_name}_vr.vcf 
#bcftools annotate --annotations $1 --header-lines $2  --columns CHROM,from,to,gene_ann_region,gene_annotation $3 > ${file_name}_ga.vcf
