#!/bin/bash

# This script takes genotyped vcf files and annotates the variants with snpEff
# $1 = vcf file, $2 = path to snpEff.jar $3 = path to snpEff.config

filename=$(basename $1 .vcf)
java -Xmx8g -jar $2  -c $3 \
    Mycobacterium_tuberculosis_h37rv $1 > "${filename}_ann.vcf"

    
mv snpEff_summary.html ${filename}.snpEff_summary.html
mv snpEff_genes.txt ${filename}.snpEff_genes.txt
