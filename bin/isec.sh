#!/bin/bash

# Parameters
# $1 = bcf.gvcf, $2 = gatk.gvcf, 

bgzip $1
bgzip $2

tabix ${1}.gz
tabix ${2}.gz

bcftools isec -p isec_dir ${1}.gz ${2}.gz
cd isec_dir
mv 0000.vcf mino_only_joint.vcf
mv 0001.vcf gatk_only_joint.vcf
mv 0002.vcf minos_gatk_joint.vcf
mv 0003.vcf gatk_minos_joint.vcf
mv *vcf ../
mv README.txt  sites.txt ../