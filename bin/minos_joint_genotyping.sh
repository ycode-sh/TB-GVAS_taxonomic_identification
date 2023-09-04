#!/bin/bash
# parameters
# $1 = minos_config, $2= minos_profile, $3 = regenotype_nf_script, $4 = fastafile, $5 = manifest

for item in "$@"; do
    if [[ $item =~ sample_[0-9]_sam\.bam$ ]]; then
        samtools index $item
    fi
done

mv ${5} manifest.tsv
nextflow run -c ${1} -profile ${2} ${3} --ref_fasta $4 --manifest manifest.tsv --outdir "OUT.joint_geno"
#nextflow run  -c ~/minos/nextflow/regenotype.config  -profile medium  ~/minos/nextflow/regenotype.nf   --ref_fasta ~/Dr_Francis/Scripts/Analysis_data/reference_data/h37rv.fasta   --manifest manifest.tsv   --outdir OUT.joint_geno

cd OUT.joint_geno
mv merged.vcf minos_joint_genotyped.vcf
mv minos_joint_genotyped.vcf ../

#sample_1_sam.bam