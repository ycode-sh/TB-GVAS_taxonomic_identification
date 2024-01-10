#!/bin/bash

# This comand takes input vcf files, fastq files, and fasta reference and output result in a specified outdir
# input_read_type = $1, sample_id = $2, reads_1 = $3, reads_2 = $4, vcf_1 = $5, vcf_2 = $6, fastafile = $7, min_frs = $8, min_gcp = $9
# input_read_type = $1, read = $2, fastafile = $3, vcf_1 = $4, vcf_2 = $5, min_frs = $6, min_gcp = 7


#vcf_array=()
#for item in "$@"; do
#    if [[ $item =~ ^[[:alnum:]]*\.vcf$ ]]; then
#        vcf_array+=("$item ")
#    fi
#done
if [[ $1 == "pe_illumina_reads" ]]; then
    minos adjudicate --filter_min_frs $8 --filter_min_gcp $9 --reads $3 --reads $4 ${2}_minos_adj_run $7 $5 $6
    cd "${2}_minos_adj_run"
    mv final.vcf "${2}_minos.vcf"
    mv "${2}_minos.vcf" ../
elif [[ $1 == "se_illumina_reads" ]]; then
    sample_name=$(basename $2 _P.fastq)
    minos adjudicate --filter_min_frs $6 --filter_min_gcp $7 --reads $2  "${sample_name}_minos_adj_run" $3 $4 $5
    cd "${sample_name}_minos_adj_run"
    mv final.vcf "${sample_name}_minos.vcf"
    mv "${sample_name}_minos.vcf" ../    

fi
