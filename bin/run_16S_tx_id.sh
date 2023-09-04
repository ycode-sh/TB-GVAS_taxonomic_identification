#!/bin/bash

# input contig = $1;  $2 = txdb path string


# Extract filename from contig extension)
filename=$(basename $1 _contig.fasta)

# Predict 16S rRNA sequence form contig 
barrnap --outseq "${filename}.fa" < $1 > "${filename}.gff"

# Extract contigs with header "16S_rRNA" using grep command
grep -A 1 "16S_rRNA" "${filename}.fa" | grep -v "^--$" > "16s_${filename}.fa"

# Introduce taxdb files into the execution context (Not figured it out yet)
#script_abs_path="$(dirname "$0")"
#taxdb.btd = $2; taxdb.bti
cp ${2}taxdb* "$PWD"



# Run blast+ analysis
db_files_rel_path="/home/dfgmrtc/Workflows/wf-taxo_id_rd/wf-module_templates/module_scripts/resources/r16SRNA_db/16S_ribosomal_RNA.dust_window_masked"
db_files_abs_path="$db_files_rel_path"
blastn -query "16s_${filename}.fa"  -task blastn -db $db_files_abs_path \
    -db_hard_mask 30 -dust yes   -outfmt "7 qseqid sseqid pident length mismatch qstart qend sstart send evalue qcovs bitscore sscinames scomnames staxids stitle" \
    -out "${filename}_16s_output.tsv"