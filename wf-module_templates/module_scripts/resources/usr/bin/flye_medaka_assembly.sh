#!/bin/bash

#This script takes $1, ONT fastq file and uses flye assmbler to generate a draft assembly. The draft assembly is passed to medaka 
# which then polishes it to generate a consensus fasta file.
# $1=fastq file

flye --nano-raw $1 --out-dir nano_output
cd nano_output
mv assembly.fasta ../

medaka_consensus -i $1 -d assembly.fasta  -o medaka_output 
cd medaka_output
mv consensus.fasta ../
