#!/bin/bash
# $1=sam file

samtools sort $1 | samtools view -b > "${1}.bam"
