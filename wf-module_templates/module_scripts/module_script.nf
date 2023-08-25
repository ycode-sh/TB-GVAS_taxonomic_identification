

// 1. Trimming with trimmomatic
// Inputs: reads, adapters; output:trimmed paired and unpaired reads
process trim_fastq {    
    publishDir "${params.out_dir}/Trimmed_reads"
    
    input:
    val(trim_sample_script)
    val(trim_fastq_se_reads)
    tuple val(trim_fastq_reads_pair_id), path(trim_fastq_reads)
    val(adapter)
    val(input_read_type) 
   

    output:
    path "*.fastq", emit: trim_out

    script:
    
    """
    if [ ! -z ${trim_fastq_se_reads} ]; then
        "./${trim_sample_script}" ${trim_fastq_se_reads} ${adapter} ${input_read_type}
    elif [ ! -z ${trim_fastq_reads_pair_id} ]; then
        "./${trim_sample_script}" ${trim_fastq_reads[0]} ${trim_fastq_reads[1]} ${trim_fastq_reads_pair_id} ${adapter} ${input_read_type}
    fi

    """
}


// 2. Kraken analysis
// Input: reads (raw or trimmed); outputs: Classified and unclassified reads, kraken reports and outputs


process kraken_contamination {
    publishDir "${params.out_dir}/Kraken_outputs"
    cpus = 8
    memory = 30.GB

    errorStrategy 'ignore'


    input:
    val (kraken_script)
    path kc_se_reads
    tuple val(kc_reads_pair_id), path(kc_reads)
    val (kraken_database_path)
    val (input_read_type)
    
    output:
        path "classified_*.fastq", emit: classified_fastq   classified_sample_2b_1.fastq
        path "unclassified_*.fastq", emit: unclassified_fastq 
        path "*_kraken_*", emit: kraken_reports
        
        
    script:
        """
        if [ ! -z ${kc_se_reads} ]; then
            "./${kraken_script}" ${kc_reads} ${kraken_database_path} ${input_read_type}
        elif [ ! -z ${kc_reads_pair_id} ]; then
            "./${kraken_script}" ${kc_reads[0]} ${kc_reads[1]} ${kc_reads_pair_id} ${kraken_database_path} ${input_read_type}
        fi
        """

}


// 3. 16S Ribosomal RNA Analysis
// Input: Fasta contig; Output: Blast output
process find_16s_hit {
    publishDir "${params.out_dir}/16S_outputs"
    errorStrategy 'ignore'


    input:
        val find_16S_hits_script
        path fasta_contigs
        val (txdb_path_str) 


    output: 
        path "*_16s_output.tsv", emit: rRNA_16S_result


    script:
        """

        "./${find_16S_hits_script}" ${fasta_contigs} ${txdb_path_str}

        """


}

// 4. De novo genome assembly
// With flye or spades
// Inputs: raw or trimmed reads

process genome_assembly {
    publishDir "${params.out_dir}/Contigs"
    errorStrategy 'ignore'
    
    input:
        val (genome_assembly_script)
        path ga_se_reads
        tuple val(ga_reads_pair_id), path(ga_reads)
        val (input_read_type)

    output:
        path "*_contig.fasta"   

    script:
        """
        if [ ! -z ${ga_se_reads} ]; then
            "./${genome_assembly_script}"  ${ga_se_reads} ${input_read_type}
        elif [ ! -z ${ga_reads_pair_id} ]; then
            "./${genome_assembly_script}" ${ga_reads[0]} ${ga_reads[1]} ${ga_reads_pair_id} ${input_read_type}
        fi
        """

}


// 5. Map reads to Reference
// Emit sam with either bwa or minimap
process emit_sam {
    publishDir "${params.out_dir}/Sam_files"
    input:
        val (emit_sam_script)
        path emit_sam_se_reads
        tuple val(emit_sam_reads_pair_id), path(emit_sam_reads)
        val (fastafile)
        val (input_read_type)
        


    output:
        path "*_sam", emit: reads_sam

    script:
        """
        if [ ! -z ${emit_sam_se_reads} ]; then
        "./${emit_sam_script}" ${emit_sam_se_reads} ${fastafile} ${input_read_type}
        elif [ ! -z ${emit_sam_reads_pair_id} ]; then
        "./${emit_sam_script}" ${emit_sam_reads[0]} ${emit_sam_reads[1]} ${fastafile} ${emit_sam_reads_pair_id} ${input_read_type}
        fi
        """
}



// Emit bam
process coordsort_sam {
        publishDir "${params.out_dir}/Bam_files"

        input:
            val (coordsort_sam_script)
            path reads_sam



        output:
            path "*_sam.bam", emit: reads_bam
            

        script:
            """

            "./${coordsort_sam_script}" ${reads_sam}

            """
}

// Reads decontamination

process bamtofastq {
    publishDir "${params.out_dir}/Decontaminated_reads"
    input:
        val (bamtofastq_script)
        path bamfile
        val (platform)

    output:
        path "*fastq*"

    script:
        """

        "./${bamtofastq_script}" ${bamfile} ${platform}

        """

}

// File parsing processes

process parse_16S {
    publishDir "${params.out_dir}/Parsed_16S_outputs"
    input:
        val (parse_16S_script)
        path output_16S_files
        val (perc_cov)
        val (perc_id) 

    output:
        path "16S_samples_output.tsv", emit: parsed_16S_out

    script:
        """

        "./${parse_16S_script}" ${output_16S_files} ${perc_cov} ${perc_id}

        """

}

process parse_kraken {
    publishDir "${params.out_dir}/Parsed_kraken_outputs"
    input:
        val (parse_kraken_script)
        path output_kraken_files


    output:
        path "parsed_kraken_report.tsv", emit: parsed_kraken_report

    script:
        """

        "./${parse_kraken_script}" ${output_kraken_files}

        """


}