

// 1. Trimming with trimmomatic
// Inputs: reads, adapters; output:trimmed paired and unpaired reads
process trim_fastq {    

    input:
    val(trim_sample_script)
    tuple val(trim_fastq_reads_pair_id), path(trim_fastq_reads)
    val(adapter) 
   

    output:
    path "${trim_fastq_reads_pair_id}*.fastq", emit: trim_out

    script:
    
    """
    "./${trim_sample_script}" ${trim_fastq_reads[0]} ${trim_fastq_reads[1]} ${trim_fastq_reads_pair_id} ${adapter}
    """
}


// 2. Kraken analysis
// Input: reads (raw or trimmed); outputs: Classified and unclassified reads, kraken reports and outputs


process kraken_contamination {
    
    cpus = 8
    memory = 30.GB

    errorStrategy 'ignore'


    input:
    val (kraken_script)
    tuple val(kc_reads_pair_id), path(kc_reads)
    val (kraken_database)
    

    output:
        path "classified_*.fastq", emit: classified_fastq   classified_sample_2b_1.fastq
        path "unclassified_*.fastq", emit: unclassified_fastq 
        path "*_kraken_*", emit: kraken_reports
        
        
    script:
        """
        "./${kraken_script}" ${kc_reads[0]} ${kc_reads[1]} ${kc_reads_pair_id} ${kraken_database}

        """

}


// 3. 16S Ribosomal RNA Analysis
// Input: Fasta contig; Output: Blast output
process find_16s_hit {
    errorStrategy 'ignore'


    input:
        val find_16S_hits_script
        path fasta_contigs 


    output: 
        path "*_16s_output.tsv", emit: rRNA_16S_result


    script:
        """

        "./${find_16S_hits_script}" ${fasta_contigs}

        """


}

// 4. De novo genome assembly
// With spades
// Inputs: raw or trimmed reads

process spades_genome_assembly {
    errorStrategy 'ignore'
    publishDir 
    
    input:
        val (spades_script)
        tuple val(ga_reads_pair_id), path(gc_reads)

    output:
        path "*_*_contig.fasta"   

    script:
        """
        "./${spades_script}" ${gc_reads[0]} ${gc_reads[1]} ${ga_reads_pair_id}

        """

}

// With flye/medaka
process flye_genome_assembly {

    input:
        val (flye_script)
        path (ont_minion_read)

    output:
        path "consensus.fasta"

    script:
        """

        "./${flye_script}" ${ont_minion_read}

        """
}

// 5. Map reads to Reference
// Emit sam
process emit_sam_bwa {
    
    input:
        val (emit_sam_bwa_script)
        tuple val(emit_sam_reads_pair_id), path(emit_sam_reads)
        val (fastafile)
        


    output:
        path "${emit_sam_reads_pair_id}_sam", emit: reads_sam

    script:
        """

        bash ${emit_sam_bwa_script} ${emit_sam_reads[0]} ${emit_sam_reads[1]} ${fastafile} ${emit_sam_reads_pair_id}

        """
}

process emit_sam_minimap {
    
    input:
        val (emit_sam_mm_script)
        path (emit_sam_reads)
        val (fastafile)
        


    output:
        path "*_sam", emit: reads_sam

    script:
        """

        "./${emit_sam_mm_script}" ${emit_sam_reads} ${fastafile} 

        """
}

// Emit bam
process coordsort_sam {


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
    
    input:
        val (bamtofastq_script)
        path bamfile
        val (platform)

    output:
        path "*fastq*"

    script:
        """

        "./${bamtofastq_script}"" ${bamfile} ${platform}

        """

}

