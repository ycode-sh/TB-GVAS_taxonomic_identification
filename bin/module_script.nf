// 1. Trimming with trimmomatic


process trim_fastq_se {    
    //publishDir "${params.out_dir_int}/Trimmed_reads", mode: 'copy', overwrite: true
    
    input:
    val(trim_sample_script)
    path(se_reads)
    val(adapter)
    val(input_read_type) 
   

    output:
    path "sample*.fastq", emit: trimmed_reads
    

    script:
    
    """
    
    bash ${trim_sample_script} ${se_reads} ${adapter} ${input_read_type}
   

    """
}

process trim_fastq_pe {    
    //publishDir "${params.out_dir_int}/Trimmed_reads", mode: 'copy', overwrite: true
    
    input:
    val(trim_sample_script)
    tuple val(reads_pair_id), path(pe_reads)
    val(adapter)
    val(input_read_type) 
   

    output:
    path "sample*.fastq", emit: trimmed_reads
    

    script:
    
    """
    bash ${trim_sample_script} ${pe_reads[0]} ${pe_reads[1]} ${reads_pair_id} ${adapter} ${input_read_type}
    
    """
}

workflow trim_fastq {
    take:
        trim_sample_script_wf
        pe_reads_wf
        adapter_seq
        input_read_type
        se_reads_wf
    main:
        if (params.cohort_joint_genotype == "ncgvcf"){
            if (params.in_data_type != "minion_ont_reads"){
                if (params.in_data_type == "pe_illumina_reads"){
                    trim_fastq_pe(trim_sample_script_wf, pe_reads_wf, adapter_seq, input_read_type)
                    trimmed_reads = trim_fastq_pe.out.trimmed_reads | flatten | filter { it =~/P.fastq/ } | map { [it.name - ~/_[12]P.fastq/, it] } | groupTuple
                } else if (params.in_data_type == "se_illumina_reads"){
                    trim_fastq_se(trim_sample_script_wf, se_reads_wf, adapter_seq, input_read_type)
                    trimmed_reads = trim_fastq_se.out
                }
            } else {
                trimmed_reads = se_reads_wf
            }
        } else if(params.cohort_joint_genotype == "cgvcf"){
            if (params.in_data_type == "pe_min"){
                trim_fastq_pe(trim_sample_script_wf, pe_reads_wf, adapter_seq, input_read_type)
                trimmed_reads = trim_fastq_pe.out.trimmed_reads | flatten | filter { it =~/P.fastq/ } | map { [it.name - ~/_[12]P.fastq/, it] } | groupTuple
                trimmed_reads_2 = se_reads_wf
            }
            else if(params.in_data_type == "pe_only"){
                trim_fastq_pe(trim_sample_script_wf, pe_reads_wf, adapter_seq, input_read_type)
                trimmed_reads = trim_fastq_pe.out.trimmed_reads | flatten | filter { it =~/P.fastq/ } | map { [it.name - ~/_[12]P.fastq/, it] } | groupTuple
            }      
            
        }
    emit:
        trimmed_reads
        //trimmed_reads_2
}


// 2. Kraken analysis
// Input: reads (raw or trimmed); outputs: Classified and unclassified reads, kraken reports and outputs

process kraken_contamination_se {
    publishDir "${params.out_dir_int}/Kraken_outputs", mode: 'copy', overwrite: true
    cpus = 8
    memory = 30.GB

    errorStrategy 'ignore'


    input:
    val (kraken_script)
    path(se_reads)
    val (kraken_database_path)
    val (input_read_type)
    
    output:
        path "classified_*.fastq", emit: classified_fastq
        path "unclassified_*.fastq", emit: unclassified_fastq 
        path "*_kraken_*", emit: kraken_reports
        
        
    script:
        """
        bash ${kraken_script} ${se_reads} ${kraken_database_path} ${input_read_type}
        
        """

}

process kraken_contamination_pe {
    publishDir "${params.out_dir_int}/Kraken_outputs", mode: 'copy', overwrite: true
    cpus = 8
    memory = 30.GB

    errorStrategy 'ignore'


    input:
    val (kraken_script)
    tuple val(reads_pair_id), path(pe_reads)
    val (kraken_database_path)
    val (input_read_type)
    
    output:
        path "classified_*.fastq", emit: classified_fastq
        path "unclassified_*.fastq", emit: unclassified_fastq 
        path "*_kraken_*", emit: kraken_reports
        
        
    script:
        """
        bash ${kraken_script} ${pe_reads[0]} ${pe_reads[1]} ${reads_pair_id} ${kraken_database_path} ${input_read_type}
        """

}

workflow kraken_contamination {
    take:
        kraken_script_wf
        reads_wf
        kraken_db_path_wf
        input_read_type


    main:
        if (params.in_data_type == "pe_illumina_reads"){
            kraken_contamination_pe(kraken_script_wf, reads_wf, kraken_db_path_wf, input_read_type)
            kraken_outputs = kraken_contamination_pe.out
        } else if (params.in_data_type == "se_illumina_reads" | params.in_data_type == "minion_ont_reads"){
            kraken_contamination_se(kraken_script_wf, reads_wf, kraken_db_path_wf, input_read_type)
            kraken_outputs = kraken_contamination_se.out
        }

    emit:
        kraken_outputs.kraken_reports
}

// 3. 16S Ribosomal RNA Analysis
// Input: Fasta contig; Output: Blast output
process find_16s_hit {
    publishDir "${params.out_dir_int}/16S_outputs", mode: 'copy', overwrite: true
    errorStrategy 'ignore'


    input:
        val find_16S_hits_script
        path fasta_contigs
        val (txdb_path_str) 


    output: 
        path "*_16s_output.tsv", emit: rRNA_16S_result


    script:
        """

        bash ${find_16S_hits_script} ${fasta_contigs} ${txdb_path_str}

        """


}

// 4. De novo genome assembly
// With flye or spades
// Inputs: raw or trimmed reads

process genome_assembly_se {
    publishDir "${params.out_dir_int}/Contigs", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    cpus = 8
    memory = 30.GB
    
    input:
        val (genome_assembly_script)
        path(se_reads)
        val (input_read_type)

    output:
        path "*_contig.fasta"   

    script:
        """
        bash ${genome_assembly_script}  ${se_reads} ${input_read_type}
        """

}

process genome_assembly_pe {
    cpus = 8
    memory = 30.GB
    publishDir "${params.out_dir_int}/Contigs", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    input:
        val (genome_assembly_script)
        tuple val(reads_pair_id), path(pe_reads)
        val (input_read_type)

    output:
        path "*_contig.fasta"   

    script:
        """
        bash ${genome_assembly_script} ${pe_reads[0]} ${pe_reads[1]} ${reads_pair_id} ${input_read_type}
       
        """

}

workflow genome_assembly {
    take:
        genome_assembly_script_wf
        reads_wf
        input_read_type


    main:
        if (params.in_data_type == "pe_illumina_reads"){
            genome_assembly_pe(genome_assembly_script_wf, reads_wf, input_read_type)
            asssembly_out = genome_assembly_pe.out
        } else if (params.in_data_type == "se_illumina_reads" | params.in_data_type == "minion_ont_reads") {
            genome_assembly_se(genome_assembly_script_wf, reads_wf, input_read_type)
            asssembly_out = genome_assembly_se.out
        }

    emit:
        asssembly_out
}


// 5. Map reads to Reference
// Emit sam with either bwa or minimap
process emit_sam_se {
    //publishDir "${params.out_dir_int}/Sam_files", mode: 'copy', overwrite: true
    cpus = 8
    memory = 30.GB
    input:
        val (emit_sam_script)
        path(se_reads)
        val (fastafile)
        val (input_read_type)
        


    output:
        path "*_sam", emit: reads_sam

    script:
        """
        
        bash ${emit_sam_script} ${se_reads} ${fastafile} ${input_read_type}
        
        """
}

process emit_sam_pe {
    //publishDir "${params.out_dir_int}/Sam_files", mode: 'copy', overwrite: true
    cpus = 8
    memory = 30.GB
    input:
        val (emit_sam_script)
        tuple val(reads_pair_id), path(pe_reads)
        val (fastafile)
        val (input_read_type)
        


    output:
        path "*_sam", emit: reads_sam

    script:
        """
        bash ${emit_sam_script} ${pe_reads[0]} ${pe_reads[1]} ${fastafile} ${reads_pair_id} ${input_read_type}
        
        """
}

process emit_sam_pe_min {
    //publishDir "${params.out_dir_int}/Sam_files", mode: 'copy', overwrite: true
    cpus = 8
    memory = 30.GB
    input:
        val (emit_sam_script)
        tuple val(reads_pair_id), path(pe_reads)
        val (fastafile)
        val (input_read_type)
        path (se_reads)
        


    output:
        path "*_sam", emit: reads_sam

    script:
        """
        bash ${emit_sam_script} ${pe_reads[0]} ${pe_reads[1]} ${fastafile} ${reads_pair_id} ${input_read_type} ${se_reads}
        
        """
}

workflow emit_sam {
    take:
        emit_sam_script_wf
        trimmed_reads
        trimmed_reads_2
        fastafile
        input_read_type

    main:
        if (params.in_data_type == "se_illumina_reads" | params.in_data_type == "minion_ont_reads"){
            emit_sam_se(emit_sam_script_wf, trimmed_reads, fastafile, input_read_type)
            emitted_sam = emit_sam_se.out
        } else if (params.in_data_type == "pe_illumina_reads" | params.in_data_type == "pe_only"){
            emit_sam_pe(emit_sam_script_wf, trimmed_reads, fastafile, input_read_type)
            emitted_sam = emit_sam_pe.out
        } else if(params.in_data_type == "pe_min"){
            //emit_sam_pe_min(emit_sam_script_wf, trimmed_reads, fastafile, input_read_type, trimmed_reads_2)
            //emitted_sam = emit_sam_pe_min.out | collect | flatten
            emit_sam_pe(emit_sam_script_wf, trimmed_reads, fastafile, input_read_type)
            emit_sam_se(emit_sam_script_wf, trimmed_reads_2, fastafile, input_read_type)

            emitted_sam = emit_sam_pe.out.concat(emit_sam_se.out)
            
        }
    emit:
        emitted_sam

}
// Emit bam
process coordsort_sam {
        publishDir "${params.out_dir_int}/Bam_files", mode: 'copy', overwrite: false

        input:
            val (coordsort_sam_script)
            path reads_sam



        output:
            path "*_sam.bam", emit: reads_bam
            

        script:
            """

            bash ${coordsort_sam_script} ${reads_sam}

            """
}

// Reads decontamination

process bamtofastq {
    publishDir "${params.out_dir_int}/Decontaminated_reads", mode: 'copy', overwrite: true
    input:
        val (bamtofastq_script)
        path bamfile
        val (input_read_type)

    output:
        path "*fastq", emit: decontaminated_reads

    script:
        """

        bash ${bamtofastq_script} ${bamfile} ${input_read_type}

        """

}

workflow pbamtofastq {
    take:
        bamtofastq_script_wf
        sorted_bamfile_wf
        input_read_type

    main:
        bamtofastq(bamtofastq_script_wf, sorted_bamfile_wf, input_read_type)
        if (params.in_data_type == "se_illumina_reads" | params.in_data_type == "minion_ont_reads") {
            bamtofastq_out = bamtofastq.out
        } else if (params.in_data_type == "pe_illumina_reads") {
            bamtofastq_out = bamtofastq.out | flatten | map {[it.name - ~/_[12].fastq/, it] } | groupTuple
        }
      

    emit:
        bamtofastq_out
}

// File parsing processes

process parse_16S {
    publishDir "${params.out_dir_int}/Parsed_16S_outputs", mode: 'copy', overwrite: true
    input:
        val (parse_16S_script)
        path output_16S_files
        val (perc_cov)
        val (perc_id) 

    output:
        path "parsed_16S_samples_output.tsv", emit: parsed_16S_out

    script:
        """

        python ${parse_16S_script} ${output_16S_files} ${perc_cov} ${perc_id}

        """

}

process parse_kraken {
    publishDir "${params.out_dir_int}/Parsed_kraken_outputs", mode: 'copy', overwrite: true
    input:
        val (parse_kraken_script)
        path output_kraken_files


    output:
        path "parsed_kraken_report.tsv", emit: parsed_kraken_report

    script:
        """

        python ${parse_kraken_script} ${output_kraken_files}

        """


}


// Bcftools variant calling

process bcf_call {
    publishDir "${params.out_dir_int}/bcf_call", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
    val bcf_call_script
    path aligned_sam
    val fastafile
    val min_BQ
    val min_MQ
    val ploidy
    val tandem_qual
    val indel_bias

    output:
        path "*_bt.vcf", emit: bcf_call_out

    script:
        """

        bash ${bcf_call_script} ${aligned_sam} ${fastafile} ${min_BQ} ${min_MQ} ${ploidy} ${tandem_qual} ${indel_bias}

        """


}


//  GATK variant calling

// Mark duplicates from coord-sort bam files (bam index created in-situ/real time)


process gatk_call {
    publishDir "${params.out_dir_int}/gatk_call", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        val gatk_call_script
        path coord_sort_bam
        val fastafile
        val min_BQ
        val ploidy

    output:
        path "*_gatk.vcf", emit: gatk_out_vcf

    script:
        """

        bash ${gatk_call_script} ${coord_sort_bam} ${fastafile} ${min_BQ} ${ploidy}

        """
}

process gatk_call_gvcf {
    publishDir "${params.out_dir_int}/gvcf_call", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        val gatk_call_script
        path coord_sort_bam
        val fastafile
        val min_BQ
        val ploidy
        val gvcf_choice

    output:
        path "*_gatk.vcf", emit: gatk_out_vcf

    script:
        """

        bash ${gatk_call_script} ${coord_sort_bam} ${fastafile} ${min_BQ} ${ploidy} ${gvcf_choice}

        """
}

process gatk_joint_genotyping {
    publishDir "${params.out_dir_int}/gatk_joint_call", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        val gatk_joint_geno_script
        path flattened_gvcfs
        val fastafile

    output:
        path "gatk_joint_genotyped.vcf"
    
    script:

    """

    bash ${gatk_joint_geno_script} ${flattened_gvcfs} ${fastafile} 

    """

}



// 9. Minos Variant adjudication
// Note "all files" variables below must contain paired reads and variants with similar prefix and suffix, only differing by sample id (integer)
// Inputs: Reference File, reads pair, vcf files pair (The latter two are grouped into a keyed-tuple)

//read1= $1 read2= $2 fastafile = $3 vcf_1 = $4 vcf_2 = $5 samplename = $6
process adjudicate_var_pe {
    publishDir "${params.out_dir_int}/Minos_call", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    input:
    val var_adj_script
    val (input_read_type)
    tuple val(sample_id), path(reads), path(vcf_1, stageAs: 'ist_vcf'), path(vcf_2)
    //tuple val(sample_id), path(items)
    val reference  
    val (frs)
    val (gcp)
    
    output:
    path  "*_minos.vcf", emit: minos_call_out_vcf

    script:
    """
    bash ${var_adj_script} ${input_read_type} ${sample_id} ${reads[0]} ${reads[1]} ${vcf_1} ${vcf_1} ${reference} ${frs} ${gcp}


    """
}

process adjudicate_var_se {
    publishDir "${params.out_dir_int}/Minos_call", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

  
    input:
    val var_adj_script
    val (input_read_type)
    path (se_reads) 
    val reference 
    tuple val(sample_id), path(vcf_1, stageAs: 'ist_vcf'), path(vcf_2)
    val (frs)
    val (gcp)
    
    
    output:
    path  "*_minos.vcf", emit: minos_call_out_vcf

    script:
    """
    bash ${var_adj_script} ${input_read_type} ${se_reads} ${reference} ${vcf_1} ${vcf_2} ${frs} ${gcp}

    """
}

workflow adjudicate_var {
    take:
        bt_vcf
        gatk_vcf
        trimmed_reads
        var_adj_script_wf
        input_read_type
        fastafile
        min_frs
        min_gcp



    main:
        if (params.in_data_type != "minion_ont_reads"){
            bcf_call_tuple = bt_vcf | map { [it.name - ~/_bt.vcf/, it] }
            gatk_call_tuple = gatk_vcf | map { [it.name - ~/_gatk.vcf/, it] }
            joined_vcfs = bcf_call_tuple.join(gatk_call_tuple)
            if (params.in_data_type == "pe_illumina_reads") {
                reads_vcf_files = trimmed_reads.combine(joined_vcfs, by:0)
                adjudicate_var_pe(var_adj_script_wf, input_read_type, reads_vcf_files, fastafile, min_frs, min_gcp)
                adj_var = adjudicate_var_pe.out
            } else if (params.in_data_type == "se_illumina_reads"){
                adjudicate_var_se(var_adj_script_wf, input_read_type, trimmed_reads, fastafile, joined_vcfs, min_frs, min_gcp)
                adj_var = adjudicate_var_se.out
            }

        }
    emit:
        adj_var
}


process minos_manifest {
    
    input:
        tuple val(tuple_id), path(file_path1), path(file_path2)

    output:
        path "${tuple_id}_manifest.tsv"
    
    script:
        """

        echo -e "${tuple_id}\t${params.out_dir_int}/Bam_files/${file_path1}\t${params.out_dir_int}/bcf_call/${file_path2}" > "${tuple_id}_manifest.tsv"
        samtools index ${params.out_dir_int}/Bam_files/${file_path1}
        """
}

process combine_manifest {

    input:
        path flattened_manifest
    
    output:
        path "combined_manifest.tsv"

    script:
        """
        echo -e "name\treads\tvcf" > combined_manifest.tsv
        for file in ${flattened_manifest}; do
            cat \$file >> combined_manifest.tsv
        done

        """
}

workflow prepare_manifest {
    take:
        bam_file
        vcf_file

    main:
        bam_file_tuple = bam_file | map { [it.name - ~/_sam.bam/, it] }
        vcf_file_tuple = vcf_file | map { [it.name - ~/_bt.vcf/, it] }
        file_paths_tuple = bam_file_tuple.join(vcf_file_tuple, by:0)
        minos_manifest(file_paths_tuple)
        flattened_manifest = minos_manifest.out | flatten | collect 
        combine_manifest(flattened_manifest)
    emit:
        combine_manifest.out
}

process minos_joint_genotyping {
    publishDir "${params.out_dir_int}/Minos_joint_call", mode: 'copy', overwrite: true
    input:
        val minos_jg_script 
        path minos_config
        val (minos_profile)
        path regenotype_nf_script
        val (fastafile)
        val manifest
        path (flattened_bam)
        path (flattened_bt_vcfs)
        
        

    output:
        path "minos_joint_genotyped.vcf"

    script:
        """
        
        bash ${minos_jg_script} ${minos_config} ${minos_profile} ${regenotype_nf_script} ${fastafile} ${manifest} ${flattened_bam} ${flattened_bt_vcfs}
        

        """
}

process isec {
    publishDir "${params.out_dir_int}/isec", mode: 'copy', overwrite: true
    input:
        val (isec_script)
        path bcf_jvcf
        path gatk_jvcf


    output:
        path "*vcf", emit: isec_vcfs
        path "*.txt"  

    script:
    """

    bash ${isec_script} ${bcf_jvcf} ${gatk_jvcf}

    """
}


workflow joint_genotyping {
    take:
        minos_config_wf
        minos_profile_wf
        regenotype_nf_script_wf
        fastafile
        manifest_wf
        gatk_joint_geno_script_wf
        flattened_gvcfs
        isec_script_wf
        

    main:
        if (params.in_data_type != "minion_ont_reads"){
            minos_joint_genotyping(minos_config_wf, minos_profile_wf, regenotype_nf_script_wf, fastafile, manifest_wf)
            gatk_joint_genotyping(gatk_joint_geno_script, flattened_gvcfs, fastafile)   
            isec(isec_script_wf, minos_joint_genotyping.out, gatk_joint_genotyping.out)
            total_jvcf = isec.out
        } else if (params.in_data_type == "minion_ont_reads") {
            minos_joint_genotyping(minos_config_wf, minos_profile_wf, regenotype_nf_script_wf, fastafile, manifest_wf)
            total_jvcf = minos_joint_genotyping.out
        }
    
    emit: 
        total_jvcf
        
}

// Bedtools manipulation

process bedtools_epi_assay {
    publishDir "${params.out_dir_int}/sub_repeat", mode: 'copy', overwrite: true
    input:
        val bedtools_script
        path vcf_file
        val variable_region_interval
        val assay_type

    output:
        path "*repeat.vcf"

    script:
        """

        bash ${bedtools_script} ${vcf_file} ${variable_region_interval} ${assay_type}

        """

}


process bedtools_clin_assay {
    publishDir "${params.out_dir_int}/intersect", mode: 'copy', overwrite: true
    input:
        val bedtools_script
        path vcf_file
        val amr_genomic_interval
        val lineage_snps
        val assay_type


    output:
        path "*intersect.vcf"


    script:
        """

        bash ${bedtools_script} ${vcf_file} ${amr_genomic_interval} ${lineage_snps} ${assay_type}

        """
}

// Annotations

process snpEff_annotation {
    publishDir "${params.out_dir_int}/snpeff_annot", mode: 'copy', overwrite: true
    input:
        val (snpEff_annot_script) 
        path (vcf_file) 
        val (snpeff_jar_path)
        val (snpeff_config_path)


    output:
        path "*.vcf" , emit: snpeff_annotated_vcf
        path "*.html" 
        path "*.txt" 


    script:
    """

    bash ${snpEff_annot_script} ${vcf_file} ${snpeff_jar_path} ${snpeff_config_path}

    """
}

process custom_annotations {
    publishDir "${params.out_dir_int}/custom_annot", mode: 'copy', overwrite: true

    input:
        val custom_annotations_script
        val gene_annotation_file
        val gene_annotation_header
        val lineage_annotation_file
        val lineage_annotation_header
        val amr_region_file
        val amr_region_header
        val variable_region_annotation_file
        val variable_region_annotation_header
        path vcf_file

    output:
        path "*c_ann.vcf", emit: custom_annotated_vcf

    script:
        """
        bash ${custom_annotations_script} ${gene_annotation_file} ${gene_annotation_header} ${lineage_annotation_file} ${lineage_annotation_header} ${amr_region_file} ${amr_region_header} ${variable_region_annotation_file} ${variable_region_annotation_header} ${vcf_file}
        """

}

process generate_consensus_fasta {
    publishDir "${params.out_dir_int}/consensus", mode: 'copy', overwrite: true
    
    input:
        val(generate_consensus_script)
        val(ref_fasta)
        path vcf_file

    output:
        path "*consensus.vcf", emit: consensus_fasta_file

    script:
        """

        bash ${generate_consensus_script} ${ref_fasta} ${vcf_file}


        """
}



process generate_tree_data {

    publishDir "${params.out_dir_int}/tree", mode: 'copy', overwrite: true
    input:
        path collected_fasta

    output:
        path "tree.nhx", emit: tree_data
        path "multi-fasta_consensus.fasta", emit: multi_fasta
        path "multi-fasta_matrix.tsv", emit: snp_matrix

    script:
        """
        cat ${collected_fasta} >  multi-fasta_consensus.fasta
        snp-dists -b multi-fasta_consensus.fasta > multi-fasta_matrix.tsv
        sreformat stockholm multi-fasta_consensus.fasta > multi-fasta_consensus.stockholm
        quicktree multi-fasta_consensus.stockholm > tree.nhx
         

        """
}

workflow phylo_tree {
    take:
        generate_consensus_script_wf
        fastafile
        vcf_file_wf

    main:
        generate_consensus_fasta(generate_consensus_script_wf, fastafile, vcf_file_wf)
        collected_consensus = generate_consensus_fasta.out | collect
        generate_tree_data(collected_consensus)

    emit:
        generate_tree_data.out.tree_data

}


process drug_res_profiling {
    publishDir "${params.out_dir_int}/Drug_resistance_profiling", mode: 'copy', overwrite: true
    input:
        val dr_profiling_script
        path collected_preped_vcfs
        path collected_sorted_dr_res_catalogue
        val variant_caller
        val dp_cov
        val dr_res_scrpt_2
        val dr_res_script_3
        val dr_formatting_script

    output:
        path "all_dr_variants.json", emit: dr_res_json
        path "intergenic_dr_variants.json", emit: dr_res_int_json
        path "dr_profiling_summary.txt"
        path "detailed_drug_resistance_profile.tsv"
        path "short_drug_resistance_profile.tsv"

    script:
        """

        python ${dr_profiling_script} ${collected_preped_vcfs} ${collected_sorted_dr_res_catalogue} ${variant_caller} ${dp_cov} > "dr_profiling_summary.txt"

        Rscript ${dr_formatting_script} --files "all_dr_variants.json" "intergenic_dr_variants.json"
        
        """
}

process snp_typing {
    publishDir "${params.out_dir_int}/Snp_typing", mode: 'copy', overwrite: true

    input:
        val snp_typing_python_script
        path collected_preped_vcfs
        val prepped_lineage_file
        val variant_caller
        val dp_cov
        val dr_res_scrpt_2
        val dr_res_script_3
        val snp_typing_R_formating_script

    
    output:
        path "lineage.json", emit: lineage_json
        path "snp_typing_summary.txt"
        path "detailed_lineage_assignments.tsv"
    
    script:
    """

    python ${snp_typing_python_script} ${collected_preped_vcfs} ${prepped_lineage_file} ${variant_caller} ${dp_cov} > snp_typing_summary.txt

    Rscript ${snp_typing_R_formating_script} "--files" "lineage.json"
    """
}


/*process format_dr_lin_in_R {
    publishDir "${params.out_dir_int}/formatted_dr_lin_results", mode: 'copy', overwrite: true
    input:
        val main_formatting_script
        path dr_res_json
        path dr_res_int_json
        path lineage_json

    output:
        path "detailed_drug_resistance_profile.tsv"
        path "short_drug_resistance_profile.tsv"
        path "detailed_lineage_assignments.tsv"

    script:
        """
        Rscript ${main_formatting_script} "--files" ${dr_res_json} ${dr_res_int_json} ${lineage_json}

        """

}
*/