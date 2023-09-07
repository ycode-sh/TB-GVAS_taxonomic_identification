'nextflow.enable=dsl2'


if(params.in_data_type == "pe_illumina_reads"){
    params.pe_reads = params.seq_path + "sample_*_{1,2}.fastq"
    params.pe_reads_channel = Channel.fromFilePairs(params.pe_reads)
    params.se_reads_channel = Channel.value("no_se_reads")
    
} else if(params.in_data_type == "se_illumina_reads" | params.in_data_type == "minion_ont_reads"){
    params.se_reads = params.seq_path + "sample_*.fastq"
    params.pe_reads_channel = Channel.value("no_pe_reads")
    params.se_reads_channel = Channel.fromPath(params.se_reads)
    
}


if(params.seq_path == "test_data_path"){
     exit 1, "No fastq files specified -- aborting"
}
if(params.adp_path == "nextera_pe_path"){
    params.adp_path_int = params.adp_path_d
}else {
    params.adp_path_int = params.adp_path
}
if(params.ref_seq_path == "h37rv.fasta"){
    params.ref_seq_path_int = params.ref_seq_path_d
} else {
    params.ref_seq_path_int = params.ref_seq_path
}
if(params.out_dir == "current_working_directory"){
    params.out_dir_int = params.out_dir_d
} else {
    params.out_dir_int = params.out_dir
}

// Import inclusions

include {trim_fastq; kraken_contamination; find_16s_hit; genome_assembly; emit_sam; coordsort_sam; pbamtofastq; parse_16S; parse_kraken} from "/home/dfgmrtc/Workflows/wf-module_templates/module_scripts/module_script.nf"

kraken_contamination
kraken_script_wf
        pe_reads_wf
        se_reads_wf
        kraken_db_path_wf
        input_read_type

trim_fastq {}
    take:
        trim_sample_script_wf
        pe_reads_wf
        adapter_seq
        input_read_type
        se_reads_wf

take:
        genome_assembly_script_wf
        pe_reads_wf
        se_reads_wf
        input_read_type

emit_sam_script_wf
        trimmed_reads
        fastafile
        input_read_type


workflow comprehensive {
    take:

    main:
        trim_fastq(trim_sample_script_wf, pe_reads_wf, se_reads_wf, kraken_db_path_wf, input_read_type)
        if (params.run_mode == "without_reads_decontamination"){
            kraken_contamination(kraken_script_wf, trim_fastq.out, kraken_db_path_wf, input_read_type)
            genome_assembly(genome_assembly_script_wf, trim_fastq.out, input_read_type)
        else if (params.run_mode == "with_reads_decontamination") {
            emit_sam(emit_sam_script_wf, trim_fastq.out, fastafile, input_read_type)
            coordsort_sam(coordsort_sam_script_wf, emit_sam.out)
            pbamtofastq(bamtofastq_script_wf, coordsort_sam.out, input_read_type)
            kraken_contamination(kraken_script_wf, pbamtofastq.out, kraken_db_path_wf, input_read_type)
            genome_assembly(genome_assembly_script_wf, pbamtofastq.out, input_read_type)
        }
        }
        find_16s_hit(find_16S_hits_script_wf, genome_assembly.out, txdb_path_str)
        flattened_16S_ouputs = find_16s_hit.out | flatten | collect
        flattened_kraken_outputs = kraken_contamination_se.out.kraken_reports | flatten | collect
        parse_kraken(parse_kraken_script_wf, flattened_kraken_outputs)
        parse_16S(parse_16S_script_wf, flattened_16S_ouputs, perc_cov, perc_id)


}

workflow comprehensive1 {
    take:
        kraken_script_wf
        genome_assembly_script_wf
        se_reads_wf
        pe_reads_wf
        kraken_database_path
        input_read_type
        trim_sample_script_wf
        adapter_seq
        emit_sam_script_wf
        fastafile
        coordsort_sam_script_wf
        bamtofastq_script_wf
        find_16S_hits_script_wf
        txdb_path_str
        parse_kraken_script_wf
        parse_16S_script_wf
        perc_cov
        perc_id

    main:
        if (params.run_mode == "without_reads_decontamination"){
            if (params.in_data_type == "minion_ont_reads"){
                kraken_contamination_se(kraken_script_wf, se_reads_wf, kraken_database_path, input_read_type)
                genome_assembly_se(genome_assembly_script_wf, se_reads_wf, input_read_type)
            } else {
                if(params.in_data_type == "pe_illumina_reads"){
                    trim_fastq_pe(trim_sample_script_wf, pe_reads_wf, adapter_seq, input_read_type)
                    trimmed_reads = trim_fastq_pe.out.trimmed_reads | flatten | filter { it =~/P.fastq/ } | map { [it.name - ~/_[12]P.fastq/, it] } | groupTuple
                    kraken_contamination_pe(kraken_script_wf, trimmed_reads, kraken_database_path, input_read_type)
                    genome_assembly_pe(genome_assembly_script_wf, trimmed_reads, input_read_type)
                } else if (params.in_data_type == "se_illumina_reads"){
                    trim_fastq_se(trim_sample_script_wf, se_reads_wf, adapter_seq, input_read_type)
                    kraken_contamination_se(kraken_script_wf, trim_fastq_se.out.trimmed_reads, kraken_database_path, input_read_type)
                    genome_assembly_se(genome_assembly_script_wf, trim_fastq_se.out.trimmed_reads, input_read_type)
                }
            }
        } else if (params.run_mode == "with_reads_decontamination"){
            if(params.in_data_type == "pe_illumina_reads"){
                emit_sam_pe(emit_sam_script_wf, pe_reads_wf, fastafile, input_read_type)
                coordsort_sam(coordsort_sam_script_wf, emit_sam_pe.out)
                bamtofastq(bamtofastq_script_wf, coordsort_sam.out, input_read_type)
                reads_decontamination_out_seq = bamtofastq.out | flatten | map {[it.name - ~/_[12].fastq/, it] } | groupTuple
                kraken_contamination_pe(kraken_script_wf, reads_decontamination_out_seq, kraken_database_path, input_read_type)
                genome_assembly_pe(genome_assembly_script_wf, reads_decontamination_out_seq, input_read_type)
            } else {
                emit_sam_se(emit_sam_script_wf, se_reads_wf, fastafile, input_read_type)
                coordsort_sam(coordsort_sam_script_wf, emit_sam_se.out)
                bamtofastq(bamtofastq_script_wf, coordsort_sam.out, input_read_type)
                kraken_contamination_se(kraken_script_wf, bamtofastq.out, kraken_database_path, input_read_type)
                genome_assembly_se(genome_assembly_script_wf, bamtofastq.out, input_read_type)
            }
        }
        if (params.in_data_type == "pe_illumina_reads"){
            find_16s_hit(find_16S_hits_script_wf, genome_assembly_pe.out, txdb_path_str)
            flattened_kraken_outputs = kraken_contamination_pe.out.kraken_reports | flatten | collect
        } else if (params.in_data_type == "se_illumina_reads" | params.in_data_type == "minion_ont_reads"){
            find_16s_hit(find_16S_hits_script_wf, genome_assembly_se.out, txdb_path_str)
            flattened_kraken_outputs = kraken_contamination_se.out.kraken_reports | flatten | collect
        }

        flattened_16S_ouputs = find_16s_hit.out | flatten | collect
        parse_kraken(parse_kraken_script_wf, flattened_kraken_outputs)
        parse_16S(parse_16S_script_wf, flattened_16S_ouputs, perc_cov, perc_id)
    emit:
        find_16s_hit.out
        
}


workflow {
    comprehensive(Channel.value(params.kraken_script), Channel.value(params.genome_assembly_script), params.se_reads_channel, params.pe_reads_channel, Channel.value(params.kraken_db_path),
    Channel.value(params.in_data_type), Channel.value(params.trim_sample_script), Channel.value(params.adp_path_int), Channel.value(params.emit_sam_script), Channel.value(params.ref_seq_path_int),
    Channel.value(params.coordsort_sam_script), Channel.value(params.bamtofastq_script), Channel.value(params.find_16S_hits_script), Channel.value(params.taxdb_path),
    Channel.value(params.parse_kraken_script), Channel.value(params.parse_16S_script), Channel.value(params.perc_cov), Channel.value(params.perc_id))
}
