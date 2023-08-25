'nextflow.enable=dsl2'

// Define parameters
params.command=""
params.run_mode=""
params.in_data_type=""
params.kraken_script="kraken_tx_id.sh"
params.genome_assembly_script="genome_assembly.sh"
params.trim_sample_script="trim_sample.sh"
params.emit_sam_script="emit_sam.sh"
params.coordsort_sam_script="coord_sort_sam.sh"
params.bamtofastq_script="bamtofastq.sh"
params.find_16S_hits_script="run_16S_tx_id.sh"
params.parse_kraken_script="parse_kraken_output.py"
params.parse_16S_script="parse_16S_output.py"
params.kraken_db_path=""
params.txdb_path=""
params.perc_cov=""
params.perc_id=""
params.seq_path=""
params.se_reads_u="nill"
//params.pe_reads=""
//params.ref_seq_path_d=""
params.ref_seq_path=""
//params.adp_path_d=""
params.adp_path=""
params.out_dir=""

if(params.in_data_type == "pe_illumina_reads"){
    params.pe_reads = params.seq_path
    params.pe_reads_u =  params.pe_reads + "/*fastq_{1,2}"
} else if(params.in_data_type == "se_illumina_reads" | params.in_data_type == "minion_ont_reads"){
    params.se_reads = params.seq_path
    params.se_reads_u = params.se_reads + "/*fastq"
}

if(params.adp_path == "Adapater_seq_file_path"){
    params.adp_path = adp_path_d
}
if(params.ref_seq_path == "Reference_file_path"){
    params.ref_seq_path = params.ref_seq_path_d
}


// Import inclusions

include {trim_fastq; kraken_contamination; find_16s_hit; genome_assembly; emit_sam; coordsort_sam; bamtofastq; parse_16S; parse_kraken} from "/home/dfgmrtc/Workflows/wf-module_templates/module_scripts/module_script.nf"


workflow comprehensive {
    take:
        kraken_script_wf
        genome_assembly_script_wf
        se_reads
        pe_reads
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
                kraken_contamination(kraken_script_wf, se_reads, pe_reads, kraken_database_path, input_read_type)
                genome_assembly(genome_assembly_script_wf, se_reads, pe_reads, input_read_type)
            } else {
                trim_fastq(trim_sample_script_wf, se_reads, pe_reads, adapter_seq, input_read_type)
                if(params.in_data_type == "pe_illumina_reads"){
                    trimmed_reads = trim_fastq.out | flatten | filter { it =~/P.fastq/ } | map { [it.name - ~/_[12]P.fastq/, it] } | groupTuple
                    kraken_contamination(kraken_script_wf, se_reads, trimmed_reads, kraken_database_path, input_read_type)
                    genome_assembly(genome_assembly_script_wf, se_reads, pe_reads, input_read_type)
                } else if (params.in_data_type == "se_illumina_reads"){
                    trimmed_reads = trim_fastq.out | filter { it =~/P.fastq/ }
                    kraken_contamination(kraken_script_wf, trimmed_reads, kc_pe_reads, kraken_database_path, input_read_type)
                    genome_assembly(genome_assembly_script_wf, trimmed_reads, pe_reads, input_read_type)
                }
            }
        } else if (params.run_mode == "with_reads_decontamination"){
            emit_sam(emit_sam_script_wf, se_reads, pe_reads, fastafile, input_read_type)
            coordsort_sam(coordsort_sam_script_wf, emit_sam.out)
            bamtofastq(bamtofastq_script_wf, coordsort_sam.out)
            if(params.in_data_type == "pe_illumina_reads"){
                reads_decontamination_out_seq = bamtofastq.out | flatten | map {[it.name - ~/.fastq_[12]/, it] } | groupTuple
                kraken_contamination(kraken_script_wf, se_reads, reads_decontamination_out_seq, kraken_database_path, input_read_type)
                genome_assembly(genome_assembly_script_wf, se_reads, reads_decontamination_out_seq, input_read_type)
            } else {
                reads_decontamination_out_seq = bamtofastq.out
                kraken_contamination(kraken_script_wf, reads_decontamination_out_seq, pe_reads, kraken_database_path, input_read_type)
                genome_assembly(genome_assembly_script_wf, reads_decontamination_out_seq, pe_reads, input_read_type)
            }
        }
        find_16s_hit(find_16S_hits_script_wf, genome_assembly.out, txdb_path_str)
        flattened_16S_ouputs = find_16s_hit.out
        //flattened_kraken_outputs = kraken_contamination.out
        parse_kraken(parse_kraken_script_wf, kraken_contamination.out)
        parse_16S(parse_16S_script_wf, flattened_16S_ouputs, perc_cov, perc_id)
        
}


workflow {
    comprehensive(Channel.value(params.kraken_script), Channel.value(params.genome_assembly_script), Channel.fromPath(params.se_reads_u), Channel.fromFilePairs(params.pe_reads_u), Channel.value(params.kraken_db_path),
    Channel.value(params.in_data_type), Channel.value(params.trim_sample_script), Channel.value(params.adp_path), Channel.value(params.emit_sam_script), Channel.value(params.ref_seq_path),
    Channel.value(params.coordsort_sam_script), Channel.value(params.bamtofastq_script), Channel.value(params.find_16S_hits_script), Channel.value(params.txdb_path),
    Channel.value(params.parse_kraken_script), Channel.value(params.parse_16S_script), Channel.value(params.perc_cov), Channel.value(params.perc_id))
}
