#!/bin/python


import argparse


import textwrap




# 1. Instantiate Parser
parser = argparse.ArgumentParser(prog = "taxo_id", 
formatter_class=argparse.RawDescriptionHelpFormatter,
description = textwrap.dedent('''\
                        dfgmrtc-labs taxonomic_id and read_decontamination workflow:
                        -------------------------------------------------------
This workflow takes sequencing reads, optionally decontaminates the reads, and runs three different but complementary algorithms to correctly
define the taxonomic identity - to species resolution - of organisms represented in the sequencing data
'''),
usage = '''
    taxo_id command [kraken_only, r16S_only, comprehensive, parse_files_only]
    --options/sub-command [all_output_files, kraken_output_files, 16S_output_files] --options  [--outdir_path]
    
Examples:
1.  To run a comprehensive taxonomic analysis on PE isolates, run: 
        taxo_id comprehensive [specify other options]
2.  To know which options are accepted by a command or sub-command, run:
        taxo_id command [sub-command] --help
3.  If the command has sub-command(s), run:
        taxo_id command sub-command --help
        to get more help about the options required by sub-command
4.  Provided you have kraken reports and you wish to parse the files without triggering a comprehensive analysis again, run:
        taxo_id parse_files_only process_kraken_outputs --kraken_output_files -options
      
      '''

)

#1b. Call argcomplete.autocomplete on parser
#argcomplete.autocomplete(parser)

# 2. Define Common argument to all commands and sub-commands
parser.add_argument("-o", "--out_dir_path",  help = "Specify an output directory path where results from this analysis will be saved. If directory already exists, hte directory and its content will be overwritten [default: %(default)s]", default="current working directory", metavar="/path/to/out_dir", dest="out_dir", type = str)
parser.add_argument('--version', action='version', version='%(prog)s 1.0')



# 3. Sub-parse parser to define positional commands
subparser = parser.add_subparsers(metavar = "commands:", dest = "command", required = True, help = "Commands must be specified after taxo_id")

# Parent parsers for subparsers
parent_parser_seq = argparse.ArgumentParser(add_help = False)
parent_parser_seq.add_argument("-s", "--seq_path", help = "Specify directory path containing input reads. Each file should end with 1/2.fastq(.gz) if PE reads or .fastq(.gz) if SE [default: %(default)s]", default="test_data_path", metavar = "path/to/fastq", dest = "fastq")
parent_parser_seq.add_argument("-d", "--in_data_type", help= "Available options include: pe_illumina_reads, se_illumina_reads, ont_minion_reads. [default: %(default)s]", default="pe_illumina_reads", metavar = "/path/to/fastq files", dest="in_data_type")
parent_parser_seq.add_argument("-m", "--run_mode", help="Avalaible options include: with_reads_decontamination, without_reads_decontamination [default: %(default)s]", default="without_reads_decontamination", metavar="specify option", dest="run_mode")
parent_parser_seq.add_argument("-o", "--outdir_path",  help = "Specify an output directory path where results from this analysis will be saved. If directory already exists, the directory and its content will be overwritten. [default: %(default)s]", default="current working directory", metavar="/path/to/out_dir", dest="out_dir", type = str)
parent_parser_seq.add_argument("-r", "--ref_seq_path", help="Specify path to directory where reference file is save. [default: %(default)s]", default="h37rv (NC_000962.3)",metavar="/path/to/ref_seq.fa", dest="ref_seq.fa")
parent_parser_adpt = argparse.ArgumentParser(add_help = False)
parent_parser_adpt.add_argument("-a", "--adpt_path", help = "Specify path where file containing adapter sequence is save. Not necessary if input data is from ONT minion [default: %(default)s]", default="nextera_pe_path", metavar = "/path/to/adapter.fa", dest = "adapters")
parent_parser_parse_files_params = argparse.ArgumentParser(add_help = False)
parent_parser_parse_files_params.add_argument("-c", "--perc_cov", default = 10, help = "Percentage coverage for blast+ algorithm  (default: %(default)s)", metavar="INT", dest = "perc_cov")
parent_parser_parse_files_params.add_argument("-i", "--perc_id", default = 98, help = "Percentage identity for blast+ algorithm (default: %(default)s)", metavar="INT", dest = "perc_id")
parent_parser_version = argparse.ArgumentParser(add_help = False)
parent_parser_version.add_argument("-v", '--version', action='version', version='%(prog)s 1.0')

## Kraken-only subparser
parser_kraken_only = subparser.add_parser("kraken_only", parents = [parent_parser_seq, parent_parser_adpt, parent_parser_version], help = "Command must be specified to run kraken2 for taxonomic identification only")
parser_kraken_only.add_argument("-db", "--kraken_db", help = "Workflow was tested with kraken microbial database. Specify path to kraken_db", metavar="/path/to/kraken_db", dest="kraken_db_path")


## 16S_only subparser
parser_16S_only = subparser.add_parser("r16S_only", parents = [parent_parser_seq, parent_parser_version], 
                                          help = "Command must be specified to run blast 16S rRNA taxonomic identification only")
parser_16S_only.add_argument("-t", "--taxdb", help="taxdb provides scientific and commnon names to taxids. Specify path to directory containing files taxdb.bti and taxdb.btd", metavar="/path/to/taxdb.*", dest="taxdb_path")
## Comprehensive subparser
parser_comprehensive = subparser.add_parser("comprehensive", parents = [parent_parser_seq, parent_parser_adpt, parent_parser_parse_files_params, parent_parser_version],
                                            help = "Default: Uses all available algorithms for taxonomic identification")
parser_comprehensive.add_argument("-db", "--kraken_db", help = "Workflow was tested with kraken microbial database. Specify path to kraken_db", metavar="/path/to/kraken_db", dest="kraken_db_path")
parser_comprehensive.add_argument("-t", "--taxdb", help="taxdb provides scientific and commnon names to taxids. Specify path to directory containing files taxdb.bti and taxdb.btd", metavar="/path/to/taxdb.*", dest="taxdb_path")

## Parse files only (Subparse subparser)

parser_parse_files = subparser.add_parser("parse_files_only",
                                            help = "Takes ouput files from kraken2 and or r16S analysis and parses the files to generate taxonomic classification").add_subparsers(metavar= "sub_commands", dest= "sub-command", help = "Sub commands are flags that don't take any value but rather consumes options")

#### Parse only kraken reports
parser_kraken_output_files = parser_parse_files.add_parser("process_kraken_outputs", help = "Flag must be specified to process kraken outputs", parents = [parent_parser_version])

parser_kraken_output_files.add_argument("-f", "--in_file_path", help="Specify path to directory containing output files from kraken analysis", metavar="path/to/in_files", dest="kraken_output_files")
parser_kraken_output_files.add_argument("-o", "--out_dir_path",  help = "Specify an output directory path where results from this analysis will be saved. If directory already exists, the directory and its content will be overwritten [default: %(default)s]", default="current working directory", metavar="/path/to/out_dir", dest="out_dir", type = str)
#parser_kraken_output_files.add_argument('--version', action='version', version='%(prog)s 1.0')

#### Parse only 16S outputs
parser_16S_output_files = parser_parse_files.add_parser("process_16S_outputs", help = "Flag must be specified to process r16S outputs", parents = [parent_parser_parse_files_params, parent_parser_version])
parser_16S_output_files.add_argument("-f", "--in_file_path", help="Specify path to directory containing output files from r16S analysis", metavar="path/to/in_files", dest="r16S_output_files")
parser_16S_output_files.add_argument("-o", "--out_dir_path",  help = "Specify an output directory path where results from this analysis will be saved. If directory already exists, the directory and its content will be overwritten [default: %(default)s]", default="current working directory", metavar="/path/to/out_dir", dest="out_dir", type = str)
### Parse all reports/outputs
parser_all_output_files = parser_parse_files.add_parser("process_all_outputs", help = "Default: Always true unless other parse_files_only flags are specified", parents = [parent_parser_parse_files_params, parent_parser_version])
parser_all_output_files.add_argument("--16S_in_file_path", help="Specify path to directory containing output files from r16S analysis", metavar="path/to/in_files", dest="r16S_output_files")
parser_all_output_files.add_argument("--kraken_in_file_path", help="Specify path to directory containing output files from kraken analysis", metavar="path/to/in_files", dest="kraken_output_files")
parser_all_output_files.add_argument("-o", "--out_dir_path",  help = "Specify an output directory path where results from this analysis will be saved. If directory already exists, the directory and its content will be overwritten [default: %(default)s]", default="current working directory", metavar="/path/to/out_dir", dest="out_dir", type = str)


#############################################################################################
## Parse argument parsers
args = parser.parse_args()





















'''
Important note: For help options, all commands and sub-commands are used as flags, i.e without -- prefix. However, when run on the
   command line for real or test data analysis, major commands are written as "-profile command", whereas all sub-commands and options are prefixed with -- 
 parent_parser_16S_str = argparse.ArgumentParser(add_help = False)
parent_parser_16S_str.add_argument("--16S_taxid_stratifier", action = "store_true", help = "Default:This flag stratifies r16S analysis by mtb taxid. To unset it, set flag to False", dest = "16S_stratifier")
 '''