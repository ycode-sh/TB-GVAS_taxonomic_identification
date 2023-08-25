#!/bin/bash

# $1 = main.nf

# Auto-display
./terminal_display.py

## Modes
I=""
M=""
H=""
Q=""

### Command line parameters defaults (Could all be replaced by user-selected/defined inputs)
sub_workflow_choice_d="Comprehensive Analysis"
sub_workflows_modes_choice_d="Without reads decontamination"
seq_dt_choice_d="Illumina paired end PE reads"
seq_reads_file_path_d="test_data_path"
ref_seq_file_path_d="h37rv.fasta"
adpt_seq_file_path_d="nextera_pe_path"
output_dir_path_d="current_working_directory"
parse_file_choice_d="parse_all_output_files"
parse_file_path_d="test_2_data_path"
percentage_cov_d=10
percentage_id_d=98

# Introductory Message
while true; do
    cat <<- _EOF_
### dfgmrtc-labs Taxonomic classification, Species identification, and Read decontamination Workflow ###
_________________________________________________________________________________________________________

This workflow takes sequencing reads, optionally decontaminates the reads, and runs three different but 
complementary algorithms to correctly define the taxonomic identity - to species resolution - of organisms 
represented in the sequencing data

    The workflow can be run in two modes: First is the interactive-dynamic mode. Here, new users with little 
    or no Bioinformatics experience can be guided through the workflow steps via an interactive prompt. 
    Responses are dynamically generated to user inputs and there is no need to manually type any command.
    The second mode requires users to specify commands, subcommands, and options  for each type of analysis. 
    This mode allows for more flexibility and user-control in adjusting default parameters to suit analysis.
    
    Please Select a mode:
    I. Interactive-dynamic 
    M. Manual command-based
    H. Display help messages
    Q. Quit workflow

_EOF_

read -ep "Select a mode (Case-insensitive): " mode
    case "$mode" in
    i|I) 
    echo "Interactive-dynamic mode initiated"
    I="Interactive" 
    
    break
    ;;
    m|M) 
    echo "Manual command-based mode initiated"
    M="Manual" 
    break
    ;;
    h|H) 
    echo "Showing a help display prompt ..."
    H="--help" 
    break
    ;;
    q|Q) 
    echo "Program terminated"
    Q="yes" 
    break
    ;;
    *)
    echo "Invalid choice: Make another choice!"
    esac
done

# Root supplied Variables/ Parameters
sub_workflows=("Comprehensive Analysis" "Uses all available algorithms for taxonomic identification and species  classification" "Kraken-only analysis" "Use only kraken analysis and custom programming for taxonomic identification and species  classification" "r16S blast+ only analysis" "Use only r16S analysis and custom programming for taxonomic identification and species  classification" "Parse files only" "Takes ouput files from kraken2 and or r16S analysis and parses the files to generate taxonomic classification" "Exit" "Exit this page")
sub_workflows_modes=("With reads decontamination" "This mode first runs a decontamination algorithm on the input sequences before species classification (Recommended choice, if contamination is suspected)" "Without reads decontamination" "Performs taxonomic identification and species classification on input sequences without attempting to decontaminate the reads")
seq_dt=("Illumina paired end PE reads" "" "Illumina single end SE reads" "" "ONT Minion reads" "")

# Associative array declarations

declare -A path_defaults_var
path_defaults_var["Fastq_files_directory_path"]="Fastq_files_directory_path"
path_defaults_var["Adapater_seq_file_path"]="Adapater_seq_file_path"
path_defaults_var["Reference_file_path"]="Reference_file_path"
path_defaults_var["Output_directory"]="Output_directory"

declare -A parse_files_params
declare -A file_paths

# Long Interactive Messages
interactive_choice_msg="An interactive dialog page will open once you click ENTER on your keyboard. Please read the information provided for each prompt carefully.
Some options in each of the prompts might have been preselected as defaults. You may change them if your analysis requirments are not compatible with the defaults"
file_parsing_params_directives="In the next prompt, you are (optionally) required to select file parsing parameters if your sub-workflow choice is either "Comprehensive" or r16S blast+ only analysis. You are free to replace the defaults to suitable values"
file_parsing_directives="Check to see that the defaults are suitable. If not, please input your choice defaults by navigating to rename and tapping your space character to edit the default. CLick ENTER/OK once done"
workflow_exit_directives="Do you wish to exit the interactive phase to view the help page?: Clicking Yes takes you to the main help page; Clicking no exits the program"
parse_file_choices_directives="Choose any one analysis type depending on the file type you you have avalaible. You may choose only option. "Process_all_outputs" is default"
sub_workflows_modes_choice_directives="Select a sub-workflow mode. Read the explanation to the right of each mode to get some context. Use arrow or coloured character key to navigate the options and click ENTER/OK to continue."
file_paths_spec_directives="The next few prompts will require you to specify paths to directories where WGS reads, adapter sequences, and reference files are stored. You will also need to specify path to output directory where analysis result should be saved. Useful defaults are provided for fastq files (if you want to run test), adapter sequence file, reference file, and output directory, You are free to select or not select any of the defaults provided based on your analysis requirements. Note that choices in this section are compulsory, that is, to run the workflow, you will have to provide your own files/file paths in the next prompt if you don't choose a default"
file_paths_checklist_directives="Select defaults. Use the arrow up or down keys to scroll across the options or press the character coloured in red to jump to any option. Once on any option of your choice, tap the space character to include such option. The former steps can also be performed once by pointing your mouse and clicking on any option of interest. You can select more than one default/options depending on your analysis requirements. Click ENTER/OK to proceed"
directory_selection_msg="In the next prompt, you will build/create absolute paths to directory where files are stored or analysis results should be saved. 
To navigate this space, use the arrow up or down keys to scroll across the files/directories and click on space character twice to select and include a file/directory in your path build.   
(DO NOT click ENTER until all the directories/files that make up your absolute path have been included. Clicking ENTER will confirm your build and 
capture whatever has been specified even if that is not the absolute paths to your file/directory). As an alternative to navigation using arrow keys, you may point 
your mouse and click on the directory/file of interest and then tap the space character to include your choice. As you increase your build, you will notice changes in the 
textarea widget below the displayed view. If you like, you may type-in the absolute paths manually using the text area widget provided. Once the absolute path to your 
files/directory has been specified click ENTER/OK to register your directory path build"
parse_file_only_sel_msg="Select or create directory where output files from this analysis should be saved. To navigate this space, use the arrow up or down keys to scroll across the files/directories and click on space character twice to select and include a file/directory in your path build.   
(DO NOT click ENTER until all the directories/files that make up your absolute path have been included. Clicking ENTER will confirm your build and 
capture whatever has been specified even if that is not the absolute paths to your file/directory). As an alternative to navigation using arrow keys, you may point 
your mouse and click on the directory/file of interest and then tap the space character to include your choice. As you increase your build, you will notice changes in the 
textarea widget below the displayed view. If you like, you may type-in the absolute paths manually using the text area widget provided. Once the absolute path to your 
files/directory has been specified click ENTER/OK to register your directory path build"
# Tags and their items
# pfc means parse file choices
pfc_tag_1="parse_all_output_files" 
pfc_item_1="Select this if you have unprocessed outputs from kraken2 and r16S analysis"
pfc_tag_2="parse_kraken_output_files"
pfc_item_2="Select this if you have unprocessed outputs from kraken2 analysis analysis only"
pfc_tag_3="parse_r16S_output_files"
pfc_item_3="Select this if you have unprocessed outputs from r16S analysis only"

# fpc means file_paths_checklist
fpc_tag_1="Fastq_files_directory_path" 
fpc_item_1="Use test fastq files" 
fpc_tag_2="Adapater_seq_file_path" 
fpc_item_2="Use default (Nextera PE)" 
fpc_tag_3="Reference_file_path" 
fpc_item_3="Use default (h37rv NC_0000926.3)" 
fpc_tag_4="Output_directory" 
fpc_item_4="Use default (Current working directory/folder)"

## Re-usable Functions
# General Functions
deplete_defaults(){
    local path_var=$1
    for chosen_default in "${option_array[@]}"; do
        if [ $path_var == "${chosen_default}" ]; then
            unset path_defaults_var[$path_var]
        fi
    done
}


message_print_function(){
    clear
    local message="$@"
    echo -e "\033[1mStatus message: $message\033[0m"
    echo "CLick ENTER/RETURN on your keyboard to continue or q to quit"
    while true; do
        read -r -n 1
        case $REPLY in
            "")
            echo "Proceeding.."
            break
            ;;
            q|Q)
            echo "uitting.."
            ./parser_python.py --help
            exit 0
            ;;
            *)
            ;;
        esac
    done
}

# Choice Processing Functions
choice_exit_wf(){
exit_choice=$(dialog --yesno "${workflow_exit_directives[@]}" 50 100  2>&1 >/dev/tty)
        clear
        ./parser_python.py --help
        exit 0
}

choice_parse_files(){
        message_print_function "Select file parsing options? Press ENTER to continue"
        parse_file_choice=$(dialog --radiolist "${parse_file_choices_directives[@]}" 100 200 200 "${pfc_tag_1[@]}" "${pfc_item_1[@]}" "on" "${pfc_tag_2[@]}" "${pfc_item_2[@]}" "off" "${pfc_tag_3[@]}" "${pfc_item_3[@]}" off 2>&1 >/dev/tty)
        clear
        message_print_function $directory_selection_msg
        if [[ $parse_file_choice == "parse_all_output_files" ]]; then
            message_print_function "Select directory paths where kraken2 outputs were saved?"
            kraken_in_file_path=$(dialog --fselect "$HOME" 100 200 2>&1 >/dev/tty)
            message_print_function "Select directory paths where r16S outputs were saved?"
            r16S_in_file_path=$(dialog --fselect "$HOME" 100 200 2>&1 >/dev/tty)
            clear
        
        elif [[ $parse_file_choice == "parse_kraken_output_files" ]]; then
            message_print_function "Select directory paths where kraken2 outputs were saved?"
            kraken_in_file_path=$(dialog --fselect "$HOME" 100 200 2>&1 >/dev/tty)
        else
            message_print_function "Select directory paths where r16S outputs were saved?"
            r16S_in_file_path=$(dialog --fselect "$HOME" 100 200 2>&1 >/dev/tty)
        fi 
        message_print_function "Select or create directory where output files from this analysis should be saved"
        parse_output_dir=$(dialog --fselect "$HOME" 100 200 2>&1 >/dev/tty)
        choice_parse_files_message="You selected ${kraken_in_file_path[@]} as directory containing kraken2 output files, ${r16S_in_file_path[@]} as directory containing r16S output files and ${parse_output_dir[@]} as directory where files from this analysis will be saved Next, You are required to select any one of the analysis types available for this workflow"
        message_print_function $choice_parse_files_message
        
        clear

        case $parse_file_choice in
            "parse_kraken_output_files")
            progress_msg_parse_file1=$(dialog --colors --msgbox "Progress report for Parse files function:
            The following information have been collected and will now be passed to the backend to run your selected analysis
                    Command option
            File parsing option: $parse_file_choice
                    File/Directory Paths
            Kraken output files: $kraken_in_file_path
            Outputs directory: $parse_output_dir
            To run this using the manual command-based mode:
            taxo_id  "parse_files_only" $parse_file_choice --kraken_in_file_path $kraken_in_file_path --out_dir_path $parse_output_dir" 100 100 2>&1 >/dev/tty)
            message_print_function "Your Analysis is running. Please wait..."
            ;;
            *)
            specify_optional_params
            percentage_cov=$(file_parsing_func "coverage" "$percentage_cov_d" 2>&1)
            percentage_id=$(file_parsing_func "identity" "$percentage_id_d" 2>&1)
            progress_msg_parse_file2=$(dialog --colors --msgbox "Progress report for Parse files function:
            
            The following information have been collected and will now be passed to the backend to run your selected analysis
                    Command option
            File parsing option: $parse_file_choice
                    File/Directory Paths
            Kraken output files:$kraken_in_file_path
            16S output files: $r16S_in_file_path
            Outputs directory: $parse_output_dir
            To run this using the manual command-based mode:
            taxo_id "parse_files_only" $parse_file_choice --kraken_in_file_path $kraken_in_file_path --16S_in_file_path $r16S_in_file_path --out_dir_path $parse_output_dir --perc_cov $percentage_cov --perc_id $percentage_id" 100 100 2>&1 >/dev/tty)
            message_print_function "Your Analysis is running. Please wait..."
            ;;
        esac
}

choice_main_workflow(){
        sub_workflow_choice_msg="You selected ${sub_workflow_choice_u} sub-workflow. This sub-workflow can be run in two modes: with or without reads decontamination. In the next prompt, you will select an option appropriate for your analysis run"
        message_print_function $sub_workflow_choice_msg
        sub_workflows_modes_choice_u=$(dialog --menu "${sub_workflows_modes_choice_directives[@]}" 100 200 200 "${sub_workflows_modes[@]}" 2>&1 >/dev/tty)
        sub_wf_choice_mode_seq_dt_intro_msg="You now have ${sub_workflow_choice_u[@]} sub-workflow ${sub_workflows_modes_choice_u[@]} mode on your buidlist. In the next prompt, you will provide some other information on your input data. Keep going"
        message_print_function $sub_wf_choice_mode_seq_dt_intro_msg
        
        seq_dt_choice_u=$(dialog --menu "Select sequencing platform from which inputs data are generated" 100 200 200 "${seq_dt[@]}" 2>&1 >/dev/tty)
        progress_msg=$(dialog --colors --msgbox "Progress Message 1:
        So far in your workflow build, the following options have been selected:
        sub_workflow_choice: ${sub_workflow_choice_u}
        sub_workflow_mode: ${sub_workflows_modes_choice_u}
        sequencing_platform and input data type: ${seq_dt_choice_u}" 100 100 2>&1 >/dev/tty)
        
        
}

# 3. File paths specification function
specify_file_paths(){
    for item in "${path_defaults_var[@]}"; do
        deplete_defaults $item    
    done
    for item in "${path_defaults_var[@]}"; do
            message_print_function "Are you ready to choose a directory/file path for $item ?"
            file_paths[$item]=$(dialog --fselect "$HOME" 100 200 2>&1 >/dev/tty)
            clear
            message_print_function "Selected filepath for $item: ${file_paths[$item]}"
    done
    message_print_function "Specify path to directory where kraken database files are saved. This is a compulsory selecion if you want to run kraken analysis. However,your choice will be ignored if r16S only is selected sub-workflow. You may press cancel to achieve a similar effect"
    kraken_db_path=$(dialog --fselect "$HOME" 100 200 2>&1 >/dev/tty)
    message_print_function "Selected filepath for kraken_db_path: ${kraken_db_path[$item]}"
    message_print_function "Specify path to directory where taxdb files are saved This is a compulsory selecion if you want to run r16S analysis. However,your choice will be ignored if r16S only is not selected sub-workflow. You may press cancel to achieve a similar effect"
    taxdb_files_path=$(dialog --fselect "$HOME" 100 200 2>&1 >/dev/tty)
    message_print_function "Selected filepath for taxdb_file_path: ${taxdb_files_path[$item]}"
}

# 4. Parameter definition function
specify_optional_params(){
    message_print_function $file_parsing_params_directives
    parse_files_params_cov=$(dialog --inputmenu "${file_parsing_directives[@]}" 100 200 200 "Percentage_coverage" 10 2>&1 >/dev/tty)
    parse_files_params["coverage"]=$parse_files_params_cov
    parse_files_params_id=$(dialog --inputmenu "${file_parsing_directives[@]}" 100 200 200 "Percentage_identity" 98 2>&1 >/dev/tty)
    parse_files_params["identity"]=$parse_files_params_id

    clear
}

# 5. Process user file paths and params inputs function
process_val_inputs(){
message_print_function $file_paths_spec_directives
        file_paths_checklist=$(dialog --checklist  "${file_paths_checklist_directives[@]}" 100 200 200 "${fpc_tag_1[@]}" "${fpc_item_1[@]}" off "${fpc_tag_2[@]}" "${fpc_item_2[@]}" "off" "${fpc_tag_3[@]}" "${fpc_item_3[@]}" "off" "${fpc_tag_4[@]}" "${fpc_item_4[@]}" "off" 2>&1 >/dev/tty)
        clear
        IFS=" " read -r -a option_array <<< "${file_paths_checklist}"
        if [ "${#option_array[@]}" -eq 0 ]; then #No defaults
            message_print_function "No default chosen"
            message_print_function $directory_selection_msg
            specify_file_paths
            specify_optional_params
        
        elif [ "${#option_array[@]}" -eq 4 ]; then #All defaults
            message_print_function "All defaults selected. Continue to see them?"
            message_print_function "Selected Defaults: ${option_array[@]}"
            specify_file_paths
            specify_optional_params
        
        else  #Some defaults
            message_print_function "Selected Defaults: ${option_array[@]}. 
            You will now select directory path for unselected default options ..."
            message_print_function $directory_selection_msg
            specify_file_paths
            specify_optional_params
        
        fi
}




## Command-line arameter Resolution

# Associative arrays
declare -A command_option_params=(
    ["Comprehensive Analysis"]="comprehensive"
    ["Kraken-only analysis"]="kraken_only"
    ["r16S blast+ only analysis"]="r16S_only"
    ["With reads decontamination"]="with_reads_decontamination"
    ["Without reads decontamination"]="without_reads_decontamination"
    ["Illumina paired end PE reads"]="pe_illumina_reads"
    ["Illumina single end SE reads"]="se_illumina_reads"
    ["ONT Minion reads"]="minion_ont_reads"
)

declare -A file_path_params=(
    ["Fastq_files_directory_path"]="seq_reads_file_path"
    ["Adapater_seq_file_path"]="adpt_seq_file_path"
    ["Reference_file_path"]="ref_seq_file_path"
    ["Output_directory"]="output_dir_path"
)

declare -A file_path_params_def=(
    ["seq_reads_file_path"]="test_data_path"
    ["ref_seq_file_path"]="h37rv.fasta"
    ["adpt_seq_file_path"]="nextera_pe_path"
    ["output_dir_path"]="current_working_directory"
)


declare -A final_array

### Resolution functions
command_option_func(){
    local set_default=$1
    local user_input=$2
    if [ -z "${user_input[@]}" ]; then
        cl_var="${command_option_params[$set_default]}"
        echo "$cl_var"
 
    else
        #echo "User supplied inputs" 
        for key in "${!command_option_params[@]}"; do
            if [[ $key == "${user_input}" ]]; then
                cl_var="${command_option_params[$key]}"
                echo "$cl_var"
                
            fi
        done
    fi
}

file_paths_func(){
for item1 in "${!file_paths[@]}"; do
    for item2 in "${!file_path_params[@]}"; do
        if [[ $item1 == $item2 ]]; then
            item3="${file_path_params[$item1]}"
            item4="${file_paths[$item1]}"
            final_array[$item3]=$item4
            unset file_path_params[$item1]
        fi
    done
done
for item5 in "${file_path_params[@]}"; do
    item6="${file_path_params_def[$item5]}"
    final_array[$item5]=$item6
done
}



file_parsing_func(){
    input=$1
    if echo "${parse_files_params[$input]}" | grep -q "RENAMED" && echo "${parse_files_params[$input]}" | grep -E -q "[0-9]{2,3}" ; then
        output=$(echo "${parse_files_params[$input]}" | awk '{print $NF}')
        
    else
        output=$2
        
    fi
    echo "$output"
}    

# Resolve 
# 1. Command Options
sub_workflow_choice=$(command_option_func "$sub_workflow_choice_d" "$sub_workflow_choice_u" 2>&1)
sub_workflows_modes_choice=$(command_option_func "$sub_workflows_modes_choice_d" "$sub_workflows_modes_choice_u" 2>&1)
seq_dt_choice=$(command_option_func "$seq_dt_choice_d" "$seq_dt_choice_u" 2>&1)

#2. File/Directory Paths
file_paths_func


## Choices logic Processing


if [ -n "$Q" ]; then # If user decides to quit program
    ./terminal_display.py
    echo "Good bye"
elif [ -n "$I" ]; then # User selects interactive promts
    message_print_function $interactive_choice_msg
    sub_workflow_choice_u=$(dialog --colors --menu "Select the sub-workflow suitable for your research question. Use the arrow up and down keys, or point your mouse and click on any choice, or press any of the character marked in red to select your option. Press the ENTER button to register your selection and proceed" 100 200 200 "${sub_workflows[@]}" 2>&1 >/dev/tty)
    case $sub_workflow_choice_u in
        "Exit")
        choice_exit_wf 
        ;;
        "Parse files only")
        choice_parse_files
        ;;
        *)
        choice_main_workflow
        process_val_inputs
        percentage_cov=$(file_parsing_func "coverage" "$percentage_cov_d" 2>&1)
        percentage_id=$(file_parsing_func "identity" "$percentage_id_d" 2>&1)
        progress_msg_alg=$(dialog --colors --msgbox "Progress report 2: 
        The following information have been collected and will now be passed to the backend to run your selected analysis
                        Command Options
        sub_workflow_choice: $sub_workflow_choice 
        Reads decontamination choice: $sub_workflows_modes_choice
        Sequencing platform/input data type: $seq_dt_choice
                        File/Directory Paths 
        Fatsq files: "${final_array["seq_reads_file_path"]}"
        Reference sequence file: "${final_array["ref_seq_file_path"]}"
        Adapter sequence file: "${final_array["adpt_seq_file_path"]}"
        Outputs directory: "${final_array["output_dir_path"]}"
        kraken_databse_dir:$kraken_db_path
        taxdb_files_dir:$taxdb_files_path
                        File-parsing parameters
        Percentage coverage: $percentage_cov
        Percentage identity: $percentage_id
        
        To run this using the manual command-based mode:
        taxo_id $sub_workflow_choice --run_mode $sub_workflows_modes_choice --in_data_type $seq_dt_choice --seq_path "${final_array["seq_reads_file_path"]}" --adp_path "${final_array["adpt_seq_file_path"]}" --ref_seq_path "${final_array["ref_seq_file_path"]}" --kraken_db_path $kraken_db_path --taxdb_path $taxdb_files_path --out_dir "${final_array["output_dir_path"]}" --perc_cov $percentage_cov --perc_id $percentage_id
        
        Note: You don't need to specify kraken database file path if your sub-workflow choice is r16S only. Likewise, don't specify taxdb if you are running only kraken analysis." 100 100 2>&1 >/dev/tty) 
        message_print_function "Your Analysis is running. Please wait..."
        nextflow run $1 --command $sub_workflow_choice --run_mode $sub_workflows_modes_choice --in_data_type $seq_dt_choice --kraken_db_path $kraken_db_path --taxdb_path $taxdb_files_path \
        --perc_cov $percentage_cov --perc_id $percentage_id --seq_path "${final_array["seq_reads_file_path"]}" --adp_path "${final_array["adpt_seq_file_path"]}" --ref_seq_path "${final_array["ref_seq_file_path"]}"
        ;;
    esac
elif [ -n "$M" ]; then
    echo "Study the help page"
    ./parser_python.py --help \n
    read -ep "Carefully specify your command-line arguments. Make sure you have checked the help page and the "interactive-dynamic" interface before using the manual command-line interface:
    cli: "


    cli_args=$REPLY
    echo "Running your specified workflow..."
    exit 0
else
    [ -n "$H" ]
    ./terminal_display.py
    ./parser_python.py --help
    exit 0
fi

