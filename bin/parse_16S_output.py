#!/bin/python

import re
import argparse
import pandas as pd

# 1. Instantiate the argument parser object and intialize an empty list where files will br storee 
parser = argparse.ArgumentParser()

parser.add_argument("file", type=str, nargs='*')
args = parser.parse_args()
files = args.file
params = []

blast_16s_output_files_list = []

# 2. Parse positional parameters and save files in a list 
for file in files:
    if re.search("[a-z]*_[0-9]*[i]?[a-z]?_16s_output.tsv", file):
        blast_16s_output_files_list.append(file)
    elif re.search("^[0-9]{2,3}$", file):
        params.append(file)
    else:
        print("no file exists")


p_id = float(params[1])
p_cov = float(params[0])

#t_str = True


# 3. define a function that filters blast output and extract important parameters like 
# evalue, taxid, subject name etc into a dictionary

# This function is mtb stratification-aware
def extract_info(blast_output: str, p_cov: float, p_id: float):
    taxa_list = []
    subjects_info = {}
    with open(blast_output, 'r') as blast:
        for line in blast:
            if line.startswith("#"):
                continue
            else:
                data = line.rstrip().rsplit('\t')
                if data[14] not in taxa_list:
                    taxa_list.append(data[14])
                if float(data[2]) >= p_id and float(data[10]) >=p_cov: #and float(data[11]) >= 60:
                    if data[12] not in subjects_info:
                        subjects_info[data[12]] = [[data[14]], [data[9]], [data[15]]]
                    else:
                        subjects_info[data[12]][0].append(data[14])
                        subjects_info[data[12]][1].append(data[9])
                        subjects_info[data[12]][2].append(data[15])

                    
    return subjects_info, taxa_list


# 4. Define a function that modifies and sorts the subjects_info and taxa_list. The function takes a single argument, and calls function extract info internally)
def update_subject_info(sample_file:str, p_cov: float, p_id: float):
    subjects_info, taxa_list = extract_info(sample_file, p_cov, p_id)
    nested_non_unique_taxids = []
    flat_non_unique_taxids = []
    for key in subjects_info.items():
        item = key[0]
        nested_non_unique_taxids.append(subjects_info[item][0])
        flat_non_unique_taxids = [item for sublist in nested_non_unique_taxids for item in sublist]
    unique_taxids = taxa_list
    taxid_frequency_dict = {taxid: flat_non_unique_taxids.count(taxid) for taxid in unique_taxids}
    for key in subjects_info:
        subjects_info[key] = [list(set(lst))[0] if len(set(lst)) == 1 else lst for lst in subjects_info[key]]

    updated_subject_info = {key:value + [taxid_frequency_dict[str(value[0])]] if str(value[0]) in taxid_frequency_dict else value for key, value in subjects_info.items()}
    updated_subject_info = {k: [v[0], sorted(v[1]) if isinstance(v[1], list) else v[1], v[2], v[3]] for k, v in updated_subject_info.items()}
    updated_subject_info = dict(sorted(updated_subject_info.items(), key=lambda x: x[1][-1], reverse=True))
    return updated_subject_info

# 5. Define a function that returns best hits from extracted data
def return_best_hits(sample_file:str, p_cov: float, p_id: float):
    updated_subject_info = update_subject_info(sample_file, p_cov, p_id)
    subject_name_list = [key for (key, value) in updated_subject_info.items()]
    
    try:
        scn_1 = subject_name_list[0]
        st_1 = updated_subject_info[subject_name_list[0]][0]
        sdht_1 = updated_subject_info[subject_name_list[0]][3]
        se_1 = updated_subject_info[subject_name_list[0]][1]
        scn_2 = subject_name_list[1]
        st_2 = updated_subject_info[subject_name_list[1]][0]
        sdht_2 = updated_subject_info[subject_name_list[1]][3]
        se_2 = updated_subject_info[subject_name_list[1]][1]
        scn_3 = subject_name_list[2]
        st_3 = updated_subject_info[subject_name_list[2]][0]
        sdht_3 = updated_subject_info[subject_name_list[2]][3]
        se_3 = updated_subject_info[subject_name_list[2]][1]
    except IndexError:
        scn_1 = "NA"
        st_1 = "NA"
        sdht_1 = "NA"
        se_1 = "NA"
        scn_2 = "NA"
        st_2 = "NA"
        sdht_2 = "NA"
        se_2 = "NA"
        scn_3 = "NA"
        st_3 = "NA"
        sdht_3 = "NA"
        se_3 = "NA"



    best_hit_params = {
        "Subject 1 scientific name":scn_1,
        "Subject 1 Taxid":st_1,
        "Subject 1 Database Hit Frequency":sdht_1,
        "Subject 1 evalue":se_1,
        "Subject 2 Scientific name":scn_2,
        "Subject 2 Taxid":st_2,
        "Subject 2 Database Hit Frequency":sdht_2,
        "Subject 2 evalue":se_2,
        "Subject 3 Scientific name":scn_3,
        "Subject 3 Taxid":st_3,
        "Subject 3 Database Hit Frequency":sdht_3,
        "Subject 3 evalue":se_3

    }
    return best_hit_params
    

# 6. Define a function that parses each sample and saves parameters as tsv file   
def parse_samples(blast_16s_output_files_list:list, p_cov: float, p_id: float):
    per_sample_parameters = {}
    for sample_file in blast_16s_output_files_list:
        if sample_file.endswith('output.tsv'):
            sample_name = sample_file.rsplit("_16s_output.tsv")[0]
            best_hit = return_best_hits(sample_file, p_cov, p_id)
            per_sample_parameters[sample_name] = best_hit

        else: 
            continue        
    per_sample_parameters_df = pd.DataFrame.from_dict(per_sample_parameters, orient='index')
    with open("parsed_16S_samples_output.tsv", 'w') as output_file:
        per_sample_parameters_df.to_csv(output_file, sep = '\t')
    return output_file


parse_samples(blast_16s_output_files_list, p_cov, p_id)
