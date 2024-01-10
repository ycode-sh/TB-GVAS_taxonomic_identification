#!/bin/python
import argparse
import re
import json
import os
import csv

from  extract_dr_res_lin_1 import process_many_vcf_files
from extract_dr_res_lin_2 import Assign_lineage
################# PARSE COMMAND LINE ARGUMENTS

parser = argparse.ArgumentParser()
parser.add_argument("all_files", type = str, nargs = "*")

args = parser.parse_args()

command_line_files = args.all_files

vcf_file_list = []
drug_resistance_file_list = []
lineage_file_list = []
sample_name_list = []
drug_name_list = []

sample_name_list_list = []

variant_caller = ""
dp_cov = ""

for any_file in command_line_files:
    if re.search("sample_[0-9]*[i]?[a-z]?_intersect.vcf", any_file): 
        vcf_file_list.append(any_file)
        file_name = os.path.basename(any_file.split("_intersect.vcf", 1)[0])  
        sample_name_list.append(file_name)
    elif re.search("lineage_snp_[0-9]*.tsv", any_file):
        lineage_file_list.append(any_file)
    elif re.search("bcftools", any_file):
        variant_caller =  any_file
    elif re.search("minos", any_file):
        variant_caller = any_file
    elif re.search("^[0-9]*$", any_file):
        dp_cov = any_file


#Pepare lineage annotation files

def define_lin_dict(lineage_file_list):
    lineage_file = lineage_file_list[0]
    lineage_reference_list = ["coll_2014", "freschi_2020", "shitikov_2017", "Napier_2020"]
    lin_dict = {lin_ref:{} for lin_ref in lineage_reference_list}

    with open(lineage_file, 'r') as lin_file:
        for line in lin_file:
            data = line.rstrip().rsplit("\t")
            lineage_main = str(data[2]).split(".", 1)[0]
            if lineage_main not in lin_dict[data[1]].keys():
                lin_dict[data[1]].setdefault(lineage_main, {})
                lin_dict[data[1]][lineage_main].setdefault(data[2], {})
                lin_dict[data[1]][lineage_main][data[2]].setdefault(str(data[0]), [str(data[3]), str(data[4]), str(data[5])])
                
            else:
                lin_dict[data[1]][lineage_main].setdefault(data[2], {})
                lin_dict[data[1]][lineage_main][data[2]].setdefault(str(data[0]), [str(data[3]), str(data[4]), str(data[5])])
          
    return lin_dict


############ RUN FUNCTIONS ##########################


general_lin_dict = define_lin_dict(lineage_file_list)
per_sample_drug_dict, per_sample_lineage_dict  = process_many_vcf_files(vcf_file_list, variant_caller, dp_cov)


################ VARIABLES ##################


main_lineages_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "La1", "asian_african_1", "asia_ancestral_2", "BOV_AFRI", "BOV", "asia_ancestral_3", "central_asia", 
                                      "asian_african_2/RD142", "pacific_RD150"]

lineage_reference_list = ["coll_2014", "freschi_2020", "shitikov_2017", "Napier_2020"]

#################### RUN CLASSES #########################

lineage_ass_instance = Assign_lineage(general_lin_dict, lineage_reference_list, main_lineages_list, per_sample_lineage_dict, sample_name_list)

##########  PRINT RESULTS TO STDOUT ##########

print(lineage_ass_instance.lineage_long_call)



##############  WRITE RESULT TO JSON  ###############

with open("lineage.json", "w") as json_lineage_file:
    json.dump(lineage_ass_instance.filtered_for_long_lin_call, json_lineage_file, indent = 4)
