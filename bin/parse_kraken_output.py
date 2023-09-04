#!/bin/python

import argparse
import re

import pandas as pd

# Intantiate the Argument Parser

parser = argparse.ArgumentParser()

parser.add_argument('files', type =str, nargs = '*')

args = parser.parse_args()

files = args.files

# Append command line files into lists

kraken_ouput_file_list = []

for file in files:
    if re.search('[a-z]*_[0-9]*[i]?[a-z]?_kraken_report', file):
        kraken_ouput_file_list.append(file)
    else:
        continue

# Parse kraken output to extract basic information

def count_keys(d):
    num_keys = len(d)
    for k, v in d.items():
        if isinstance(v, dict):
            num_keys = num_keys + count_keys(v)
    return num_keys

def extract_family_to_species_info(kraken_output:str):
    unclassified_reads = 0
    perc_unclassified_reads = 0
    result_dict = {}
    current_f_key = None
    current_g_key = None
    current_s_key = None

    with open(kraken_output, 'r') as file:
        for line in file:
            data = line.rstrip().rsplit('\t')
            if data[3] == 'U':
                 unclassified_reads = int(data[1])
                 #print(unclassified_reads)
                 perc_unclassified_reads = float(data[0])
                 #print(perc_unclassified_reads)
            if data[3] == 'R':
                total_reads = unclassified_reads + int(data[1])

            if float(data[0]) != 0.00:
                if re.search("F", data[3]):
                        current_f_key = data[5]
                        current_g_key = None
                        current_s_key = None
                        result_dict[current_f_key] = [data[3], float(data[0]), {}]
                elif re.search("G[0-9]?", data[3]):
                        current_g_key = data[5]
                        current_g_pcent = float(data[0])
                        current_s_key = None
                        try:
                            result_dict[current_f_key][2][current_g_key] = [data[3], float(data[0]), {}]
                        except KeyError:
                             pass
                elif re.search("^S$", data[3]):
                    if current_g_key is not None:
                                current_s_key = data[5]
                                try:
                                    result_dict[current_f_key][2][current_g_key][2][current_s_key] = [data[3], float(data[0]), (float(data[0])/current_g_pcent)*100]
                                except KeyError:
                                     pass
                    
    return total_reads, unclassified_reads, perc_unclassified_reads, result_dict


def sort_result_dict_recursively(result_dict:dict):
    #print(result_dict)
    try:
        sorted_dict = sorted(result_dict.items(), key = lambda x: x[1][1], reverse = True)
    except AttributeError:
        sorted_dict = result_dict
    if sorted_dict is dict:
        for k, v in sorted_dict:
            if isinstance(v[2], dict):
                re_sorted_dict = v[2]
                try:
                    for key, val in re_sorted_dict:
                        sorted_dict = sort_result_dict_recursively(val)
                except ValueError:
                    pass
    else:
         pass


    return sorted_dict

def provide_family_to_species_result(kraken_output:str):
     total_reads, unclassified_reads, perc_unclassified_reads, result_dict = extract_family_to_species_info(kraken_output)
     sorted_dict = dict(sort_result_dict_recursively(result_dict))
     family_keys = [key for key, value in sorted_dict.items()]
     
     genus_level = sorted_dict[family_keys[0]][2]
     sorted_genus_level_dict = dict(sort_result_dict_recursively(genus_level))
     
     family_1_genus_keys = [key for key, value in sorted_genus_level_dict.items()]
     
     try:
        genus_1_species_level = sorted_genus_level_dict[family_1_genus_keys[0]][2]
     except IndexError:
          genus_1_species_level = "N/A"
     try:
         sorted_genus_1_species_level_dict = dict(sort_result_dict_recursively(genus_1_species_level))
         family_1_genus_1_species_keys = [key for key, value in sorted_genus_1_species_level_dict.items()]
     except ValueError:
         sorted_genus_1_species_level_dict = genus_1_species_level
         family_1_genus_1_species_keys = sorted_genus_1_species_level_dict

     try:
        genus_2_species_level = sorted_genus_level_dict[family_1_genus_keys[1]][2]
     except IndexError:
        genus_2_species_level = "N/A"

     try:     
        sorted_genus_2_species_level_dict = dict(sort_result_dict_recursively(genus_2_species_level))
     except ValueError:
          sorted_genus_2_species_level_dict = genus_2_species_level
     try:
        family_1_genus_2_species_level_keys = [key for key, value in sorted_genus_2_species_level_dict.items()]
     except AttributeError:
        family_1_genus_2_species_level_keys = sorted_genus_2_species_level_dict
    
     try:
          try:
            msg_1_f_1 = "/".join([family_1_genus_1_species_keys[0],sorted_genus_1_species_level_dict[family_1_genus_1_species_keys[0]][0]])
            psmsg_1_f_1 = sorted_genus_1_species_level_dict[family_1_genus_1_species_keys[0]][1]
            #g_2_f_1 = family_1_genus_keys[1], sorted_genus_level_dict[family_1_genus_keys[1]][0]
            psg_2_f_1 = sorted_genus_level_dict[family_1_genus_keys[1]][1]
          except TypeError:
            msg_1_f_1 = "N/A" 
            psmsg_1_f_1 = "N/A" 
            #g_2_f_1 = "N/A" 
            psg_2_f_1 = "N/A" 
     except IndexError:
          msg_1_f_1 = "N/A" 
          psmsg_1_f_1 = "N/A" 
          #g_2_f_1 = "N/A" 
          psg_2_f_1 = "N/A"  
          
    
     try:
          try:
             msg_2_f_1 = "/".join([family_1_genus_2_species_level_keys[0],sorted_genus_2_species_level_dict[family_1_genus_2_species_level_keys[0]][0]])
             psmsg_2_f_1 = sorted_genus_2_species_level_dict[family_1_genus_2_species_level_keys[0]][1]
          except TypeError:
              msg_2_f_1 = "N/A"
              psmsg_2_f_1 = "N/A"
     except IndexError:
          msg_2_f_1 = "N/A" 
          psmsg_2_f_1= "N/A" 
          
           
     try:
         g_2_f_1 = "/".join([family_1_genus_keys[1], sorted_genus_level_dict[family_1_genus_keys[1]][0]])
         pg_2_f_1 = sorted_genus_level_dict[family_1_genus_keys[1]][1]
     except IndexError:
         g_2_f_1 = "N/A"
         pg_2_f_1 = "N/A"
         
     try:
         g_1_f_1 = "/".join([family_1_genus_keys[0], sorted_genus_level_dict[family_1_genus_keys[0]][0]])
         pg_1_f_1 = sorted_genus_level_dict[family_1_genus_keys[0]][1]
     except IndexError:
         g_1_f_1 = "N/A"
         pg_1_f_1 = "N/A"


     
     reports = {
        "Total number of reads": total_reads,
        "Unclassified_reads":unclassified_reads,
        "percent of unclassified reads ": perc_unclassified_reads, 
        "Main Family 1": family_keys[0],
        "Percent supporting main family 1": sorted_dict[family_keys[0]][1],
        "Genus 1 in family 1": g_1_f_1,
        "Percent supporting genus 1 in family 1":pg_1_f_1,
        "Main species in genus 1 of family one" : msg_1_f_1,
        "Percent supporting main species in genus 1 of family one": psmsg_1_f_1,
        "Genus 2 in family 1":g_2_f_1,
        "Percent supporting genus 2 in family 1": pg_2_f_1,
        "Main species in genus 2 of family one": msg_2_f_1,
        "Percent supporting main species in genus 2 of family one": psmsg_2_f_1,
        "Main family 1 genus list":family_1_genus_keys
     }

     return reports




#print(provide_family_to_species_result(kraken_ouput_file_list[0]))

def parse_many_files(kraken_ouput_file_list: list):
     sample_level_dict = {}
     for each_file in kraken_ouput_file_list:
          sample_name = each_file.rsplit("_kraken_report")[0]
          sample_level_dict[sample_name] = provide_family_to_species_result(each_file)
     
     with open("parsed_kraken_report.tsv", 'w') as out_file:
          all_samples_data_df = pd.DataFrame.from_dict(sample_level_dict, orient = 'index')
          all_samples_data_df.to_csv(out_file, sep = '\t')


parse_many_files(kraken_ouput_file_list)
