#!/bin/python
from tabulate import tabulate
from collections import defaultdict
import os
import json

def remove_empty_conf_grad_drug_name(data):
    if isinstance(data, dict):
        filtered_dict = {}
        for key, value in data.items():
            if isinstance(value, dict):
                filtered_value = remove_empty_conf_grad_drug_name(value)
                if filtered_value:
                    filtered_dict[key] = filtered_value
            else:
                filtered_dict[key] = value
        return filtered_dict if filtered_dict else None
    return data

def merge_filtered_dict(*filtered_dicts):
    merged_dict = defaultdict(list)
    for dictionary in filtered_dicts:
        for key, value in dictionary.items():
            merged_dict[key].append(value)
    return dict(merged_dict)
                
    

        

def convert_filtered_dict_to_table(filtered_drug_res_dict, headers=None):
    table = []
    if isinstance(filtered_drug_res_dict, dict):
        for key, value in filtered_drug_res_dict.items():
            if isinstance(value, dict):
                subtable = convert_filtered_dict_to_table(value, headers + [key] if headers else [key])
                if subtable:
                    table.extend(subtable)
            
            else:
                table.append(headers + [key, value[:]] if headers else [key, value])

    return table


def create_table_from_filtered_dict(filtered_drug_res_dict):
    table = convert_filtered_dict_to_table(filtered_drug_res_dict)
    if not table:
        return None
    
    headers = ["Sample_name"] + ["Drug"] + ["Confidence_grading"] + ["Gene"] + ["Variant"] + ["Variant_attributes"]
    
    return tabulate(table, headers, tablefmt="grid")      

def create_table_from_filtered_dict2(filtered_drug_res_dict):
    table = convert_filtered_dict_to_table(filtered_drug_res_dict)
    if not table:
        return None
    
    headers = ["Sample_name"] + ["Reference"] + ["Lineage"] + ["Sub_lineage"] + ["position"] + ["Allele_change, lsp, spoligptype"]
    return tabulate(table, headers, tablefmt="grid")


def create_table_from_filtered_dict3(filtered_drug_res_dict):
    table = convert_filtered_dict_to_table(filtered_drug_res_dict)
    if not table:
        return None
    
    headers = ["Sample_name"] + ["Drug"] + ["Gene"] + ["variant"] + ["var_eff", "var_pos", "allele_change", "strand"]

    return tabulate(table, headers, tablefmt="grid")



class Drug_name_List(): # Converts user-specified drug-name list into a class object - A dictionary with the drug name as keys and mutation catalogue as values
    def __init__(self, my_per_drug_res_dict, drug_name_list):
        self.my_per_drug_res_dict = my_per_drug_res_dict
        self.drug_name_list = drug_name_list
        database_drug_list = ["AMI", "BDQ", "CAP", "CFZ", "CFZ", "DLM", "EMB", "ETH", "INH", "KAN", "LEV", "LZD", "MXF", "PZA", "RIF", "STM"]
        for drug_name in self.drug_name_list:
            
            if drug_name not in database_drug_list:
                assert "Drug name %s does not exist in our database" %drug_name
            else:
                self.drug_name_led_attr = {drug_name:self.my_per_drug_res_dict[drug_name] for drug_name in self.drug_name_list}

class Confidence_grading(Drug_name_List): # Convert the drug-name led class object derived above into another class object - A dictionary with drug names as keys, and another dictionary as values. 
    # The value dictionary have confidence grading as keys and a dictionary of genes, variants, positins as values
    def __init__(self, my_per_drug_res_dict, drug_name_list, confidence_grading_list):
        Drug_name_List.__init__(self, my_per_drug_res_dict, drug_name_list)
        self.confidence_grading_list = confidence_grading_list
        self.drugname_confidence_grading_led_attr = {drug_name:{} for drug_name in self.drug_name_list}
        for key in self.drugname_confidence_grading_led_attr:
            self.drugname_confidence_grading_led_attr[key] = {confidence_grading:self.drug_name_led_attr[key][confidence_grading] for confidence_grading in self.confidence_grading_list if confidence_grading in self.drug_name_led_attr[key].keys()}


class Sample_class(): # Coverts sample information into a class object
    def __init__(self, per_sample_drug_dict, sample_name_list):
        self.per_sample_drug_dict = per_sample_drug_dict
        self.sample_name_list = sample_name_list
        ## Process intergenic variants 
        for sample_name in self.sample_name_list:
            self.sample_led_int_variants_only = {sample_name: {gene_name: {variant:self.per_sample_drug_dict[sample_name][gene_name][variant]
                                                for variant in self.per_sample_drug_dict[sample_name][gene_name].keys()
                                                if self.per_sample_drug_dict[sample_name][gene_name][variant][0][0] == "intragenic_variant" or 
                                                self.per_sample_drug_dict[sample_name][gene_name][variant][0][0] == "intergenic_region"} 
                                                for gene_name in self.per_sample_drug_dict[sample_name].keys()}
                                                for sample_name in self.sample_name_list}
        if self.sample_led_int_variants_only:
                self.copied_int_variants_attr = self.sample_led_int_variants_only.copy()
                self.filtered_int_variants_attr =  remove_empty_conf_grad_drug_name(self.copied_int_variants_attr)
        else:
            pass
    
        
        ### Process sample VCFs with other variants

        for sample_name in self.sample_name_list:
            if sample_name in per_sample_drug_dict.keys():
                self.sample_led_dr_res_prof_attr = {sample_name:per_sample_drug_dict[sample_name] for sample_name in self.sample_name_list}
            else:
                print("sample %s vcf file not provided" %sample_name)


        ###### inhA-fabG1 relationship
        # Pending
        


class Assign_drug_res(Confidence_grading, Sample_class): # Sample name, drug name, and confidence gradings can be determined by 
    def __init__(self, my_per_drug_res_dict, drug_name_list, confidence_grading_list, per_sample_drug_dict, sample_name_list):
        Confidence_grading.__init__(self, my_per_drug_res_dict, drug_name_list, confidence_grading_list)
        Sample_class.__init__(self, per_sample_drug_dict, sample_name_list)
        self.sample_name_list = sample_name_list
        
        #### Process intergenic variants only
        self.samplename_drug_led_int_var_attr = {}
        self.novel_int_var = {}
        if self.filtered_int_variants_attr is not None:
            for int_sample_name in self.filtered_int_variants_attr.keys():
                for int_gene_name in self.filtered_int_variants_attr[int_sample_name].keys():
                    for int_var in self.filtered_int_variants_attr[int_sample_name][int_gene_name].keys():
                        for drug_name in self.drugname_confidence_grading_led_attr.keys():
                            for confidence_grading in self.drugname_confidence_grading_led_attr[drug_name].keys():
                                if int_gene_name in self.drugname_confidence_grading_led_attr[drug_name][confidence_grading].keys():
                                    for key, value in self.drugname_confidence_grading_led_attr[drug_name][confidence_grading][int_gene_name].items():
                                            for sel_var_attr in self.filtered_int_variants_attr[int_sample_name][int_gene_name][int_var]:
                                                if "/".join([key[0].upper(), key[-1].upper()]) == sel_var_attr[2]:        
                                                    if value[0] == sel_var_attr[1]:
                                                            self.samplename_drug_led_int_var_attr.setdefault(int_sample_name, {})
                                                            self.samplename_drug_led_int_var_attr[int_sample_name].setdefault(drug_name, {})
                                                            self.samplename_drug_led_int_var_attr[int_sample_name][drug_name].setdefault(confidence_grading, {})
                                                            self.samplename_drug_led_int_var_attr[int_sample_name][drug_name][confidence_grading].setdefault(int_gene_name, {})
                                                            self.samplename_drug_led_int_var_attr[int_sample_name][drug_name][confidence_grading][int_gene_name].setdefault(key, sel_var_attr)
        
                                                elif  "/".join([key[0].upper(), key[-1].upper()]) != sel_var_attr[2]:
                                                    if value[0] != sel_var_attr[1]:
                                                        self.novel_int_var.setdefault(int_sample_name, {})
                                                        self.novel_int_var[int_sample_name].setdefault(drug_name, {})
                                                        self.novel_int_var[int_sample_name][drug_name].setdefault(int_gene_name, {})
                                                        self.novel_int_var[int_sample_name][drug_name][int_gene_name].setdefault(sel_var_attr[2], sel_var_attr)
                                                        


                                                        #self.novel_int_var[int_sample_name][drug_name][int_gene_name][sel_var_attr[2]].append(sel_var_attr[0])
                                                        #self.novel_int_var[int_sample_name][drug_name][int_gene_name][sel_var_attr[2]].append(sel_var_attr[1])
                                                        #self.novel_int_var[int_sample_name][drug_name][int_gene_name][sel_var_attr[2]].append(sel_var_attr[2])
                                                        
                                                   

        self.long_table_dr_int_var_call = create_table_from_filtered_dict(self.samplename_drug_led_int_var_attr)

        self.long_novel_int_variant_call =  create_table_from_filtered_dict3(self.novel_int_var)
       

        ### Process other variants (p. and c. variants)
        self.samplename_drug_led_dr_res_prof_attr = {sample_name:{drug_name: {confidence_grading:{} 
                                                    for confidence_grading in self.confidence_grading_list} 
                                                    for drug_name in self.drug_name_list} 
                                                    for sample_name in self.sample_name_list}
        # Remove confidence grading
        self.drugname_led_no_conf_grading_attr = {}

        for drug_name in self.drugname_confidence_grading_led_attr.keys():
                for confidence_grading in self.drugname_confidence_grading_led_attr[drug_name].keys():
                        for gene_string in self.drugname_confidence_grading_led_attr[drug_name][confidence_grading].keys():
                            for variant_string in self.drugname_confidence_grading_led_attr[drug_name][confidence_grading][gene_string].keys():
                               #for pos in self.drugname_confidence_grading_led_attr[drug_name][confidence_grading][gene_string][variant_string]:
                                    self.drugname_led_no_conf_grading_attr.setdefault(drug_name, {})
                                    self.drugname_led_no_conf_grading_attr[drug_name].setdefault(gene_string, {})
                                    self.drugname_led_no_conf_grading_attr[drug_name][gene_string].setdefault(variant_string, self.drugname_confidence_grading_led_attr[drug_name][confidence_grading][gene_string][variant_string])
                                    
                                    

        
        #print(self.drugname_led_no_conf_grading_attr)
        self.samplename_drug_led_novel_variants = {sample_name:{drug_name: {}
                                                    for drug_name in self.drug_name_list} 
                                                    for sample_name in self.sample_name_list}
        

        for sample_name in self.samplename_drug_led_dr_res_prof_attr.keys():
                for drug_name in self.samplename_drug_led_dr_res_prof_attr[sample_name].keys():
                        for confidence_grading in self.samplename_drug_led_dr_res_prof_attr[sample_name][drug_name].keys():
                            #print(self.drugname_confidence_grading_led_attr[drug_name][confidence_grading])
                            
                            # Extract novel variants
                            self.samplename_drug_led_novel_variants[sample_name][drug_name] = {matched_gene_name:{matched_variant:self.sample_led_dr_res_prof_attr[sample_name][matched_gene_name][matched_variant]
                            
                            for matched_variant in self.sample_led_dr_res_prof_attr[sample_name][matched_gene_name].keys()
                            if self.sample_led_dr_res_prof_attr[sample_name][matched_gene_name][matched_variant][0][0] != "synonymous_variant"
                            if  matched_variant not in set(self.drugname_led_no_conf_grading_attr[drug_name][matched_gene_name].keys())
                            if len(matched_variant) != 1 
                            }
                            for matched_gene_name in self.sample_led_dr_res_prof_attr[sample_name].keys()
                            if matched_gene_name in self.drugname_confidence_grading_led_attr[drug_name][confidence_grading].keys()
                            }
                          
                            # Extract catalogued variants
                            self.samplename_drug_led_dr_res_prof_attr[sample_name][drug_name][confidence_grading] = {matched_gene_name:{matched_variant:self.sample_led_dr_res_prof_attr[sample_name][matched_gene_name][matched_variant][0]
                            for matched_variant in self.sample_led_dr_res_prof_attr[sample_name][matched_gene_name].keys()
                            #if any(list(map(lambda est_var: matched_variant in est_var, list(self.drugname_confidence_grading_led_attr[drug_name][confidence_grading][matched_gene_name].keys())))) is True
                            #if self.sample_led_dr_res_prof_attr[sample_name][matched_gene_name][matched_variant][1] in list(self.drugname_confidence_grading_led_attr[drug_name][confidence_grading][matched_gene_name].values())}
                            if  matched_variant in self.drugname_confidence_grading_led_attr[drug_name][confidence_grading][matched_gene_name].keys()
                            }
                            for matched_gene_name in self.sample_led_dr_res_prof_attr[sample_name].keys()
                            if matched_gene_name in self.drugname_confidence_grading_led_attr[drug_name][confidence_grading].keys()}
                            
                            
        # if self.per_sample_drug_dict[sample_name][matched_gene_name][matched_variant] != "intragenic_variant" or self.per_sample_drug_dict[sample_name][matched_gene_name][matched_variant] != "intergenic_region"}

        ### Debug             

        # Filter off Empty cells
        
        # 1. Catalogued variants
        self.copied_attr_for_long_dr_call =  self.samplename_drug_led_dr_res_prof_attr.copy()
        self.filtered_attr_for_long_dr_call = remove_empty_conf_grad_drug_name(self.copied_attr_for_long_dr_call)
        self.long_table_dr_call = create_table_from_filtered_dict(self.filtered_attr_for_long_dr_call)
        
        
        # 2. Novel variants
        self.copied_novel_variants = self.samplename_drug_led_novel_variants.copy()
        self.filtered_novel_variants = remove_empty_conf_grad_drug_name(self.copied_novel_variants)
        self.long_table_novel_var_call = create_table_from_filtered_dict3(self.filtered_novel_variants)
        
        
        


        
        #self.register = json.dumps(self.filtered_attr_for_long_dr_call)
        
        
        
        



    def __call__(self, drug_names:list, sample_names:list,confidence_list:list, return_only_p_c_call=False):
        self.filtered_attr_for_short_dr_call = self.filtered_attr_for_long_dr_call.copy()
        self.filtered_attr_for_short_dr_int_var_call = self.samplename_drug_led_int_var_attr.copy()
        self.filtered_novel_variants_for_call=self.filtered_novel_variants.copy()

        sample_names = tuple(sample_names)
        drug_names = tuple(drug_names)
        confidence_list = tuple(confidence_list)
        for sample in sample_names:
            for drug in drug_names:
                for confidence in confidence_list:
                    self.filtered_called_dr = {sample:{drug: {confidence:self.filtered_attr_for_short_dr_call[sample][drug][confidence]
                                    for confidence in confidence_list
                                    if confidence in self.filtered_attr_for_short_dr_call[sample][drug].keys()}
                                    for drug in drug_names
                                    if drug in self.filtered_attr_for_short_dr_call[sample].keys()}
                                    for sample in sample_names}
                    self.filered_int_var_called_dr = {sample:{drug: {confidence:self.filtered_attr_for_short_dr_int_var_call[sample][drug][confidence]
                                    for confidence in confidence_list
                                    if confidence in self.filtered_attr_for_short_dr_int_var_call[sample][drug].keys()}
                                    for drug in drug_names
                                    if drug in self.filtered_attr_for_short_dr_int_var_call[sample].keys()}
                                    for sample in sample_names}
                    self.filtered_novel_variant_called = {sample:{drug:self.filtered_novel_variants_for_call[sample][drug]
                                    for drug in drug_names
                                    if drug in self.filtered_novel_variants_for_call[sample].keys()}
                                    for sample in sample_names}
        
        
        self.short_table_dr_call =  create_table_from_filtered_dict(self.filtered_called_dr)
        self.short_table_int_var_dr_call =  create_table_from_filtered_dict(self.filered_int_var_called_dr)
        self.novel_variant = create_table_from_filtered_dict3(self.filtered_novel_variant_called)

        if return_only_p_c_call == True:
            return self.short_table_dr_call
        else:
            return self.short_table_dr_call, self.short_table_int_var_dr_call, self.novel_variant



                               

class Prepare_lineage_dict():
    def __init__(self, general_lin_dict, lineage_reference_list, main_lineages_list):
        self.general_lin_dict = general_lin_dict
        self.lineage_reference_list = lineage_reference_list
        self.main_lineages_list = main_lineages_list
        self.lin_ref_main_lin_led_attr = {lineage_ref: {main_lineage: self.general_lin_dict[lineage_ref][main_lineage]
                                        for main_lineage in self.main_lineages_list if main_lineage in self.general_lin_dict[lineage_ref]} for  lineage_ref in self.lineage_reference_list}
        #print(self.lin_ref_main_lin_led_attr)

class Assign_lineage(Prepare_lineage_dict):
    def __init__(self, general_lin_dict, lineage_reference_list, main_lineages_list, per_sample_lineage_dict, sample_name_list):
        Prepare_lineage_dict.__init__(self, general_lin_dict, lineage_reference_list, main_lineages_list)
        self.per_sample_lineage_dict = per_sample_lineage_dict
        self.sample_name_list = sample_name_list
    
        self.sample_ref_main_lin_led_attr = {sample_name:{lineage_ref:{main_lineage:{sub_lineage:{position:[self.per_sample_lineage_dict[sample_name][position], self.general_lin_dict[lineage_ref][main_lineage][sub_lineage][position][1], self.general_lin_dict[lineage_ref][main_lineage][sub_lineage][position][2]] 
                                            for position in  self.lin_ref_main_lin_led_attr[lineage_ref][main_lineage][sub_lineage].keys()
                                            if position in self.per_sample_lineage_dict[sample_name].keys()} 
                                            for sub_lineage in self.lin_ref_main_lin_led_attr[lineage_ref][main_lineage].keys()} 
                                            for main_lineage in main_lineages_list
                                            if main_lineage in self.general_lin_dict[lineage_ref]} 
                                            for lineage_ref in self.lineage_reference_list} 
                                            for sample_name in self.sample_name_list}
        
        self.sample_ref_main_lin_led_attr_long_call = self.sample_ref_main_lin_led_attr.copy()
        self.filtered_for_long_lin_call = remove_empty_conf_grad_drug_name(self.sample_ref_main_lin_led_attr_long_call)
        self.lineage_long_call = create_table_from_filtered_dict2(self.filtered_for_long_lin_call)

