#!/bin/python

import re
from extract_dr_res_lin_1 import convert_p_str, convert_c_str, convert_n_string


def process_variant_string(string, allele_change, variant_annotation):
            function_level_dict = {"frameshift_variant":"run_c_string_func", "missense_variant":"run_p_string_func",
                                    "synonymous_variant":"run_p_string_func", "upstream_gene_variant": "run_c_string_func",
                                    "downstream_gene_variant":"run_c_string_func",
                                    "intergenic_region":"run_n_string_func", "intragenic_variant":"run_n_string_func"}

            if function_level_dict[variant_annotation] == "run_p_string_func": 
                 var = convert_p_str(string)
            elif function_level_dict[variant_annotation] == "run_c_string_func":
                 var = convert_c_str(string, allele_change, variant_annotation)
            elif function_level_dict[variant_annotation] == "run_n_string_func":
                 var = convert_n_string(string)
            else:
                 print("This variant_annotation:%s is not captured by our software" %variant_annotation)
            return var
    
def process_bcf_vcf(data, extracted_info:dict):
    info_tags = {"IMF":"", "DP":"", "MQ":"", "DP4":""}
    variant_info = data[7].split("ANN")[0].split(";")[:]
    for info_tag_idx in range(len(variant_info[:])):
        if variant_info[info_tag_idx].split("=")[0] in info_tags.keys():
            variant_info[info_tag_idx]
            specific_tag = variant_info[info_tag_idx].split("=")[0]
            tag_value = variant_info[info_tag_idx].split("=")[1]
            info_tags[specific_tag] = tag_value
            #print(info_tags)
    FRS = (float(info_tags["DP4"].split(",")[2]) + float(info_tags["DP4"].split(",")[3]))/float(info_tags["DP"])
    # float(info_tags["IMF"]) < 0.75   or FRS < 0.45
    if float(info_tags["DP"]) < 20  or float(info_tags["MQ"]) < 20:
        pass
    else:
        annotation_info = data[7].split("ANN=")[1]
        extracted_info.setdefault(data[1], {})
        extracted_info[data[1]]["allele_change"] = "/".join([data[3], data[4]])
        #extracted_info[data[1]]["allele_balance"] = round(float(info_tags["DP4"].split(",")[2])/float(info_tags["DP4"].split(",")[3]), 2)
        # Extract Lineage information
        if re.search(r'(?:lineage=)', annotation_info):
            lineage_info = {"lineage":"", "exp_allele":"", "lsp":"", "spoligotype":"", "lin_ref":""}
            lineage_info_regexp = {"lineage":re.compile(r'(?:lineage)=((\d\.?)*)'), "exp_allele":re.compile(r'(?:exp_allele=)((\w)*)'), "lsp":re.compile(r'(?:lsp=)((\w_?-?\s?,?)*)'), "spoligotype":re.compile(r'(?:spoligotype=)((\w_?-?\s?,?)*)'), "lin_ref":re.compile(r'(?:lin_ref=)((\w_?,?)*)')}
            for l_item in lineage_info_regexp.keys():
                
                match = lineage_info_regexp[l_item].search(annotation_info)
                if match:
                    if "NA" in match.group(1):
                        trimmed_lineage_info = match.group(1).replace("NA", ",").replace(",", "")
                        lineage_info[l_item] = trimmed_lineage_info
                    else:
                        lineage_info[l_item] = match.group(1)
            if lineage_info["exp_allele"] == extracted_info[data[1]]["allele_change"].split("/")[1]:
                extracted_info[data[1]]["lineage_info"] = lineage_info
            else:
                pass
                
        # Extract AMR information
        if re.search(r'(?:antibiotics_gene=)', annotation_info):
            amr_info = {"antibiotics_gene":"", "related_antibiotics":""}
            amr_info_regexp = {"antibiotics_gene":re.compile(r'(?:antibiotics_gene=)((\w_?,?)*)'), "related_antibiotics":re.compile(r'(?:related_antibiotics=)((\w,?_?)*)')}
            for amr_item in amr_info_regexp.keys():
                match = amr_info_regexp[amr_item].search(annotation_info)
                if match:
                    amr_info[amr_item] = match.group(1)
            
            extracted_info[data[1]]["amr_info"] = amr_info
        # Extract variable region information
        if re.search(r'(?:vr_locus_tag=)', annotation_info):
            vr_info = {"vr_locus_tag":"", "comment":""}
            vr_info_regexp = {"vr_locus_tag":re.compile(r'(?:vr_locus_tag=)(\w*)'), "comment":re.compile(r'(?:comment=)((\w_?)*)')}
            for vr_item in vr_info_regexp.keys():
                match = vr_info_regexp[vr_item].search(annotation_info)
                if match:
                    vr_info[vr_item] = match.group(1)
            extracted_info[data[1]]["vr_info"] = vr_info
        # Extract gene annotation information
        if re.search(r'(?:gene_ann_region=)', annotation_info):
            gene_annotation_dict = {}
            gene_annotation_id = {"gene":[2, 3], "CDS":[5, 6, 7], "repeat_region":[1, 6], "pseudogene":[2, 3, 6], "ncRNA":[1, 2], "exon":[5, 6], "sequence_feature":[1, 2, 6], "mobile_genetic_element":[1, 4, 6], "tRNA":[1, 2, 5, 6], "transcript":[1, 2, 5, 6], "rRNA":[1, 2, 5, 6]}
            gene_annotation_regexp = [re.compile(r'(?:Name=)((\w_?)*)'), re.compile(r'(?:gbkey=)((\w_?)*)'), re.compile(r'(?:gene=)(\w*)'), re.compile(r'(?:gene_biotype=)((\w_?)*)'), 
                                  re.compile(r'(?:mobile_element_type=)((\w_?\:?)*)'), re.compile(r'(?:product=)((\w_?-?\.?\(?\)?\/?)*)'), 
                                  re.compile(r'(?:Note=)((\w_?\w?\'?-?\.?\(?\)?\.?\/?\%?\:?\'?)*)'), re.compile(r'(?:inference=protein_motif\:PROSITE\:)(\w*)')]
            available_gene_annotation_regions = annotation_info.split("gene_ann_region=")[1].split(";")[0].split(",")[:]
            for reg in available_gene_annotation_regions:
                match_regexp_list = [gene_annotation_regexp[idx] for idx in gene_annotation_id[reg]]
                gene_annotation_dict.setdefault(reg, [re_pattern.search(annotation_info).group(1) for re_pattern in match_regexp_list
                                                    if re_pattern.search(annotation_info) is not None])

            extracted_info[data[1]]["gene_annotation_info"] = gene_annotation_dict
            pmid_pattern = re.compile(r'(?:PMID:)(\d*)')
            if pmid_pattern.findall(annotation_info):
                 extracted_info[data[1]]["gene_annotation_info"]["CDS"].append(pmid_pattern.findall(annotation_info))
        # Extract info on snpEff annotation
        #snpeFF_ann = {"variant_allele":None, "variant_position":"", "variant_annotation":"", "variant_effect":"", "affected_locus":None, "variant_type":""}         
        snpeFF_ann = {"variant_allele":None, "variant_position":None, "variant_description":{}}
        allele_change = "/".join([data[3], data[4]])
        try:
            annotation_info.split(";")[1]
            all_snpEff_annot = annotation_info.split(";")[0].split(",")[:]
            
        except IndexError:
             all_snpEff_annot = annotation_info.split(",")

        for ann in all_snpEff_annot:
            variant_type = ""
            ann_var_list = ann.split("|")[:]
            variant_annotation = ann_var_list[1]
            if ann_var_list[1] in ["frameshift_variant", "upstream_gene_variant", "downstream_gene_variant", "intergenic_region", "intragenic_variant"]:
                 if ann_var_list[9].split(".", 1)[0] == "c":
                      hgvs_c = ann_var_list[9]
                      variant_type = process_variant_string(hgvs_c, allele_change, variant_annotation)
                 elif ann_var_list[9].split(".", 1)[0] == "n":
                      hgvs_n = ann_var_list[9]
                      variant_type = process_variant_string(hgvs_n, allele_change, variant_annotation)
                 else:
                      hgvc_id = ann_var_list[9]
                      print("%s is not captured by my code" %hgvc_id)
            elif ann_var_list[1] in ["missense_variant", "synonymous_variant"]:
                        if ann_var_list[10].split(".", 1)[0] == "p":
                            hgvs_p = ann_var_list[10]
                            variant_type = process_variant_string(hgvs_p, allele_change, variant_annotation)

                        else:
                            hgvs_id = ann_var_list[10]
                            print("%s is not captured by my code" %hgvs_id)
            if snpeFF_ann["variant_allele"] is None and snpeFF_ann["variant_position"] is None:
                snpeFF_ann["variant_allele"] = ann_var_list[0]
                snpeFF_ann["variant_position"] = data[1]
                snpeFF_ann["variant_description"].setdefault(ann_var_list[3], [variant_type, ann_var_list[1], ann_var_list[2]])
            else:
                if ann_var_list[0] != snpeFF_ann["variant_allele"]:
                    pass
                else:
                    snpeFF_ann["variant_description"].setdefault(ann_var_list[3], [variant_type, ann_var_list[1], ann_var_list[2]])

            extracted_info[data[1]]["snpeFF_ann_info"] = snpeFF_ann
    return extracted_info


def combine_similar_key_values(input_dict):
    value_keys = {}
    key_values = {}
    for key, value in input_dict.items():
        if value not in value_keys:
            value_keys.setdefault(value, [key])
        else:
            value_keys[value].append(key)

    for keyx, valuex in value_keys.items():
        if len(valuex) == 1:
            key_values[valuex[0]] = keyx
        else:
            new_keys_ranges = f"{valuex[:]}"
            key_values[new_keys_ranges] = keyx
    return key_values


def process_a_vcf(a_vcf_file, variant_caller = "bcftools"):
    extracted_info = {}
    with open(a_vcf_file, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue
            else:
                data = line.rstrip().rsplit("\t")
                if variant_caller == "bcftools":
                    extracted_info = process_bcf_vcf(data, extracted_info)
                elif variant_caller == "minos":
                    pass
    var_pos_not_in_gene_annotation = []
    var_pos_in_gene_annotation = []
    var_pos_in_lineage_info = []
    var_pos_in_amr_info = []
    var_pos_in_vr_info = []
    var_pos_in_gene_annotation = []
    for var_pos in extracted_info.keys():
         if "gene_annotation_info" not in extracted_info[var_pos].keys():
              var_pos_not_in_gene_annotation.append(var_pos)
         elif "lineage_info" in extracted_info[var_pos].keys():
            var_pos_in_lineage_info.append(var_pos)
         elif "amr_info" in extracted_info[var_pos].keys():
            var_pos_in_amr_info.append(var_pos)
         elif "vr_info" in extracted_info[var_pos].keys():
            var_pos_in_vr_info.append(var_pos)
         elif "gene_annotation_info" in extracted_info[var_pos].keys():
            var_pos_in_gene_annotation.append(var_pos)            

    
    extracted_gene_annot = {extracted_info[pos]["gene_annotation_info"]["gene"][0]:{pos: [extracted_info[data[1]]["allele_change"], extracted_info[pos]["snpeFF_ann_info"]["variant_description"][extracted_info[pos]["gene_annotation_info"]["gene"][0]][0], extracted_info[pos]["gene_annotation_info"]["CDS"]]} for pos in var_pos_in_gene_annotation
                            if "gene" in extracted_info[pos]["gene_annotation_info"].keys()
                            if extracted_info[pos]["gene_annotation_info"]["gene"][0] != "protein_coding"
                            if "CDS" in extracted_info[pos]["gene_annotation_info"].keys()
                            if extracted_info[pos]["gene_annotation_info"]["gene"][0] in extracted_info[pos]["snpeFF_ann_info"]["variant_description"].keys()
                            if extracted_info[pos]["snpeFF_ann_info"]["variant_description"][extracted_info[pos]["gene_annotation_info"]["gene"][0]][1] == "missense_variant"}
    
    
    extracted_gene_annot2 = {pos: extracted_info[pos]["gene_annotation_info"]["gene"][0] for pos in var_pos_in_gene_annotation
                             if "gene" in extracted_info[pos]["gene_annotation_info"].keys()
                             if extracted_info[pos]["gene_annotation_info"]["gene"][0] != "protein_coding"
                             if extracted_info[pos]["gene_annotation_info"]["gene"][0] in extracted_info[pos]["snpeFF_ann_info"]["variant_description"].keys()
                             if extracted_info[pos]["snpeFF_ann_info"]["variant_description"][extracted_info[pos]["gene_annotation_info"]["gene"][0]][1] == "missense_variant"}
    
    #print(combine_similar_key_values(extracted_gene_annot2))
    
    extracted_vr_region = {pos: extracted_info[pos]["vr_info"]["vr_locus_tag"] for pos in var_pos_in_vr_info
                           if "vr_info" in extracted_info[pos]}
    
    
    
    extracted_repeat_gene_annot = {pos:extracted_info[pos]["gene_annotation_info"]["repeat_region"][1] for pos in var_pos_in_gene_annotation
                                   if "repeat_region" in extracted_info[pos]["gene_annotation_info"].keys()
                                   }

    
    
    extracted_amr_annot = {pos:extracted_info[pos]["amr_info"]["antibiotics_gene"] for pos in var_pos_in_amr_info
                           if "amr_info" in extracted_info[pos].keys()}

    
    
    extracted_lineage_annot = {pos:extracted_info[pos]["lineage_info"] for pos in var_pos_in_lineage_info}

    extracted_pseudo_gene_annot = {pos:extracted_info[pos]["gene_annotation_info"]["pseudogene"][0] for pos in var_pos_in_gene_annotation
                                   if "pseudogene" in extracted_info[pos]["gene_annotation_info"]
                                   if extracted_info[pos]["gene_annotation_info"]["pseudogene"][0] != "unknown"}

    
    
    extracted_seq_feature_gene_annot = {pos:extracted_info[pos]["gene_annotation_info"]["sequence_feature"][1] for pos in var_pos_in_gene_annotation
                                        if "sequence_feature" in extracted_info[pos]["gene_annotation_info"]
                                        }
    
    extracted_mob_elem_gene_annot = {pos:extracted_info[pos]["gene_annotation_info"]["mobile_genetic_element"][1] for pos in var_pos_in_gene_annotation
                                    if "mobile_genetic_element" in extracted_info[pos]["gene_annotation_info"]}
   
    extracted_transcript_gene_annot = {pos:extracted_info[pos]["gene_annotation_info"]["transcript"][1] for pos in var_pos_in_gene_annotation
                                        if "transcript" in extracted_info[pos]["gene_annotation_info"]}

    extracted_ncRNA_gene_annot = {pos:extracted_info[pos]["gene_annotation_info"]["ncRNA"][1] for pos in var_pos_in_gene_annotation
                                        if "ncRNA" in extracted_info[pos]["gene_annotation_info"]}

    #extracted_transcript_gene_annot
    #print(combine_similar_key_values(extracted_ncRNA_gene_annot))
    
    return extracted_info



a_vcf_file = "/home/dfgmrtc/Workflows/wf-vdg/New_results/custom_annot/sample_2_bt_c_ann.vcf"

process_a_vcf(a_vcf_file)
#print(process_a_vcf(a_vcf_file))