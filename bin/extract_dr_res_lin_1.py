#!/bin/python

import re
import os

# Step 2: Prepare annotation files
# 2A
# This function is working as expected. It takes a list containing dr_resistance catalogue files
def define_drug_res_dict(drug_resistance_file):
    confidence_grading_summary_list = ["Assoc_w_R", "Assoc_w_R_Interim", "combo", "Not_assoc_w_R", "Not_assoc_w_R_Interim", "Uncertain_significance"]
    with open(drug_resistance_file, 'r') as drug_file:
        drug_res_dict = {conf_grade:{} for conf_grade in confidence_grading_summary_list}
        for line in drug_file:
            if line.startswith("#"):
                pass
            line_data = line.rstrip().rsplit("\t")
            gene_string = str(line_data[1]).split("_", 1)[0]
            var_string =  str(line_data[1]).split("_", 1)[1]
            conf_string = str(line_data[3])
            pos_string = str(line_data[2])
            if gene_string not in drug_res_dict[conf_string].keys():
                drug_res_dict[conf_string].setdefault(gene_string,{})
                if var_string not in drug_res_dict[conf_string][gene_string].keys():
                    drug_res_dict[conf_string][gene_string].setdefault(var_string, [])
                    drug_res_dict[conf_string][gene_string][var_string].append(pos_string)
                else:
                    drug_res_dict[conf_string][gene_string][var_string].append(pos_string)
            else:
                drug_res_dict[conf_string][gene_string].setdefault(var_string, [])
                drug_res_dict[conf_string][gene_string][var_string].append(pos_string)
    return drug_res_dict 

# This function is working as expected. It takes each dr catalogue files and calls the funcion above (define_drug_res_dict). Information from the 
# output of the function above becomes the value of a key:value dictionary while drug name is used as key
def make_per_drug_res_dict(drug_resistance_file_list):
    per_drug_dr_res_dict = {}
    for file in drug_resistance_file_list:
        drug_acronym = file.rsplit("sorted_")[1].rsplit(".tsv")[0]
        drug_res_dict = define_drug_res_dict(file)
        if drug_acronym not in per_drug_dr_res_dict.keys():
            per_drug_dr_res_dict.setdefault(drug_acronym, {})
            per_drug_dr_res_dict[drug_acronym] = drug_res_dict
    return per_drug_dr_res_dict


# 2B: 

# 3. Parse SnpEff-annotated VCF files

# The "conver_p_str" function will convert synonymous SNPs from HGVS format to WHO mutation catalogue format while "convert_c_str" 
# will convert cdna SNPs and INDELs from HGVS format to WHO mutation catalogue format

def convert_p_str(string):  # TO extract synonynus, missense, and "stop gained" variants
    fp = ''
    nfp = ''
    snp = ''
    vn = ''
    string = string.split(".", 1)[1].replace("*", "!")
    amino_acids_dict = {"Ala": "A", "Arg": "R", "Asn" : "N", "Asp" : "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", 
                    "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", 
                    "Tyr": "Y", "Val": "V", "del": "del", "!": "!"}
    pattern = r'^([a-zA-Z]{1,3})(\d+)([a-zA-Z]{1,3})?(!)?$'

    match = re.match(pattern, string)
    if match:
        if match.group(3) != "dup" and match.group(3) != "fs":
            fp = match.group(1)
            vn = match.group(2)
            if match.group(3) is None and match.group(4) == "!":
                sp = match.group(4)
            else:   
                sp = match.group(3)
            
            nfp = amino_acids_dict[fp]
            snp = amino_acids_dict[sp]
            
    
    return "".join([nfp, vn, snp])

def convert_c_str(string, allele_change, variant_eff): # To extract frameshift, upstream_gene and downstream_gene variants
    snp = ""
    sn_indel = ""
    indels = ""
    dup_type = ""
    dup_sn_indel = ""
    spl_str = string.split(".", 1)[1].replace("*", "")
    ref = allele_change.split("/", 1)[0].lower()
    alt = allele_change.split("/", 1)[1].lower()

    if "dup" in spl_str:  # Extract single nucleotide indel involving a duplication   
        numb = str(int(re.findall(r'\d+', spl_str)[0]) + 1) # Added 1 bcos I realized snpEff miscalculates var_pos in frameshift mutations

        if int(len(ref)) - int(len(alt)) == 1:
            dup_type = "del"
        elif int(len(ref)) - int(len(alt)) == -1:
            dup_type = "ins"
        dup_sn_indel = "_".join([numb, dup_type, "1", ref, alt])
        
    elif re.search("^[-]?[0-9]*[a-zA-Z][>][a-zA-Z]$", spl_str):  # Extract single nucleotide upstream and dowstream variants
        pattern = r'^([-]?[0-9]*)([a-zA-Z])(>)([a-zA-Z])$'
        match = re.match(pattern, spl_str)
        #print(match.group(1))
        if match:
            snp =  "".join([match.group(2).lower(), match.group(1), match.group(4).lower()])

    elif re.search("^[-]?[0-9]*(del)?(ins)?[a-zA-Z]$", spl_str):   # Extract single nucleotide Insertions or deletions (may be a frameshift variant or not)
        pattern = r'^([-]?[0-9]*)((ins)?(del)?)([a-zA-Z])$'
        
        match = re.match(pattern, spl_str)
        if match:
            if variant_eff == "frameshift_variant": # If sn_INDEL leads to a frameshift on Open Reading frame
                print(match.group(1))
                if match.group(2) == "ins":
                    modified_pos = str(int(match.group(1)) - 1) # snpeffs mistankenly miscalculates to +1 when it's frameshifts insertion so I subtracted 1
                elif match.group(2) == "del":
                    modified_pos = str(int(match.group(1)) + 1) # Added 1 bcos I realized snpEff miscalculates var_pos in frameshift mutations
                
                sn_indel = "_".join([modified_pos, match.group(2), "1", ref, alt])
            elif variant_eff == "upstream_gene_variant" and "-" in match.group(1): # snpeffs mistankenly miscalculates to -1 when it's upstream_gene_variant and it's a deletion, so I added +1 to remove the "-1" added
                modified_pos = str(int(match.group(1)) + 1)
                sn_indel = "_".join([modified_pos, match.group(2), "1", ref, alt])
            else:
                sn_indel = "_".join([match.group(1), match.group(2), "1", ref, alt])
            
    else:        # Extract mnp INDELS 
        re.search("^[-]?[0-9]*_[-]?[0-9]*(del)?(ins)?[a-zA-Z]*$", spl_str)
        spl_str1 = spl_str.split("_", 1)[1]   # Verified: The second item is picked
        spl_str1 = re.findall(r'\d+', spl_str1)[0]
        pattern = r'^([-]?[0-9]*_[-]?[0-9]*)((ins)?(del)?)([a-zA-Z]*)$'
        match = re.match(pattern, spl_str)
        if match:
            if match.group(2) == "del":
                spl_str2 = spl_str.split("l", 1)[1]
            elif match.group(2) == "ins":
                spl_str2 = spl_str.split("s", 1)[1]
            else:
                raise TypeError("Variant %s type not known" %match.group(2))
            nN = len(spl_str2)
            if variant_eff == "frameshift_variant":
                modified_pos = str(int(spl_str1))   # Verified: No need to add +1 even though it's a frameshifts mutation positions (Left it the way it is)
                indels = "_".join([modified_pos, match.group(2), str(nN), ref, alt])
            else:
                indels = "_".join([spl_str1, match.group(2), str(nN), ref, alt])

    if snp:
        return snp
    elif sn_indel:
        return sn_indel
    elif dup_sn_indel:
        return dup_sn_indel
    else:
        return indels

def convert_n_string(string):  #Extract inter and intragenic variants (only snps for now)
    string = string.split(".", 1)[1].replace("*", "")
    int_snps = ""
    if re.search("^[-]?[0-9]*[a-zA-Z][>][a-zA-Z]$", string):
        pattern = r'^([-]?[0-9]*)([a-zA-Z])(>)([a-zA-Z])$'
        match = re.match(pattern, string)
        if match:
            int_snps =  "".join([match.group(2).lower()])
    return int_snps

def modify_drug_res_dict(dr_res_variants, var, allele_change, Gene_name, variant_eff, var_pos, exp_drugs):  # Include mata-data on every variant. Make all metadata (and the variant they are associated with) into a list
    details_list = []
    details_list.append(variant_eff)
    details_list.append(var_pos)
    details_list.append(allele_change)
    details_list.append(exp_drugs)
    dr_res_variants[Gene_name].setdefault(var, [])
    dr_res_variants[Gene_name][var].append(details_list)
    
    return dr_res_variants

def process_any_string(dr_res_variants:dict, string, allele_change, Gene_name, variant_eff, var_pos, exp_drugs):
    function_level_dict = {"frameshift_variant":"run_c_string_func", "missense_variant":"run_p_string_func",
                            "synonymous_variant":"run_p_string_func", "upstream_gene_variant": "run_c_string_func",
                            "downstream_gene_variant":"run_c_string_func", "stop_gained": "run_p_string_func",
                            "intergenic_region":"run_n_string_func", "intragenic_variant":"run_n_string_func"}

    if function_level_dict[variant_eff] == "run_p_string_func": 
        var = convert_p_str(string)
        dr_res_variants = modify_drug_res_dict(dr_res_variants, var, allele_change, Gene_name, variant_eff,  var_pos, exp_drugs)
    elif function_level_dict[variant_eff] == "run_c_string_func":
        var = convert_c_str(string, allele_change, variant_eff)
        dr_res_variants = modify_drug_res_dict(dr_res_variants, var, allele_change, Gene_name, variant_eff,  var_pos, exp_drugs)
    elif function_level_dict[variant_eff] == "run_n_string_func":
        var = convert_n_string(string)
        dr_res_variants = modify_drug_res_dict(dr_res_variants, var, allele_change, Gene_name, variant_eff,  var_pos, exp_drugs)
    
    else:
        print("This varient_effect:%s is not captured by our software" %variant_eff)
    

    return dr_res_variants


def minos_vcf(data, lineage_positions_dict, dr_res_variants, dp_cov):
    if data[9].split(":")[0] == "0/0":  # Only screen/process hetrozygous and homozygous alt variants further 
        pass
    elif float(data[9].split(":")[1].replace('.', '0')) < float(dp_cov):   # Only screen variants with DP > 10
        pass
    else:
        if data[10] == "lineage_snps":
            if data[4] == data[17]:
                lineage_positions_dict.setdefault(data[1], "/".join([data[3], data[4]]))
        elif data[10] == "amr_regions":
            for item in range(len(data[7].split("=")[1].split(",")[:])):
                allele_change = "/".join([data[3], data[4]])
                ANN = data[7].split("=")[1].split(",")[item].split("|")[:]
                if ANN[3] in data[14]:
                    dr_res_variants.setdefault(ANN[3], {})
                    if ANN[1] in ["frameshift_variant", "upstream_gene_variant", "downstream_gene_variant", "intergenic_region", "intragenic_variant"]:
                        if ANN[9].split(".", 1)[0] == "c":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[9], allele_change, ANN[3], ANN[1], data[1], data[16])
                        elif ANN[9].split(".", 1)[0] == "n":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[9], allele_change, ANN[3], ANN[1], data[1], data[16])
                        else:
                            print("%s is not captured by c and n string code" %{ANN[9].split(".", 1)[0]})
                    elif ANN[1] in ["missense_variant", "synonymous_variant", "stop_gained"]:
                        if ANN[10].split(".", 1)[0] == "p":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[10], allele_change, ANN[3], ANN[1], data[1], data[16])
                        else:
                            string_id = ANN[10].split(".", 1)[0]
                            print("%s is not captured by my code" %string_id)


def bcftools_vcf(data, lineage_positions_dict, dr_res_variants, dp_cov):
    if data[9].split(":")[0] == "0":
        pass
    elif float(data[7].split("DP=")[1].split(";")[0]) < float(dp_cov):
        pass
    else:
        if data[10] == "lineage_snps":
            # if "/".join([data[3], data[4]]) == "/".join([data[16], data[17]]): This was deprecated because of new lineage additions in order to ensure continuity
            if data[4] == data[17]:
                lineage_positions_dict.setdefault(data[1], "/".join([data[3], data[4]]))
        elif data[10] == "amr_regions":
            for item in range(len(data[7].split("ANN=")[1].split(",")[:])):
                allele_change = "/".join([data[3], data[4]])
                ANN = data[7].split("ANN=")[1].split(",")[item].split("|")[:]
                if ANN[3] in data[14]:
                    dr_res_variants.setdefault(ANN[3], {})
                    if ANN[1] in ["frameshift_variant", "upstream_gene_variant", "downstream_gene_variant", "intergenic_region", "intragenic_variant"]:
                        if ANN[9].split(".", 1)[0] == "c":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[9], allele_change, ANN[3], ANN[1], data[1], data[16])
                        elif ANN[9].split(".", 1)[0] == "n":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[9], allele_change, ANN[3], ANN[1], data[1], data[16])
                        else:
                            print("%s is not captured by c and n string code" %{ANN[9].split(".", 1)[0]})
                    elif ANN[1] in ["missense_variant", "synonymous_variant", "stop_gained"]:
                        if ANN[10].split(".", 1)[0] == "p":
                            dr_res_variants = process_any_string(dr_res_variants, ANN[10], allele_change, ANN[3], ANN[1], data[1], data[16])
                        else:
                            string_id = ANN[10].split(".", 1)[0]
                            print("%s is not captured by my code" %string_id)

def proces_a_vcf_file(a_vcf_file, variant_caller, dp_cov):    
    lineage_positions_dict = {}
    dr_res_variants = {}
    with open(a_vcf_file, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue
            data = line.rstrip().rsplit("\t")

            if variant_caller == "bcftools":
                bcftools_vcf(data, lineage_positions_dict, dr_res_variants, dp_cov)
            elif variant_caller == "minos":
                minos_vcf(data, lineage_positions_dict, dr_res_variants, dp_cov)
                            
    return lineage_positions_dict, dr_res_variants
                            
                        
def process_many_vcf_files(vcf_file_list, variant_caller, dp_cov):
    per_file_dr_res_dict = {}
    per_file_lineage_dict = {}
    for a_vcf_file in vcf_file_list:
        vcf_file_name = os.path.basename(a_vcf_file.split("_intersect.vcf")[0])
        lineage_positions_list, dr_res_variants =  proces_a_vcf_file(a_vcf_file, variant_caller, dp_cov)
        per_file_dr_res_dict.setdefault(vcf_file_name, {})
        per_file_dr_res_dict[vcf_file_name] = dr_res_variants
        per_file_lineage_dict.setdefault(vcf_file_name, {})
        per_file_lineage_dict[vcf_file_name] = lineage_positions_list
    return per_file_dr_res_dict, per_file_lineage_dict
