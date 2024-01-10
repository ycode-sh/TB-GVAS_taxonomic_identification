#!/bin/Rscript
#library(rjson)
#library(tidyverse)


## Load json files into R
#dr_res_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_dr_res.json")
#dr_res_int_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_dr_res_int.json")
#lineage_profile_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_lineage.json")

########## PROCESS RESISTANCE DATA ###############


### Extract candidates from dr_res_df and append to empty vectors

fill_res_vectors <- function(dr_res_df){
    ## Instantiate empty vectors for all variable/column names in empty dataframe
    sample_name_vec <- character()
    drug_name_vec <- character()
    conf_grading_vec <- character()
    gene_name_vec <- character()
    variant_vec <- character()
    variant_type_vec <- character()
    variant_pos_vec <- character()
    allele_change_vec <- character()
    for (sample_name in names(dr_res_df)){
        for (drug_name in names(dr_res_df[[sample_name]])){
            for (conf_grading in names(dr_res_df[[sample_name]][[drug_name]])){
                for (gene_name in names(dr_res_df[[sample_name]][[drug_name]][[conf_grading]])){
                    for (variant in names(dr_res_df[[sample_name]][[drug_name]][[conf_grading]][[gene_name]])){
                        sample_name_vec <- c(sample_name_vec, sample_name)
                        drug_name_vec <- c(drug_name_vec, drug_name)
                        conf_grading_vec <- c(conf_grading_vec, conf_grading)
                        gene_name_vec <- c(gene_name_vec, gene_name)
                        variant_vec <- c(variant_vec, variant)
                        variant_type_vec <- c(variant_type_vec, dr_res_df[[sample_name]][[drug_name]][[conf_grading]][[gene_name]][[variant]][1])
                        variant_pos_vec <- c(variant_pos_vec, dr_res_df[[sample_name]][[drug_name]][[conf_grading]][[gene_name]][[variant]][2])
                        allele_change_vec <- c(allele_change_vec, dr_res_df[[sample_name]][[drug_name]][[conf_grading]][[gene_name]][[variant]][3])
                        }
                    } 
                }
            }
        }
    ### Construct a vector of lengths each filled vector elements
    vector_length <- c(length(sample_name_vec), length(gene_name_vec), length(variant_vec), length(variant_type_vec),
                   length(allele_change_vec), length(conf_grading_vec), length(drug_name_vec))
    
    dr_res_result_list <- list(vector_length, sample_name_vec, drug_name_vec, conf_grading_vec, gene_name_vec, variant_vec, variant_type_vec, variant_pos_vec, allele_change_vec)
    return(dr_res_result_list)
    }


### Check that all lengths value in vector length are equal. If true, fill the empty dataframe with columns with corresponding vectors
fill_res_data_frame <- function(dr_res_df){
    
    ## Instantiate an empty dataframe
    empty_res_dataframe <- tibble(sample_name = character(), gene_name = character(), variant = character(),
                         variant_type = character(), variant_pos = character(), allele_change = character(), 
                         conf_grading = character(), associated_drug = character())

    dr_res_result_list <- fill_res_vectors(dr_res_df)
    if (length(unique(dr_res_result_list[[1]])) == 1){
    
        empty_res_dataframe <- add_row(empty_res_dataframe, .rows = length(dr_res_result_list[[2]]))
        empty_res_dataframe$sample_name <- dr_res_result_list[[2]]
        empty_res_dataframe$gene_name <- dr_res_result_list[[5]] 
        empty_res_dataframe$variant <- dr_res_result_list[[6]]
        empty_res_dataframe$variant_type <- dr_res_result_list[[7]]
        empty_res_dataframe$variant_pos <- dr_res_result_list[[8]]
        empty_res_dataframe$allele_change <- dr_res_result_list[[9]]
        empty_res_dataframe$conf_grading <- dr_res_result_list[[4]]
        empty_res_dataframe$associated_drug <- dr_res_result_list[[3]]
    } 
    filled_res_data_frame <- empty_res_dataframe
    return(filled_res_data_frame)

}


## dr_res_mathcing_func

any_match <- function(pattern, ass_drugs){
    grepl(pattern, ass_drugs, perl=TRUE)
}

# Step 4: Create WHO dr_res category matching function
dr_res_matching_func <- function(unnamed_pattern_list, ass_drugs, named_pattern_list){
    res_def_msg <- character()
    if(all(!unlist(lapply(unnamed_pattern_list, any_match, ass_drugs)))){
        res_def_msg_nill <- "Nill"
        return(res_def_msg_nill)
    } else {
        
        for (res_def in names(named_pattern_list)){
            if (grepl(named_pattern_list[[res_def]], ass_drugs, perl=TRUE)){
            res_def_msg <- c(res_def_msg, res_def)
            } 
        }
        
        return(res_def_msg)
    }

    
}

merge_dr_files <- function(dr_res_table, dr_res_int_table){
    # Merge and arrange drug resistance files
    dr_res_all_df <- fill_res_data_frame(dr_res_table)
    dr_res_int_df <- fill_res_data_frame(dr_res_int_table)
    dr_res_full_df <- rbind(dr_res_all_df, dr_res_int_df)
    dr_res_full_df <- arrange(dr_res_full_df, sample_name)

    return(dr_res_full_df)
}

create_dr_res_named_list <- function(dr_res_full_df){
    # Filter full res_df into a df with conf grading "Assoc_w_R" and "Assoc_w_R_Interim"
    res_only_dr_res_full_df <-  dr_res_full_df %>%
        filter(conf_grading == "Assoc_w_R" | conf_grading == "Assoc_w_R_Interim")
    # Instantiate a named_list
    sample_name_ass_drug_named_list <- list()
    # Fill up the named list
    for (no in 1:nrow(res_only_dr_res_full_df)){
        key <- res_only_dr_res_full_df$sample_name[no]
        value <- res_only_dr_res_full_df$associated_drug[no]
        if(is.null(sample_name_ass_drug_named_list[[key]])){
            sample_name_ass_drug_named_list[[key]] <- value
            
        } else {
            sample_name_ass_drug_named_list[[key]] <- c(sample_name_ass_drug_named_list[[key]], value)
            }

    }

    return(sample_name_ass_drug_named_list) 
}

uniq_sample_names <- function(dr_res_full_df){
    
    all_sample_names <- character()

    for (no in 1:nrow(dr_res_full_df)){
        all_sample_names <- c(all_sample_names, dr_res_full_df$sample_name[no])
    }

    all_sample_names <- unique(all_sample_names)

    return(all_sample_names)
}

sort_ass_drug_str <- function(ass_drugs){
    vector_ass_dr <- unlist(strsplit(ass_drugs, ","))
    sorted_indices <- order(sapply(vector_ass_dr, function(x) substr(x, 1, 1)))
    sorted_values <- vector_ass_dr[sorted_indices]
    sorted_ass_drugs <- paste(sorted_values, collapse = ",")
    return(sorted_ass_drugs)
}
curate_short_dr_profile <- function(all_sample_names, sample_name_ass_drug_named_list, named_pattern_list, unnamed_pattern_list){
    

    dr_res_short_df <- tibble(sample_name = character(), ass_res_drugs = character(), who_drug_res_profile = character())
    
    dr_res_short_df <- add_row(dr_res_short_df, .rows = length(all_sample_names))
    
    for (no in 1:length(all_sample_names)){
        if (all_sample_names[no] %in% names(sample_name_ass_drug_named_list)){ # i.e if a sample name is present in the keys of our drug res names list
            dr_res_short_df$sample_name[no] <- all_sample_names[no]
            ass_drugs <- paste(unique(sample_name_ass_drug_named_list[[all_sample_names[no]]]), collapse = ",")  # Collaps all drugs associated with that sample_name in our drug res named list to a comma separated list
            # Sort the ass_drugs string
            ass_drugs <- sort_ass_drug_str(ass_drugs)
            dr_res_short_df$ass_res_drugs[no] <- sprintf("Resistant to %s", ass_drugs) # Include the ass comma separated list into the row corresponding to the ass_drugs column
            dr_res_short_df$who_drug_res_profile[no] <- dr_res_matching_func(unnamed_pattern_list, ass_drugs, named_pattern_list)
        } else {
            dr_res_short_df$sample_name[no] <- all_sample_names[no]
            dr_res_short_df$ass_res_drugs[no] <- c("Susceptible to all drugs")
            dr_res_short_df$who_drug_res_profile[no] <- c("Drug-Susceptible")
        }
    }

    return(dr_res_short_df)
}

create_short_dr_profile <- function(dr_res_table, dr_res_int_table){
    # Merge and arrange drug resistance files
    dr_res_full_df <- merge_dr_files(dr_res_table, dr_res_int_table)
    
    # Create drug resistance named list
    sample_name_ass_drug_named_list <- create_dr_res_named_list(dr_res_full_df)

    # Create a vector of unique sample names using full_res_dataframe
    all_sample_names <- uniq_sample_names(dr_res_full_df)

    named_pattern_list <- list("Rifampicin-Resistant" = "^(?!.*INH.*RIF).*RIF", "Isoniazid-Monoresistant" = "^(?!(.*EMB.*PZA.*RIF)|(.*PZA.*RIF)|(.*RIF)).*INH", 
    "Ethambutol-Monoresistant" = ".*EMB(?!(.*INH.*RIF)|(.*INH)|(.*PZA)|(.*RIF))", "Pyrazinamide-Monoresistant" = "^(?!(.*INH.*RIF)|(.*INH)|(.*RIF)|(.*EMB)).*PZA", "Polydrug-Resistant" = "^(?!(.*INH.*RIF)|(.*INH)|(.*RIF)).*EMB.*PZA",
    "Multidrug-Resistant " = "^(?!(.*AMI)|(.*CAP)|(.*CFZ)|(.*KAN)|(.*LEV)|(.*MXF)).*INH.*RIF", "Pre-Extensively Drug-Resistant" = "^(?=.*INH.*RIF)(.*MXF|.*CFZ|.*LEV)|(.*AMI|.*CAP|.*KAN)", 
    "Extensively Drug-Resistant" = "^(?=(.*CFZ.*INH.*RIF)|(.*INH.*MXF.*RIF)|(.*INH.*LEV.*RIF))(.*AMI|.*CAP|.*KAN)")


    unnamed_pattern_list <- list("^(?!.*INH.*RIF).*RIF", "^(?!(.*EMB.*PZA.*RIF)|(.*PZA.*RIF)|(.*RIF)).*INH", ".*EMB(?!(.*INH.*RIF)|(.*INH)|(.*PZA)|(.*RIF))", "^(?!(.*INH.*RIF)|(.*INH)|(.*RIF)|(.*EMB)).*PZA", "^(?!(.*INH.*RIF)|(.*INH)|(.*RIF)).*EMB.*PZA", "^(?!(.*AMI)|(.*CAP)|(.*CFZ)|(.*KAN)|(.*LEV)|(.*MXF)).*INH.*RIF", "^(?=.*INH.*RIF)(.*MXF|.*CFZ|.*LEV)|(.*AMI|.*CAP|.*KAN)", "^(?=(.*CFZ.*INH.*RIF)|(.*INH.*MXF.*RIF)|(.*INH.*LEV.*RIF))(.*AMI|.*CAP|.*KAN)")

    #Step 4: Create a new empty dataframe that contains a short drug resistance summary of all samples, one sample per row 
    #Function is working as it should
    
    dr_res_short_df <- curate_short_dr_profile(all_sample_names, sample_name_ass_drug_named_list, named_pattern_list, unnamed_pattern_list)

    return(dr_res_short_df)
    

}



## PROCESS LINEAGE DATA
### Extract candidates from lineage_profile_table and append to empty vectors

fill_lin_vectors <- function(lineage_profile_table){
    #print(lineage_profile_table)
    lin_sample_names <- names(lineage_profile_table)
    lin_names_vec <- character()
    reference_vec <- character()
    lineage_vec <- character()
    sub_lineage_vec <- character()
    var_position_vec <- character()
    lin_allele_change_vec <- character()
    lsp_vec <- character()
    spoligotype_vec <- character()
    for (sample_name in lin_sample_names){
        for (reference in names(lineage_profile_table[[sample_name]])){
            for (lineage in names(lineage_profile_table[[sample_name]][[reference]])){
                for (sub_lineage in names(lineage_profile_table[[sample_name]][[reference]][[lineage]])){
                    for (var_pos in names(lineage_profile_table[[sample_name]][[reference]][[lineage]][[sub_lineage]])){
                            allele_change <- lineage_profile_table[[sample_name]][[reference]][[lineage]][[sub_lineage]][[var_pos]][1]
                            lsp <- lineage_profile_table[[sample_name]][[reference]][[lineage]][[sub_lineage]][[var_pos]][2]
                            spoligotype <- lineage_profile_table[[sample_name]][[reference]][[lineage]][[sub_lineage]][[var_pos]][3]
                            lin_names_vec <- c(lin_names_vec, sample_name)
                            reference_vec <- c(reference_vec, reference)
                            lineage_vec <- c(lineage_vec, lineage)
                            sub_lineage_vec <- c(sub_lineage_vec, sub_lineage)
                            var_position_vec <- c(var_position_vec, var_pos)
                            lin_allele_change_vec <- c(lin_allele_change_vec, allele_change)
                            lsp_vec <- c(lsp_vec, lsp)
                            spoligotype_vec <- c(spoligotype_vec, spoligotype)
                        
                    }
                }
            }
        }
    }


    lineage_vector_length <- c(length(lin_names_vec), length(reference_vec), length(lineage_vec), length(sub_lineage_vec), 
                           length(var_position_vec), length(lsp_vec), length(spoligotype_vec), length(lin_allele_change_vec))
    lin_result_list <- list(lineage_vector_length, lin_names_vec, reference_vec, lineage_vec, sub_lineage_vec, var_position_vec, lin_allele_change_vec, lsp_vec, spoligotype_vec)
    return(lin_result_list)
}

### Check that all lengths value in vector length are equal. If true, fill the empty dataframe with columns with corresponding vectors
fill_lin_dataframe <- function(lineage_profile_table){
    
    # Instatiate an empty dataframe
    empty_lin_dataframe <- tibble(sample_name = character(), reference = character(), lineage = character(), 
                                   sub_lineage = character(), allele_change = character(), lsp = character(), spoligotype = character(), var_position = character())

    lin_result_list <- fill_lin_vectors(lineage_profile_table)
    if (length(unique(lin_result_list[[1]])) == 1){
        empty_lin_dataframe <- add_row(empty_lin_dataframe, .rows = length(lin_result_list[[2]]))
        empty_lin_dataframe$sample_name <- lin_result_list[[2]]
        empty_lin_dataframe$reference <- lin_result_list[[3]]
        empty_lin_dataframe$lineage <- lin_result_list[[4]]
        empty_lin_dataframe$sub_lineage <- lin_result_list[[5]]
        empty_lin_dataframe$allele_change <- lin_result_list[[7]]
        empty_lin_dataframe$lsp <- lin_result_list[[8]]
        empty_lin_dataframe$spoligotype <- lin_result_list[[9]]
        empty_lin_dataframe$var_position <- lin_result_list[[6]]
    }

    filled_lin_dataframe <- empty_lin_dataframe
    
    return(filled_lin_dataframe)
}


### Use lineage_profile_table to derive lineage short df

create_short_lin_profile <- function(lineage_profile_table){
    short_lin_profile_df <- tibble(sample_name = character(), lineage_definition = character())
    sample_name_vec <- character()
    lineage_defintion_vec <- character()
    

    for (sample_name in names(lineage_profile_table)){
        bov_most_comm_sub_lin <- character()
        per_sample_sub_lin_vec <- character()
        sample_name_vec <- c(sample_name_vec, sample_name)
        for (reference in names(lineage_profile_table[[sample_name]])){
           if ("BOV_AFRI" %in% names(lineage_profile_table[[sample_name]][[reference]]) & "BOV" %in% names(lineage_profile_table[[sample_name]][[reference]])){
                bov_most_comm_sub_lin <- paste("BOV_AFRI", "BOV", collapse = "_")
           } else {
                lin_frequency_table <- table(names(lineage_profile_table[[sample_name]][[reference]]))
                most_comm_lin <- names(lin_frequency_table)[which.max(lin_frequency_table)]
                sub_lineages <- names(lineage_profile_table[[sample_name]][[reference]][[most_comm_lin]])
                ## Append sub_lineages to existing per_sample_sub_lin_vec
                per_sample_sub_lin_vec <- c(per_sample_sub_lin_vec, sub_lineages)
           }
        }
        # Create sub_lineages_table for each table and extract the most comman sub_lineage per sample
        per_sample_sub_lin_table <- table(per_sample_sub_lin_vec)
        most_comm_sub_lin <- names(per_sample_sub_lin_table)[which.max(per_sample_sub_lin_table)]
        if (length(bov_most_comm_sub_lin) == 0){ # i.e if there are no bovine lineages
            lineage_defintion_vec <- c(lineage_defintion_vec, most_comm_sub_lin)

        } else {
            lineage_defintion_vec <- c(lineage_defintion_vec,bov_most_comm_sub_lin)
            
        }

    }
    ### Fill up the short lineage dataframe
    if (length(sample_name_vec) == length(lineage_defintion_vec)){
        short_lin_profile_df <- add_row(short_lin_profile_df, .rows = length(sample_name_vec))
        short_lin_profile_df$sample_name <- sample_name_vec
        short_lin_profile_df$lineage_definition <- lineage_defintion_vec
    }

    return(short_lin_profile_df)
}
