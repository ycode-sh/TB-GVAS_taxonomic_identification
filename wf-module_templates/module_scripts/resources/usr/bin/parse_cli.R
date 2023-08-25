#!/bin/Rscript

suppressPackageStartupMessages({
    library(rjson)
    library(tidyverse)
})


source("/home/dfgmrtc/Workflows/wf-module_templates/module_scripts/resources/usr/bin/parse_JSON.R")


## Load json files into R
dr_res_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_dr_res.json")
dr_res_int_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_dr_res_int.json")
lineage_profile_table <- rjson::fromJSON(file="/home/dfgmrtc/Practice/run_folder/n_animal_lineage.json")


dr_res_full_df <- merge_dr_files(dr_res_table, dr_res_int_table)
dr_res_short_df <- create_short_dr_profile(dr_res_table, dr_res_int_table)
lineage_full_df <- fill_lin_dataframe(lineage_profile_table)
lineage_short_df <- create_short_lin_profile(lineage_profile_table)
#print(as.data.frame(create_short_dr_profile(dr_res_table, dr_res_int_table)))

