#!/bin/Rscript

suppressPackageStartupMessages({
    library(rjson)
    library(tidyverse)
    library(argparser)
})


source("/home/dfgmrtc/Workflows/wf-vdg/bin/filter_dr_res_lin_1.R")


## Load json files into R
parser <- arg_parser("Parse command-line arguments")
parser <- add_argument(parser, "--files", type = "character", nargs = "*", help = "Files not correctly parsed")

arguments <- parse_args(parser)


if (!is.null(arguments$files)){
  for (item in unlist(strsplit(arguments$files, split = ","))){
    #print(item)
    if(item == "all_dr_variants.json"){
      dr_res_table <- rjson::fromJSON(file=item)
    } else if(item == "intergenic_dr_variants.json"){
      dr_res_int_table <- rjson::fromJSON(file=item)
    } 
  }
}

dr_res_full_df <- merge_dr_files(dr_res_table, dr_res_int_table)
#print(dr_res_full_df)
dr_res_short_df <- create_short_dr_profile(dr_res_table, dr_res_int_table)

write_tsv(dr_res_full_df, "detailed_drug_resistance_profile.tsv")
write_tsv(dr_res_short_df, "short_drug_resistance_profile.tsv")
