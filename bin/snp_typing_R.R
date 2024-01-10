#!/bin/Rscript
suppressPackageStartupMessages({
    library(rjson)
    library(tidyverse)
    library(argparser)
})

source("/home/dfgmrtc/Workflows/wf-vdg/bin/filter_dr_res_lin_1.R")


## Load json files into R
parser <- arg_parser("Parse command-line arguments")
parser <- add_argument(parser, "--files", type = "character", nargs = "*", help = "json files to parse")

arguments <- parse_args(parser)


if (!is.null(arguments$files)){
  for (item in arguments$files){
    #print(item)
    if(item == "lineage.json"){
      lineage_profile_table <- rjson::fromJSON(file=item)
    } 
  }
}

lineage_full_df <- fill_lin_dataframe(lineage_profile_table)
lineage_short_df <- create_short_lin_profile(lineage_profile_table)



write_tsv(lineage_full_df, "detailed_lineage_assignments.tsv")
