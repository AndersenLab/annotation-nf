library(dplyr)
library(tidyr)
library(data.table)
library(readr)

#Read-in input files
args <- commandArgs(trailingOnly = T)

scores <- data.table::fread(args[1])

parse_score <- function(df){
  parsed <- df %>%
    dplyr::rename("CHROM" = "V1", "POS" = "V2", "REF" = "V3", "ALT" = "V4", "ANNOTATION" = "V5", "BLOSUM" = "V6", "Grantham" = "V7", "Percent_Protein" = "V8") %>%
    tidyr::separate_rows(ANNOTATION, BLOSUM, Grantham, Percent_Protein, sep=",")
  
}

scored_with_annotation <- parse_score(scores)

readr::write_tsv(scored_with_annotation, "BCSQ_scores_parsed.tsv")