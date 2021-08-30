#Libraries
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

args <- commandArgs(trailingOnly = T)
# args1 - WI-BCSQ.tsv
# args2 - AA scores
# args3 - transcript AA lengths


#Read-in Extracted Table Files
BCSQ <- data.table::fread(args[1])

#Parse Files
new_parse_VCF <- function(df){
  parsed <- df %>% 
    dplyr::rename("CHROM" = "V1", "POS" = "V2", "ID" = "V3", "REF" = "V4", "ALT" = "V5", "QUAL" = "V6", "FILTER" = "V7", "ANNOTATION" = "V8") %>%
    tidyr::separate_rows(ANNOTATION, sep=",")%>%
    tidyr::separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|") 
}

parsed_BCSQ <- new_parse_VCF(BCSQ)

#Add ref_alt column for AA substitution

BCSQ_p_translate <- function(df){ #Fails for AA changes that go more than one position 
  translated <- df %>%  
    dplyr::mutate("AA" = stringr::str_extract(df$AMINO_ACID_CHANGE, "[A-Z]"))%>% #Grabs the first AA 
    dplyr::mutate("ALT_AA" = stringr::str_sub(df$AMINO_ACID_CHANGE,-1)) %>% #Grabs the ALT AA
    dplyr::mutate("AA_POS" = stringr::str_extract(df$AMINO_ACID_CHANGE, "([0-9])+")) %>% #Grabs the AA position
    tidyr::unite("REF_ALT_AA", "AA", "ALT_AA", sep= "|") #Unites AA, ALT
  return(translated)
} #Translate BCSQ Amino Acids 

clean_BCSQ <- BCSQ_p_translate(parsed_BCSQ)

# data.table::fwrite(clean_BCSQ, "clean_BCSQ.csv") #Write out clean BCSQ

#Adding BLOSUM and Grantham Scores
# clean_BCSQ <- data.table::fread("clean_BCSQ.csv") #Read in clean BCSQ

AA_Scores <- data.table::fread(args[2]) #read in Amino Acid scores

clean_BCSQ <- dplyr::left_join(clean_BCSQ, AA_Scores, by = "REF_ALT_AA")

#Adding Protein Length and Calculating Percent Protein
gff_AA_Length <- data.table::fread(args[3]) #read in transcript AA lengths from the gff file

scored_BCSQ <- dplyr::left_join(clean_BCSQ, gff_AA_Length, by = "TRANSCRIPT")

scored_BCSQ <- scored_BCSQ %>%
  dplyr::mutate(AA_POS = as.numeric(AA_POS)) %>%
  dplyr::mutate("PerProtein" = (AA_POS/AA_Length)*100) %>%
  dplyr::mutate_if(is.numeric, round, digits=2) #Round percent protein to 2 decimal places.

#data.table::fwrite(scored_BCSQ, "scored_BCSQ.csv") #Write out for analysis

#Creating a "fake VCF" File to Annotate VCF

#scored_BCSQ <- data.table::fread("scored_BCSQ.csv")

#Grab the desired columns required for VCF format
pre_VCF <- scored_BCSQ %>%
  dplyr::group_by(CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>%
  dplyr::summarise(BLOSUM = paste(BSCORE, collapse=","), GRANTHAM = paste(GSCORE, collapse = ",") , PERCENT_PROTEIN = paste(PerProtein, collapse = ","))

#Format annotation columns to be in the 'INFO' field of a VCF
pre_VCF$BLOSUM <- paste('BLOSUM=', pre_VCF$BLOSUM, sep="")
pre_VCF$GRANTHAM <- paste('GRANTHAM=', pre_VCF$GRANTHAM, sep="")
pre_VCF$PERCENT_PROTEIN <- paste('PERCENT_PROTEIN=', pre_VCF$PERCENT_PROTEIN, sep="")

pre_VCF <- pre_VCF %>%
  tidyr::unite("INFO", "BLOSUM", "GRANTHAM","PERCENT_PROTEIN", sep = ";")

readr::write_tsv(pre_VCF,"pre_VCF.tsv")

#Read back in file to add an extra header so that it will appear to be a VCF
pre_VCF <- data.table::fread("pre_VCF.tsv", header = FALSE)

anno_vcf <- pre_VCF %>%
  dplyr::rename("##fileformat=VCFv4.2" ="V1") #Allows file to be compatible with bcftools

anno_vcf[anno_vcf=="CHROM"]<-"#CHROM" #Makes table header part of the VCF header

readr::write_tsv(anno_vcf,"anno_vcf.vcf")
