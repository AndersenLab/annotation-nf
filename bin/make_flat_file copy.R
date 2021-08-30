library(data.table)
# library(hablar)
library(rebus)
library(tidyverse)
library(fuzzyjoin)
library(valr)

args <- c("BCSQ_samples_parsed.tsv", "AA_Scores", "AA_Length", "wormbase_name_key.txt", "divergent_regions_strain.bed.gz")

#Read-in input files
args <- commandArgs(trailingOnly = T)
# args[1] - parsed samples
# args[2] - AA scores
# args[3] - AA lengths
# args[4] - wbgene conversion
# args[5] - divergent

# read in parsed samples
parsed_sample_BCSQ <- data.table::fread(args[1])

#Parse BCSQ Annotation
new_parse_VCF <- function(df){
  parsed <- df %>% 
    tidyr::separate_rows(ANNOTATION, sep=",")%>%
    tidyr::separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|") 
}

parsed_ann_sample_BCSQ <- new_parse_VCF(parsed_sample_BCSQ)

#Format Amino Acid Change for AA_Score table
BCSQ_p_translate <- function(df){ #Fails for AA changes that go more than one position 
  translated <- df %>%  
    dplyr::mutate("AA" = stringr::str_extract(df$AMINO_ACID_CHANGE, "[A-Z]"))%>% #Grabs the first AA 
    dplyr::mutate("ALT_AA" = stringr::str_sub(df$AMINO_ACID_CHANGE,-1)) %>% #Grabs the ALT AA
    dplyr::mutate("AA_POS" = stringr::str_extract(df$AMINO_ACID_CHANGE, "([0-9])+")) %>% #Grabs the AA position
    tidyr::unite("REF_ALT_AA", "AA", "ALT_AA", sep= "|") #Unites AA, ALT
  return(translated)
}

translated_sample_BCSQ <- BCSQ_p_translate(parsed_ann_sample_BCSQ)  

#Add scores
AA_Scores <- data.table::fread(args[2]) #read in Amino Acid scores

strain_variant_score <- dplyr::left_join(translated_sample_BCSQ, AA_Scores, by = "REF_ALT_AA")

#Add Percent Protein

gff_AA_Length <- data.table::fread(args[3]) #read in transcript AA lengths from the gff file

strain_variant_score_pro <- dplyr::left_join(strain_variant_score, gff_AA_Length, by = "TRANSCRIPT")

strain_variant_score_pro <- strain_variant_score_pro %>%
  dplyr::mutate(AA_POS = as.numeric(AA_POS),
                AA_Length = as.numeric(AA_Length)) %>%
  dplyr::mutate("PerProtein" = (AA_POS/AA_Length)*100) %>%
  dplyr::mutate_if(is.numeric, round, digits=2) #Round percent protein to 2 decimal places

#Format Flat File for next steps

flat_file <- strain_variant_score_pro %>%
  dplyr::select(CHROM,POS,REF,ALT,CONSEQUENCE,GENE,TRANSCRIPT,BIOTYPE,STRAND,AMINO_ACID_CHANGE,DNA_CHANGE,Strains,BSCORE,GSCORE,PerProtein) %>%
  dplyr::rename("BLOSUM" = "BSCORE", "Grantham" = "GSCORE", "PERCENT_PROTEIN" = "PerProtein")


#Remove non-single AA substitutions. Make columns integers
clean_flat_file <- function(df){
  clean_flat_file <- df %>%
    dplyr::mutate(BLOSUM = ifelse(CONSEQUENCE != "missense" & CONSEQUENCE !="*missense" &
                                    CONSEQUENCE != "stop_gained" & CONSEQUENCE !="*stop_gained" & CONSEQUENCE !="start_lost&splice_region", NA, BLOSUM)) %>%
    dplyr::mutate(Grantham = ifelse(CONSEQUENCE != "missense" & CONSEQUENCE !="*missense" &
                                      CONSEQUENCE != "stop_gained" & CONSEQUENCE !="*stop_gained" & CONSEQUENCE !="start_lost&splice_region", NA, Grantham)) %>%
    dplyr::mutate(BLOSUM = as.numeric(BLOSUM),
                  Grantham = as.numeric(Grantham),
                  PERCENT_PROTEIN = as.numeric(PERCENT_PROTEIN))
    # hablar::convert(num(BLOSUM, Grantham, Percent_Protein)) #Replace with hablar conversion to numeric
  #dplyr::mutate(Grantham = sapply(clean_flat_file$Grantham, as.integer)) %>%
  #dplyr::mutate(BLOSUM = sapply(clean_flat_file$BLOSUM, as.integer)) %>%
  #dplyr::mutate(Percent_Protein = sapply(clean_flat_file$Percent_Protein, as.integer))
}

cleaned_flat_file <- clean_flat_file(flat_file)

# readr::write_tsv(cleaned_flat_file, "~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile.tsv")

##################
# add gene name
name_key <- data.table::fread(args[4]) %>%
  dplyr::select(wbgeneID, public_name)

add_gene <- cleaned_flat_file %>%
  dplyr::left_join(name_key, by = c( "GENE" = "wbgeneID")) %>%
  dplyr::rename("WORMBASE_ID" = "GENE") %>% 
  dplyr::rename("GENE" = "public_name")


# data.table::fwrite(add_gene, "~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile-gene.tsv")

###################
# impact scoring

# order_annotation <- function(df){
#   df[with(df, order(CHROM, POS)), ]
# }


impact_scoring <- function(df) {  #Score 3 or Less
  
  clean <- df %>%
    dplyr::mutate(multi_con = ifelse(grepl("&", CONSEQUENCE), T, F)) %>% # mark multi consequences
    dplyr::mutate(new_consq = CONSEQUENCE) %>%
    tidyr::separate_rows(CONSEQUENCE, sep = "&") %>%
    dplyr::mutate(CONSEQUENCE = sapply(CONSEQUENCE, impact_numeric)) %>%
    dplyr::group_by_at(setdiff(names(df), "CONSEQUENCE")) %>% # does this work?
    dplyr::mutate(VARIANT_IMPACT = max(CONSEQUENCE, na.rm = TRUE)) %>% # select the highest impact to be variant impact
    dplyr::ungroup() %>%
    dplyr::mutate(CONSEQUENCE = new_consq) %>%
    dplyr::select(-new_consq) %>%
    dplyr::distinct() %>%
    dplyr::mutate(VARIANT_IMPACT = sapply(VARIANT_IMPACT, impact_numeric_tocharacter)) %>% # convert from numeric to character
    # as.data.table() %>%
    # order_annotation() %>% #Re-order
    dplyr::mutate(VARIANT_IMPACT = sapply(VARIANT_IMPACT, as.character))

  return(clean)
}

#Functions called above
impact <- function(x) {  #Converts consequence impact valye
  dplyr::recode(x,
                "missense" ="HIGH",
                "synonymous"="LOW",
                "stop_lost"= "HIGH",
                "stop_gained"="HIGH",
                "inframe_deletion"="HIGH",
                "inframe_insertion"="HIGH",
                "frameshift"="HIGH",
                "splice_acceptor"="HIGH",
                "splice_donor"="HIGH",
                "start_lost"="HIGH",
                "splice_region"="HIGH",
                "stop_retained"="LOW",
                "5_prime_utr"="LOW",
                "3_prime_utr"="LOW",
                "non_coding"="LOW",
                "intron"="LOW",
                "intergenic"="LOW",
                "inframe_altering"="HIGH",
                "coding_sequence"="LOW",
                "feature_elongation"="LOW",
                "start_retained"="LOW" ,
                "*missense" ="HIGH",
                "*synonymous"="LOW",
                "*stop_lost"= "HIGH",
                "*stop_gained"="HIGH",
                "*inframe_deletion"="HIGH",
                "*inframe_insertion"="HIGH",
                "*frameshift"="HIGH",
                "*splice_acceptor"="HIGH",
                "*splice_donor"="HIGH",
                "*start_lost"="HIGH",
                "*splice_region"="HIGH",
                "*stop_retained"="LOW",
                "*5_prime_utr"="LOW",
                "*3_prime_utr"="LOW",
                "*non_coding"="LOW",
                "*intron"="LOW",
                "*intergenic"="LOW",
                "*inframe_altering"="HIGH",
                "*coding_sequence"="LOW",
                "*feature_elongation"="LOW",
                "*start_retained"="LOW")
  
} #Converts Consequence to Impact
impact_numeric <- function(x) {
  dplyr::recode(x,
                "missense" =3,
                "synonymous"=1,
                "stop_lost"= 4,
                "stop_gained"=4,
                "inframe_deletion"=3,
                "inframe_insertion"=3,
                "frameshift"=4,
                "splice_acceptor"=4,
                "splice_donor"=4,
                "start_lost"=4,
                "splice_region"=3,
                "stop_retained"=1,
                "5_prime_utr"=2,
                "3_prime_utr"=2,
                "non_coding"=2,
                "intron"=2,
                "intergenic"=2,
                "inframe_altering"=3,
                "coding_sequence"=2,
                "feature_elongation"=2,
                "start_retained"=1 ,
                "*missense" =3,
                "*synonymous"=1,
                "*stop_lost"= 4,
                "*stop_gained"=4,
                "*inframe_deletion"=3,
                "*inframe_insertion"=3,
                "*frameshift"=4,
                "*splice_acceptor"=4,
                "*splice_donor"=4,
                "*start_lost"=4,
                "*splice_region"=3,
                "*stop_retained"=1,
                "*5_prime_utr"=2,
                "*3_prime_utr"=2,
                "*non_coding"=2,
                "*intron"=2,
                "*intergenic"=2,
                "*inframe_altering"=3,
                "*coding_sequence"=2,
                "*feature_elongation"=2,
                "*start_retained"=1)
  
}
impact_numeric_tocharacter <- function(x){
  dplyr::recode(x,
                "4"="HIGH",
                '3'="HIGH",
                "2"="LOW",
                '1'="LOW")
}

impacted <- impact_scoring(add_gene) %>% #Perfrom Impact Scoring
  dplyr::mutate(linker = stringr::str_detect(.$CONSEQUENCE, pattern = "@" )) #Linker check for Linker impact score


impacted$VARIANT_IMPACT = if_else(impacted$linker == TRUE, "Linker", impacted$VARIANT_IMPACT ) #Add linker variant impact



table_ready <- impacted %>%
  dplyr::select(-c(linker, multi_con))

###################
# add divergent region flat

divergent <- data.table::fread(args[5], 
                         col.names = c("chrom", "start", "end", "strain"))

#Condense Overlapping Regions (If portion of regions shared by >1 strian
condensed <- divergent %>%
  dplyr::mutate(DIVERGENT = "D") # Add marker to Divergent Region

#Join the data - if a position is within a divergent region Divergent tag is added

join <- fuzzyjoin::genome_join(table_ready, condensed,  by = c(
  "CHROM" = "chrom",
  "POS" = "start",
  "POS" = "end"),
  mode = "left") %>%
  select(-chrom, -start, -end)

data.table::fwrite(join, "WI-BCSQ-flatfile.tsv" )