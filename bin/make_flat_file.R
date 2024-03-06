library(data.table)
library(rebus)
library(dplyr)
library(tidyr)
library(readr)
library(fuzzyjoin)
library(IRanges)
library(stringr)

# args <- c("BCSQ_samples_parsed.tsv", "AA_Scores", "AA_Length", "wormbase_name_key.txt", "divergent_regions_strain.bed.gz")

#Read-in input files
args <- commandArgs(trailingOnly = T)
# args[1] - parsed samples
# args[2] - parsed_score
# args[3] - wbgene conversion
# args[4] - divergent
# args[5] - snpeff

parsed_sample_BCSQ <- data.table::fread(args[1])

#Add scores
parsed_score <- data.table::fread(args[2])

strain_variant_score <- dplyr::left_join(parsed_sample_BCSQ, parsed_score)

strain_variant_score <- strain_variant_score %>%
  tidyr::separate("ANNOTATION", into = c("CONSEQUENCE", "GENE", "TRANSCRIPT", "BIOTYPE", "STRAND", "AMINO_ACID_CHANGE", "DNA_CHANGE"), sep = "\\|") %>%
  dplyr::na_if("NA")


#Remove non-single AA substitutions. Make columns integers
clean_flat_file <- function(df){
  clean_flat_file <- df %>%
    dplyr::mutate(BLOSUM = ifelse(CONSEQUENCE != "missense" & CONSEQUENCE !="*missense" &
                                    CONSEQUENCE != "stop_gained" & CONSEQUENCE !="*stop_gained" & CONSEQUENCE !="start_lost&splice_region", NA, BLOSUM)) %>%
    dplyr::mutate(Grantham = ifelse(CONSEQUENCE != "missense" & CONSEQUENCE !="*missense" &
                                      CONSEQUENCE != "stop_gained" & CONSEQUENCE !="*stop_gained" & CONSEQUENCE !="start_lost&splice_region", NA, Grantham)) %>%
    # hablar::convert(num(BLOSUM, Grantham, Percent_Protein)) #Replace with hablar conversion to numeric
    dplyr::mutate(BLOSUM = as.numeric(BLOSUM),
                  Grantham = as.numeric(Grantham),
                  Percent_Protein = as.numeric(Percent_Protein))
}

cleaned_flat_file <- clean_flat_file(strain_variant_score)

# readr::write_tsv(cleaned_flat_file, "~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile.tsv")

##################
# add gene name - CE has a unique name key that we need to convert to
if(args[5] == "c_elegans") {
  name_key <- data.table::fread(args[3]) %>%
    dplyr::select(wbgene = wbgeneID, gene_name = public_name)  %>% 
    dplyr::mutate(gene_name = ifelse(str_detect(gene_name, "\\."), #Fix transcript name being used as gene name 
        str_extract(gene_name, pattern = "\\w+\\.\\d"), 
        gene_name)) %>% 
    dplyr::distinct()

  add_gene <- cleaned_flat_file %>%
    dplyr::left_join(name_key, by = c( "GENE" = "wbgene")) %>%
    dplyr::rename("WORMBASE_ID" = "GENE") %>% 
    dplyr::rename("GENE" = "gene_name")
} else {
  name_key <- data.table::fread(args[3],
              col.names = c("index", "id", "alt_id", "dupe", "cat1", "cat2"),
              na.strings = c(NA_character_, ""))%>%
    dplyr::select(id, alt_id)%>%
    #if the alt_id is NA, use the qx_id as the gene name if not use the alt_id
    dplyr::mutate(GENE = ifelse(is.na(alt_id),id, alt_id))

  add_gene <- cleaned_flat_file %>%
    dplyr::left_join(name_key %>% dplyr::rename(genie = GENE), by = c("GENE" = "id")) %>%
    dplyr::select(-alt_id, -GENE)%>%
    dplyr::rename(GENE = genie)
}

# data.table::fwrite(add_gene, "~/projects/b1059/projects/Ryan/csq/flat_file/WI.20210121.hard-filter.isotype.bcsq.20210401.pre.flatfile-gene.tsv")

###################
# impact scoring

# order_annotation <- function(df){
#   df[with(df, order(CHROM, POS)), ]
# }


impact_scoring <- function(df) {  #Score 3 or Less
  
  # high_cons <- c("missense", "stop_lost", "stop_gained", "inframe_deletion", "inframe_insertion", "frameshift", "splice_acceptor", 
  #                "splice_donor", "start_lost", "splice_region", "inframe_altering", "*missense", "*stop_lost", "*stop_gained", 
  #                "*inframe_deletion", "*inframe_insertion", "*frameshift", "*splice_acceptor", "*splice_donor", "*start_lost", 
  #                "*splice_region", "*inframe_altering")
  # low_cons <- c("synonymous", "stop_retained", "5_prime_utr", "3_prime_utr", "non_coding", "intron", "intergenic", "coding_sequence", 
  #               "feature_elongation", "start_retained", "*synonymous", "*stop_retained", "*5_prime_utr", "*3_prime_utr", "*non_coding", 
  #               "*intron", "*intergenic", "*coding_sequence", "*feature_elongation", "*start_retained")
  
  one <- c("synonymous", "stop_retained", "start_retained","*synonymous","*stop_retained","*start_retained")
  two <- c("5_prime_utr","3_prime_utr","non_coding","intron","intergenic","coding_sequence","feature_elongation",
           "*5_prime_utr","*3_prime_utr","*non_coding","*intron","*intergenic","*coding_sequence","*feature_elongation")
  three <- c("missense","inframe_deletion","inframe_insertion","splice_region","inframe_altering","*missense",
             "*inframe_deletion","*inframe_insertion","*splice_region","*inframe_altering")
  four <- c("stop_lost", "stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost",
            "*stop_lost", "*stop_gained", "*frameshift", "*splice_acceptor", "*splice_donor", "*start_lost", "*inframe_altering")
  
  clean <- df %>%
    dplyr::mutate(multi_con = ifelse(grepl("&", CONSEQUENCE), T, F)) %>% # mark multi consequences
    dplyr::mutate(new_consq = CONSEQUENCE) %>%
    tidyr::separate_rows(CONSEQUENCE, sep = "&") %>%
    # dplyr::mutate(CONSEQUENCE = sapply(CONSEQUENCE, impact_numeric),
    dplyr::mutate(CONSEQUENCE = dplyr::case_when(CONSEQUENCE %in% one ~ 1,
                                                 CONSEQUENCE %in% two ~ 2,
                                                 CONSEQUENCE %in% three ~ 3,
                                                 CONSEQUENCE %in% four ~ 4,
                                                 TRUE ~ 0)) %>% # for some reason na.rm is not working
    dplyr::group_by_at(setdiff(names(df), "CONSEQUENCE")) %>% # does this work?
    dplyr::mutate(VARIANT_IMPACT = max(CONSEQUENCE),
                  VARIANT_IMPACT = ifelse(VARIANT_IMPACT == 0, NA, VARIANT_IMPACT)) %>%
    # dplyr::mutate(VARIANT_IMPACT = max(CONSEQUENCE, na.rm = TRUE)) %>% # select the highest impact to be variant impact
    dplyr::ungroup() %>%
    dplyr::mutate(CONSEQUENCE = new_consq) %>%
    dplyr::select(-new_consq) %>%
    dplyr::distinct() %>%
    dplyr::mutate(VARIANT_IMPACT = dplyr::case_when(VARIANT_IMPACT %in% c(3, 4) ~ "HIGH",
                                                    VARIANT_IMPACT %in% c(1, 2) ~ "LOW",
                                                    TRUE ~ "NA")) %>%
    dplyr::mutate(VARIANT_IMPACT = ifelse(VARIANT_IMPACT == "NA", NA, VARIANT_IMPACT))

  return(clean)
}


impacted <- impact_scoring(add_gene) %>% #Perfrom Impact Scoring
  dplyr::mutate(linker = stringr::str_detect(.$CONSEQUENCE, pattern = "@" )) #Linker check for Linker impact score


impacted$VARIANT_IMPACT = if_else(impacted$linker == TRUE, "Linker", impacted$VARIANT_IMPACT ) #Add linker variant impact



table_ready <- impacted %>%
  dplyr::select(-c(linker, multi_con))

###################
# add divergent region flat

if(args[5] == "c_elegans") {
    divergent <- data.table::fread(args[4], 
                                   col.names = c("chrom", "start", "end"))
    
    #Condense Overlapping Regions (If portion of regions shared by >1 strian
    condensed <- divergent %>%
      dplyr::mutate(DIVERGENT = "D") # Add marker to Divergent Region
    
    #Join the data - if a position is within a divergent region Divergent tag is added
    
    test <- table_ready %>%
      dplyr::mutate(start = POS, end = POS)
    join <- fuzzyjoin::genome_join(test, condensed, by = c(
      "CHROM" = "chrom", 
      "start" = "start", 
      "end" = "end"), 
      mode = "left") %>%
      dplyr::select(-start.x, -start.y, -end.x, -end.y, -chrom)
    
    # join <- fuzzyjoin::genome_join(table_ready, condensed,  by = c(
    #   "CHROM" = "chrom",
    #   "POS" = "start",
    #   "POS" = "end"),
    #   mode = "left") %>%
    #   select(-chrom, -start, -end)
} else {
    join <- table_ready %>%
        dplyr::mutate(DIVERGENT = NA)
}


##################
# add snpeff annotation
#snpeff <- data.table::fread(args[5], sep = " ") %>%
  #dplyr::rename("CHROM" = "V1", "POS" = "V2", "ID" = "V3", "REF" = "V4", "ALT" = "V5", "QUAL" = "V6", "FILTER" = "V7", "ANNOTATION" = "V8") %>%
  #tidyr::separate_rows(ANNOTATION, sep = ",") %>%
  #tidyr::separate("ANNOTATION", into = c("ALLELE", "CONSEQUENCE", "SNPEFF_IMPACT", "GENE", "WORMBASE_ID", "FEATURE_TYPE", "TRANSCRIPT", "RANK", "HGVS.c", "HGVS.p", "cDNA_POS", "cds_POS", "PROTEIN_POS", "DISTANCE", "ERROR"), sep = "\\|") %>%
  #dplyr::select(CHROM, POS, REF, ALT, WORMBASE_ID, TRANSCRIPT, SNPEFF_IMPACT)

# combine with rest of annotation
final_flat <- join #%>%
  #dplyr::left_join(snpeff) %>%
  #dplyr::select(names(join)[names(join) != "DIVERGENT"], "SNPEFF_IMPACT", "DIVERGENT")


readr::write_tsv(final_flat, "WI-annotated-flatfile.tsv" )
