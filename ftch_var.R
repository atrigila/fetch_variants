library(httr)
library(jsonlite)
library(xml2)
library(tidyverse)

# 1. Read the following table
table_variants_per_gene <- read.table("variants_per_gene.txt")

# 2. Load the following 

consequence_type <- c("stop_gained", "splice_donor_variant",
                        "splice_acceptor_variant", "missense_variant", 
                        "frameshift_variant") 
  
final_variants_all <- data.frame(gene_symbol = character(),
                                   rsid = character(),
                                   chrom = character(),
                                   position = numeric(),
                                   ref_allele = character(),
                                   alt_allele = character(),
                                   clinvar = character(),
                                   type = character(),
                                   cadd = character(),
                                   revel = character(),
                                   meta = character(),
                                   mutassess = character(),
                                   combo_id = character(),
                                   ancestry = character(),
                                   c_effect = character(),
                                   c_classification = character(),
                                   frequency_z = character(),
                                   frequency_h = character())
  
  server <- "https://rest.ensembl.org"
  middle <- "/vep/human/id/"
  variants <- table_variants_per_gene$id  
  end <- "?CADD=1;dbNSFP=MutationAssessor_rankscore,MetaLR_rankscore,REVEL_rankscore"

# 3. Read the following function

fetch_variant_info <- function(variant) {
  request_variant_extra_info <- GET(paste0(server, middle, variant, end), content_type("application/json"))
  variant_table <- (fromJSON(toJSON(content(request_variant_extra_info)))) 
  subtable_variants <- variant_table$transcript_consequences[[1]] 
  
  if(is.null(subtable_variants) == TRUE) {
    final_table_variant_to_add <- data.frame(gene_symbol = NA,
                                             rsid = NA,
                                             chrom = NA,
                                             position = NA,
                                             ref_allele = NA,
                                             alt_allele = NA,
                                             clinvar = NA,
                                             type = NA,
                                             cadd = NA,
                                             revel = NA,
                                             meta = NA,
                                             mutassess = NA,
                                             combo_id = NA,
                                             ancestry = NA,
                                             c_effect = NA,
                                             c_classification = NA,
                                             frequency_z = NA,
                                             frequency_h = NA)
    return(final_table_variant_to_add)
    
  } else {

    subtable_variants_filtered <- subtable_variants%>% 
      unnest(cols = c(consequence_terms)) %>%
      filter(consequence_terms %in% consequence_type)
    
    
    subtable_variants_filtered <- subtable_variants_filtered%>% 
      as_tibble(rownames=NA) %>% 
      select_if(names(.) %in% c("revel_rankscore", "metalr_rankscore", "mutationassessor_rankscore", "cadd_phred", "consequence_terms", "transcript_id", "gene_symbol"))  %>%
      mutate(id =  ifelse("variant_table$id" %in% names(.), variant_table$id, variant)) %>%
      rownames_to_column() %>% 
      relocate(gene_symbol, .after = last_col()) %>%
      group_by(rowname) %>% 
      summarise(across(1:gene_symbol, ~map_chr(., ~paste(., collapse=",")))) %>%
      unique() %>% rename(gene_id = gene_symbol)
    
    
    final_table_variant <- table_variants_per_gene %>%
      filter(id == variant) %>% full_join(subtable_variants_filtered, by = "id") %>%
      select_if(names(.) %in% c("gene_id", "id", "chrom", "position", "alleles", "clinical_significance",
                                "consequence_type", "cadd_phred", "revel_rankscore", "metalr_rankscore",
                                "mutationassessor_rankscore", "transcript_id") ) # it's possible to drop transcript id if undesired.
    
    final_table_variant_to_add <- final_table_variant %>%
      rownames_to_column() %>%
      rowwise() %>%
      mutate(gene_symbol = ifelse("gene_id" %in% names(.), gene_id, NA),
             rsid = ifelse("id" %in% names(.), id, NA),
             chrom = ifelse("chrom" %in% names(.), chrom, NA),
             position = ifelse("position" %in% names(.), position, NA),
             ref_allele = ifelse("alleles" %in% names(.), str_remove(alleles, ",.+"), NA),
             alt_allele = ifelse("alleles" %in% names(.), str_remove(alleles, paste0(str_remove(alleles, ",.+"),",")), NA),
             #   alleles = ifelse("alleles" %in% names(.), alleles, NA),
             clinvar = ifelse("clinical_significance" %in% names(.), clinical_significance, NA),
             type = ifelse("consequence_type" %in% names(.), consequence_type, NA),
             cadd = ifelse("cadd_phred" %in% names(.), cadd_phred, NA),
             revel = ifelse("revel_rankscore" %in% names(.), revel_rankscore, NA),
             meta = ifelse("metalr_rankscore" %in% names(.), metalr_rankscore, NA),
             mutassess = ifelse("mutationassessor_rankscore" %in% names(.), mutationassessor_rankscore, NA),
             combo_id = paste0(id, "_", alt_allele),
             ancestry = NA,
             c_effect = NA,
             c_classification = NA,
             frequency_z = NA,
             frequency_h = NA) %>%
      select(everything()) %>% select(gene_symbol, rsid, chrom, position,ref_allele, alt_allele,  clinvar,
                                      type, cadd, revel, meta, mutassess, combo_id, ancestry, c_effect,
                                      c_classification, frequency_z, frequency_h) %>% unique()
    
    
  return(final_table_variant_to_add)
    
    
  }
  
}
  
# 4.  Run the function in many cores
num_cores <- detectCores() - 1 

#start.time <- Sys.time()

results <- mclapply(variants[13262:49575], fetch_variant_info, mc.cores = num_cores)

#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken 

# 4. When the execution above ends successfully, concatenate results in a single file
final_table <- do.call(rbind.data.frame, results)

#write.table(final_table, "final_tabletxt", row.names = F, quote = F)



