###===Insect Characterisation of Bellingham Creeks======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 23 Dec 2021   Modified: 23 Dec 2021
###=====================

###====Dependencies====
library(tidyverse)
library(here)

#####====Read in Data ==========
COI.hash.annotated <- readRDS(here("Input/COI.all.previous.hashes.annotated.rds"))
COI_index_output_df <- read_csv(here("Input/COI_index_output_df.csv"))

#####==============

## ASV Table: aggregate for march/april (average index across mo)

COI_index_MarAprOnly_df <- COI_index_output_df[-grep("^0821", COI_index_output_df$Sample, value = FALSE),] 

COI_index_spring <- COI_index_MarAprOnly_df %>% 
  select(Site, Hash, Normalized.reads) %>% 
  distinct()  %>% # remove duplicate rows 
  group_by(Site, Hash) %>% summarise(mean.Normalized.reads = mean(Normalized.reads))

## ADDING TAXONOMY TO ASV Table: 

assign.taxa <- function(asv_table, annotations, taxa_col = annotations$order){
  
  asv_taxa_df <- asv_table %>% mutate(class = NA, taxon = NA) 
  
  for(i in 1:nrow(asv_taxa_df)){
     taxa_Row <- which(annotations$representative %in% asv_taxa_df$Hash[i])
     
     if(length(taxa_Row) > 0){
       asv_taxa_df$class[i] <- annotations$class[taxa_Row]
       asv_taxa_df$taxon[i] <- taxa_col[taxa_Row]
     }
     else{print(c(i,"unidentified Hash"))}
   }
   return(asv_taxa_df)
}

assign.taxa(COI_index_spring, COI.hash.annotated)


