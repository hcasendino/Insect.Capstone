###===Insect Taxonomy of Bellingham Creeks======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 23 Dec 2021   Modified: 4 Jan 2022
###=====================

###====Dependencies====
library(tidyverse)
library(here)

#####====Read in Data ==========
COI.hash.annotated <- readRDS(here("Input/COI.all.previous.hashes.annotated.rds"))
COI_index_output_df <- read_csv(here("Input/COI_index_output_df.csv"))

#####=====Adding Taxonomy=========

## ASV Table: aggregate for march/april (average index across mo)

COI_index_MarAprOnly_df <- COI_index_output_df[-grep("^0821", COI_index_output_df$Sample, value = FALSE),] 

COI_index_spring <- COI_index_MarAprOnly_df %>% 
  select(Site, Hash, Normalized.reads) %>% 
  distinct()  %>% # remove duplicate rows 
  group_by(Site, Hash) %>% summarise(mean.Normalized.reads = mean(Normalized.reads))

## ADDING TAXONOMY TO ASV Table: 

# first, look at order  
assign.taxa <- function(asv_table, annotations, taxa_col = annotations$order){
  
  asv_taxa_df <- asv_table %>% mutate(class = NA, taxon = NA) 
  identified_hashes <- rep(NA, nrow(asv_taxa_df))
  
  for(i in 1:nrow(asv_taxa_df)){
     taxa_Row <- which(annotations$representative %in% asv_taxa_df$Hash[i])
     
     if(length(taxa_Row) > 0){
       asv_taxa_df$class[i] <- annotations$class[taxa_Row]
       asv_taxa_df$taxon[i] <- taxa_col[taxa_Row]
       identified_hashes[i] <- i
     }
    
     if((i %% 1000) == 0){
       print(i/nrow(asv_taxa_df))
       }
  }
   return(asv_taxa_df)
}

output_order <- assign.taxa(COI_index_spring, COI.hash.annotated, taxa_col = annotations$order) # 98% of hashes don't have order or class ID(only "") :( 
write_csv(output_order, "Input/COI_class_order_taxonomy_df.csv") 

# second, look at family  
output_family <- assign.taxa(COI_index_spring, COI.hash.annotated, taxa_col = annotations$family) 
write_csv(output_family, "Input/COI_class_family_taxonomy_df.csv") 

# genus
output_genus <- assign.taxa(COI_index_spring, COI.hash.annotated, taxa_col = annotations$genus) 
write_csv(output_family, "Input/COI_class_genus_taxonomy_df.csv") 



#####=====Exploring Insect Subset=========

# order level
insects <- output_order %>% filter(class == "Insecta") %>% 
            select(-class) %>%
            rename("order" = "taxon") %>% 
            filter(order != "") # 436 rows 

# how many orders per site? (5 - 8)
insects %>% group_by(Site) %>% 
        summarise(nOrders = length(unique(order)))
  
# family level
insects_family <- output_family %>% filter(class == "Insecta") %>% 
  select(-class) %>%
  rename("family" = "taxon") %>% 
  filter(family != "") # 331 rows 

# how many families per site? (14 - 28)
insects_family %>% group_by(Site) %>% 
  summarise(nFam = length(unique(family)))

