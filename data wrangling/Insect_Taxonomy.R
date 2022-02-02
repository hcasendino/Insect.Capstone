###===Insect Taxonomy of Bellingham Creeks - Indexed & Non-Indexed======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 23 Dec 2021   Modified: 30 Jan 2022
###=====================

###====Dependencies====
library(tidyverse)
library(here)

#####====Read in Data ==========
COI.hash.annotated <- readRDS(here("Input/COI.all.previous.hashes.annotated.rds"))
COI_asv_reads_df <- read_csv(here("Input/COI_asv_reads_df.csv"))
COI_index_output_df <- read_csv(here("Input/COI_index_output_df.csv"))

#####==== Adding Taxonomy ==========

COI.hash.annotated <- COI.hash.annotated %>% rename("Hash" = "representative")
asv_reads_annotated <- left_join(COI_asv_reads_df, COI.hash.annotated, by = "Hash")
write_csv(asv_reads_annotated, "Input/COI_reads_taxonomy.csv") 

#####===== Adding Taxonomy to eDNA Indexed Data (March-April only) =========
## ASV Table: aggregate for march/april (average index across mo)

COI_index_MarAprOnly_df <- COI_index_output_df[-grep("^0821", COI_index_output_df$Sample, value = FALSE),] 
118388
COI_index_spring <- COI_index_MarAprOnly_df %>% 
  select(Site, Hash, Normalized.reads) %>% 
  distinct()  %>% # remove duplicate rows 
  group_by(Site, Hash) %>% summarise(mean.Normalized.reads = mean(Normalized.reads))

asv_index_annotated <- left_join(COI_index_spring, COI.hash.annotated, by = "Hash")

write_csv(asv_index_annotated, "Input/COI_MarApr_index_taxonomy.csv") 
