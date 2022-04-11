###===Insect Taxonomy of Bellingham Creeks - ASV reads======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 23 Dec 2021   Modified: 4 Apr 2022
###=====================

###====Dependencies====
library(tidyverse)
library(here)

#####==== Read in Data & Clean ==========

COI.hash.annotated <- readRDS(here("Input/COI.all.previous.hashes.annotated.rds"))
COI_asvs <- read.csv(here("Input/20220401.combined.COI.ASV.table.csv"))

### COI ASV TABLE : Cleanup
COI_asvs_ID_cols <- COI_asvs %>% 
  separate(Sample_name, into = c("Locus", "mmyy","Site","Reach","Biological.replicate"), "[.]")

# For rows with missing values, which locus? (mostly delta & MBT, and some COI)(COI ones are all Kangaroo control) 
COI_asvs_ID_cols %>% 
  filter_all(any_vars(is.na(.))) %>% 
  group_by(Locus) %>% 
  summarise(Count = n())

# Remove the rows with missing values
rows_to_remove_df <- COI_asvs_ID_cols %>% 
  mutate(index = (1:nrow(COI_asvs_ID_cols))) %>% 
  filter_all(any_vars(is.na(.)))

clean_COI_asvs <- COI_asvs_ID_cols %>% 
  slice(-rows_to_remove_df$index) %>%  # remove rows
  select(-Locus) %>% unite("Sample", c(mmyy, Site), sep = "_", remove = F) %>% 
  filter(Site != "Kangaroo") # remove kangaroo

### COI ASV TABLE : eDNA Reads (just renaming)

COI_asv_reads_df <- clean_COI_asvs

#####==== Adding Taxonomy ==========

COI.hash.annotated <- COI.hash.annotated %>% rename("Hash" = "representative")
asv_reads_annotated <- left_join(COI_asv_reads_df, COI.hash.annotated, by = "Hash")
write_csv(asv_reads_annotated, "Input/COI_reads_taxonomy.csv") 

