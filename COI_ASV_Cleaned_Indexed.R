###===Bellingham Creeks eDNA Data: Cleaning & Calculating eDNA Index======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 7 Nov 2021   Modified: 4 Jan 2022
###=====================

###====Dependencies====
library(tidyverse)
library(here)

#####====Read in Data - From Ezza & Eily ============
COI_asvs <- read.csv(here("Input/20211120.combined.COI.ASV.table.csv"))
# MiFish_asvs  <- read.csv(here("Input/20211120.combined.MiFish.ASV.table.csv"))
# all_seq <- read.csv(here("Input/master_sequencing_datasheet_20211026.csv"))

#####==============

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


### COI ASV TABLE : eDNA Index

source("eDNA_index_simple_func.r")
COI_index_output_df <- eDNAindex(clean_COI_asvs, reach = FALSE) 

write_csv(COI_index_output_df, "Input/COI_index_output_df.csv")
