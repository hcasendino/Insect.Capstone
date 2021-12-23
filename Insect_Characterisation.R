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








