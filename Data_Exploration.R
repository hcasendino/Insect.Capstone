###===Insect Characterisation of Bellingham Creeks======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 7 Nov 2021   Modified: 21 Dec 2021
###=====================

###====Dependencies====
library(here)
library(tidyverse)


#####====Read in Data - From Ezza & Eily ============
COI_asvs <- read.csv("20211120.combined.COI.ASV.table.csv")
MiFish_asvs  <- read.csv("20211120.combined.MiFish.ASV.table.csv")
all_seq <- read.csv("master_sequencing_datasheet_20211026.csv")
COI.all.previous.hashes.annotated <- readRDS("~/Desktop/Insect.Capstone/rds_files/COI.all.previous.hashes.annotated.rds")


#####==============

### COI ASV TABLE : Cleanup
COI_asvs_ID_cols <- COI_asvs %>% 
  separate(Sample_name, into = c("Locus", "MonthDay","Site","Reach","Bottle"), "[.]")

# For rows with missing values, which locus? (mostly delta & MBT, and some COI)
COI_asvs_ID_cols %>% 
  filter_all(any_vars(is.na(.))) %>% 
  group_by(Locus) %>% 
  summarise(Count = n())

# For the COI, which had missing values? (COI ones are all Kangaroo control) 
COI_asvs_ID_cols %>% 
  filter_all(any_vars(is.na(.))) %>% 
  filter(Locus == "COI")

COI_asvs_ID_cols %>% 
  mutate(index = (1:nrow(COI_asvs_ID_cols))) %>% 
  filter_all(any_vars(is.na(.)))




COI.all.previous.hashes.annotated

# COI: Which insects? 
COI_annotated<- COI.all.previous.hashes.annotated
COI_annotated %>% filter(order == "Diptera")%>% group_by(family) %>% summarise(sum = n())
COI_annotated %>% filter(order == "Ephemeroptera")%>% group_by(family) %>% summarise(sum = n())
COI_annotated %>% filter(order == "Plecoptera") %>% group_by(family) %>% summarise(sum = n())
COI_annotated %>% filter(order == "Trichoptera") %>% group_by(family) %>% summarise(sum = n())



