###=== Characterizing Creek Communities among Bellingham Creeks======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 28 Jan 2022   Modified: 30 Jan 2022


asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
asv_reads_annotated <- asv_reads_annotated %>% filter(mmyy != "521") %>% # FOR NOW, we'll remove may june and july because messed up sequencing runs
                        filter(mmyy != "621" ) %>% 
                         filter(mmyy != "721" )

###====Dependencies====
library(tidyverse)
library(here)
library(vegan)

###====Gross Insecta Richness across Creeks======

asv_reads_annotated %>% group_by(Sample, Biological.replicate) %>% 
              summarise(n()) # summaarise # of species in insecta for each creek, plot 






















