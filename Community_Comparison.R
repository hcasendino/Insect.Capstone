###===Insect Community Comparison among Bellingham Creeks======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 4 Jan 2022   Modified: 4 Jan 2022
###=====================

###====Dependencies====
library(tidyverse)
library(here)
library(vegan)

###============

insect.communities.comparison <- function(tax){
    
    readcsv_call <- paste0("Input/COI_class_",tax, "_taxonomy_df.csv") # read in csv corresponding to specific taxon level

    df <- read.csv(readcsv_call, col.names = c("Site","Hash","mean.Normalized.reads", "class", paste0(tax))) %>%
      filter(class == "Insecta") %>%
      select(-c(Hash,class)) %>% 
      group_by(Site, order) %>% 
      summarise(grand_mean_index = mean(mean.Normalized.reads)) %>%  # summarise index for each taxon by site? seems bad but can workshop
      select(grand_mean_index, everything())

    df <- df[-which(df[,ncol(df)] == ""),]  # remove empty rows
                        
    df<- df %>% pivot_wider(colnames(df), names_from = Site, values_from = grand_mean_index) %>% 
            replace(is.na(.), 0)
      
    # bray curtis (use centroids?); permanova 
    dis <- vegdist(df[,-1])
    # grouping <- factor(rep(1,5), labels= colnames(df[,-1]) )
    # centroid<- betadisper(dis, group = grouping) (grouping is wrong)
    # adonis(formula = subset ~ test_groups) 

    return(df)
}

output <- insect.communities.comparison(tax = "order")
