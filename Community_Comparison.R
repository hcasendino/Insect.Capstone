###===Insect Community Comparison among Bellingham Creeks======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 4 Jan 2022   Modified: 18 Jan 2022
###=====================

families <- read.csv(here("Input","COI_class_family_taxonomy_df.csv"))
orders <- read.csv(here("Input","COI_class_order_taxonomy_df.csv"))

###====Dependencies====
library(tidyverse)
library(here)
library(vegan)

###=====Bray Curtis, PERMANOVA & SIMPER=======

# getting data into dist form, and creating separate df with sites only 
mean_taxa_index <- families %>% 
  filter(class == "Insecta") %>% 
  select(-class) %>% 
  filter(taxon !="") %>% 
  group_by(Site, taxon) %>% 
  summarise(mean_index = mean(mean.Normalized.reads)) %>% 
  pivot_wider(id_cols = Site, names_from = taxon, values_from = mean_index) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  column_to_rownames(var = "Site")
  
Sites <- data.frame(Site = unique(families$Site))

# getting Bray Curtis Dissimilarities 

brays <- vegdist(mean_taxa_index)
hist(brays)
 
# Check Multivariate Dispersions among sites (for PERMANOVA) 
# dispersion should be mostly equal among sites

groups <- factor(Sites$Site)
mod <- betadisper(brays, groups)
anova(mod)
permutest(mod, pairwise= T, permutations = 99)

# PERMANOVA (any difference, then, which specific sites are different)
# are adonis and betadisper failing because I haven't clustered at all? 

adonis(mean_taxa_index ~ Site, data = Sites, method = "bray", permutations = 999)

# SIMPER Analysis (b/w different sites (or site clusters), which families/orders driving difference)

sim <- with(Sites, simper(mean_taxa_index, Site))
summary(sim)
 



###=====Figs=======

# NMDS plot

plotbrays <- metaMDS(mean_taxa_index, k =2) # tress is (nearly) zero: you may have insufficient data
ordiplot(plotbrays, type = "n")
orditorp(plotbrays, display = "species", col= "red", cex = 0.5, air = 0.1, pch=".")
orditorp(plotbrays, display = "sites", cex = 1, air = 0.25)



