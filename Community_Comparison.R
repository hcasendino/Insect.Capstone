###===Insect Community Comparison among Bellingham Creeks======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 4 Jan 2022   Modified: 18 Jan 2022


families <- read.csv(here("Input","COI_class_family_taxonomy_df.csv"))
orders <- read.csv(here("Input","COI_class_order_taxonomy_df.csv"))
genuses <- read.csv(here("Input","COI_class_genus_taxonomy_df.csv"))

###====Dependencies====
library(tidyverse)
library(here)
library(vegan)

###=====Bray Curtis, PERMANOVA & SIMPER on all sites (no clustering)=======

# getting data into dist form, and creating separate df with sites only (choose which taxon level)

create.df.for.comparison <- function(total_taxa_index_df){
    return(total_taxa_index_df %>% 
    filter(class == "Insecta") %>% 
    select(-class) %>% 
    filter(taxon !="") %>% 
    group_by(Site, taxon) %>% 
    summarise(mean_index = mean(mean.Normalized.reads)) %>% 
    pivot_wider(id_cols = Site, names_from = taxon, values_from = mean_index) %>% 
    mutate_all(~replace(., is.na(.), 0)) %>% 
    column_to_rownames(var = "Site"))
}

mean_taxa_index <- create.df.for.comparison(families)
Sites <- data.frame(Site = rownames(mean_taxa_index))

# getting Bray Curtis Dissimilarities 

brays <- vegdist(mean_taxa_index)

# Check Multivariate Dispersions among sites (for PERMANOVA), dispersion should be mostly equal among sites

groups <- factor(Sites$Site)
mod <- betadisper(brays, groups) # Average distance to median is 0 for each site ?? 
anova(mod)

# PERMANOVA (any difference, then, which specific sites are different)

adonis2(mean_taxa_index ~ Site, data = Sites, method = "bray", permutations = 999) # messing up because there's exactly one row for each factor??

# SIMPER Analysis (b/w different sites, which families/orders driving difference)

sim <- with(Sites, simper(mean_taxa_index, Site))
output <- summary(sim) 


###=====Bray Curtis, PERMANOVA & SIMPER on all sites (with clustering by restoration status)=======

mean_taxa_index <- create.df.for.comparison(orders)
brays <- vegdist(mean_taxa_index)

# Check Multivariate Dispersion 

groups <- factor(c(rep(1,1), rep(2,4)), labels = c("restored","unrestored"))
mod <- betadisper(brays, groups) 
anova(mod) 

# PERMANOVA 

Sites_status <- Sites %>% 
  mutate(
    Status = case_when(
      Site == "1Prt" ~ "restored",
      TRUE  ~  "unrestored"
    ))


adonis2(mean_taxa_index ~ Status, data = Sites_status, method = "bray", permutations = 999)
# p = 0.6 @ family level, 0.8 at order level

# SIMPER Analysis (doesn't apply if adonis doesn't reveal significant differences)

# sim <- with(Sites_status, simper(mean_taxa_index, Status))
# output <- summary(sim) 



###=====Figs=======

# Fig 1(a-c): NMDS plot of Bray dissim at Order, Family, Genus level with site factor (no cluster)

mean_taxa_index_o <- create.df.for.comparison(orders)
mean_taxa_index_f <- create.df.for.comparison(families)
# mean_taxa_index_g <- create.df.for.comparison(genuses)

plotbrays <- metaMDS(mean_taxa_index_o, k =2) 
ordiplot(plotbrays, type = "n")
orditorp(plotbrays, display = "species", col= "red", cex = 0.5, air = 0.01, pch=".") 
orditorp(plotbrays, display = "sites", cex = .75, air = 0.01)
ggsave(file = here("Figures", "1a_order_NMDS_noncluster.png"), width = 20, height = 8)

plotbrays <- metaMDS(mean_taxa_index_f, k =2) 
ordiplot(plotbrays, type = "n")
orditorp(plotbrays, display = "species", col= "red", cex = 0.5, air = 0.1, pch=".", main = "Family")
orditorp(plotbrays, display = "sites", cex = .75, air = 0.01)
ggsave(file = here("Figures", "1b_family_NMDS_noncluster.png"), width = 20, height = 8)

plotbrays <- metaMDS(mean_taxa_index_g, k =2) 
ordiplot(plotbrays, type = "n")
orditorp(plotbrays, display = "species", col= "red", cex = 0.5, air = 0.1, pch=".", main = "Genus")
orditorp(plotbrays, display = "sites", cex = .75, air = 0.01)
ggsave(file = here("Figures", "1c_genus_NMDS_noncluster.png"), width = 20, height = 8)


# Fig 2: 

par(mfrow = c(1,2))
# hist(vegdist(mean_taxa_index_g), main = "Genus BCDs", xlab = "", ylab = "")
hist(vegdist(mean_taxa_index_f), main = "Family BCDs", xlab = "", ylab = "")
hist(vegdist(mean_taxa_index_o), main = "Order BCDs", xlab = "", ylab = "")
ggsave(file = here("Figures", "2_taxa_bcds.png"), width = 16, height = 8)

# Fig 3: Spring - August NMDS polygon plots (padden and chuckanut)


# Table 1: Summary Statistics by site for spring (orders and genus )







