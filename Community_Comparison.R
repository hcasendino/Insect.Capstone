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
library(formattable)

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

mean_taxa_index <- create.df.for.comparison(genuses)
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

mean_taxa_index <- create.df.for.comparison(genuses)
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
# p = 0.6 @ genus level, p = 0.6 @ family level, 0.8 @ order level

# SIMPER Analysis (doesn't apply if adonis doesn't reveal significant differences)

# sim <- with(Sites_status, simper(mean_taxa_index, Status))
# output <- summary(sim) 



###=====Comparing Spring to August (Padden and Chuckanut)=======
###=====Comparing Spring to August (Padden and Chuckanut)=======
COI.hash.annotated <- readRDS(here("Input/COI.all.previous.hashes.annotated.rds"))
COI_index_output_df <- read_csv(here("Input/COI_index_output_df.csv"))
annotations <- COI.hash.annotated 

COI_index_all <- COI_index_output_df %>% 
  select(Site, mmyy, Hash, Normalized.reads) %>% 
  distinct()  %>% # remove duplicate rows 
  group_by(Site, mmyy, Hash) %>% summarise(mean.Normalized.reads = mean(Normalized.reads)) %>% filter(Site == "4Pad" | Site == "3Chk")
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

output_spring_summer <- assign.taxa(COI_index_all, annotations, taxa_col = annotations$genus)


Pad <- 
  output_spring_summer %>% 
  filter(class == "Insecta") %>% 
  select(-class) %>% 
  filter(taxon !="") %>% 
  filter(Site == "4Pad") %>% 
  select(-Site) %>% 
  group_by(mmyy, taxon) %>% 
  summarise(mean_index = mean(mean.Normalized.reads)) %>% 
  pivot_wider(id_cols = mmyy, names_from = taxon, values_from = mean_index) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  column_to_rownames(var = "mmyy")


Chk <- 
  output_spring_summer %>% 
  filter(class == "Insecta") %>% 
  select(-class) %>% 
  filter(taxon !="") %>% 
  filter(Site == "3Chk") %>% 
  select(-Site) %>% 
  group_by(mmyy, taxon) %>% 
  summarise(mean_index = mean(mean.Normalized.reads)) %>% 
  pivot_wider(id_cols = mmyy, names_from = taxon, values_from = mean_index) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  column_to_rownames(var = "mmyy"))



plotbrays <- metaMDS(Pad, k =2) 
ordiplot(plotbrays, type = "n")
orditorp(plotbrays, display = "species", col= "red", cex = 0.5, air = 0.01, pch=".") 
orditorp(plotbrays, display = "sites", cex = .75, air = 0.01)


plotbrays <- metaMDS(Chk, k =2) 
ordiplot(plotbrays, type = "n")
orditorp(plotbrays, display = "species", col= "red", cex = 0.5, air = 0.01, pch=".") 
orditorp(plotbrays, display = "sites", cex = .75, air = 0.01)
]
###=====Figs=======

# Fig 1(a-c): NMDS plot of Bray dissim at Order, Family, Genus level with site factor (no cluster)

mean_taxa_index_o <- create.df.for.comparison(orders)
mean_taxa_index_f <- create.df.for.comparison(families)
mean_taxa_index_g <- create.df.for.comparison(genuses)

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

par(mfrow = c(1,3))
hist(vegdist(mean_taxa_index_g), main = "Genus BCDs", xlab = "", ylab = "", xlim = c(0, 0.9), ylim = c(0,5))
hist(vegdist(mean_taxa_index_f), main = "Family BCDs", xlab = "",xlim = c(0, 0.9),ylab = "")
hist(vegdist(mean_taxa_index_o), main = "Order BCDs", xlab = "",xlim = c(0, 0.9), ylab = "")
ggsave(file = here("Figures", "2_taxa_bcds.png"), width = 16, height = 8)

# Fig 3: Spring - August NMDS polygon plots (padden and chuckanut)


# Table 1: Summary Statistics by site for spring (genus )

a<-   mean_taxa_index_g %>%  rownames_to_column(var = "Site") %>% 
  pivot_longer(!Site,names_to = "Genus", values_to = "mean_index")  %>% group_by(Site) %>% 
  arrange("Genus", desc(mean_index)) %>% group_by(Site) %>% 
  slice(1:10) %>% select(-mean_index) %>%  # a bunch (more than 10) have mean index value of 1...
 pivot_wider(names_from = Site, values_from = Genus)

taxa_summary <- formattable(data.frame("Portage" = a$`1Prt`[[1]], 
                   "Barnes" = a$`2Brn`[[1]],
                   "Chuckanut" = a$`3Chk`[[1]],
                   "Padden" = a$`4Pad`[[1]],
                   "Squalicum" = a$`5Sqm`[[1]]))
ggsave(file = here("Figures", "Tab1_genus_spring_summary.png"), width = 16, height = 8)


