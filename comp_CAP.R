###=== Comparing Creek Communities among Bellingham Creeks======
# Looking at differences in insect communities across creek, month, and reach (using updated eDNA index)
# Written by Helen Casendino (hcas1024@uw.edu) & Ezza 
# Created: 28 Jan 2022   Modified: 13 Feb 2022

# Dependencies 
library(tidyverse)
library(here)
library(vegan)
library(ggpubr)
library(wesanderson)

asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
asv_reads_annotated <- asv_reads_annotated %>% filter(mmyy != "521" & mmyy != "621" & mmyy != "721") # FOR NOW, we'll remove may june and july because messed up sequencing runs

###====Data wrangling======
# first, index by max order proportion, not hash 

mean_hash_props <- asv_reads_annotated %>% 
  group_by(Sample,Biological.replicate) %>% 
  mutate (Tot = sum(nReads)) %>% 
  group_by(Sample,Biological.replicate, Hash) %>% 
  mutate(Row.prop = nReads / Tot) %>% # This creates the proportion on each biological replicate
  group_by(Sample, Hash) %>% 
  summarise (mean.prop = mean(Row.prop))  #"Averaging ratios between Biological replicates

order_index_insects <- left_join(mean_hash_props, asv_reads_annotated, by = c("Hash", "Sample")) %>% 
  select(Sample,Hash, mean.prop, class, order) %>% distinct() %>% 
  filter(class == "Insecta" & order != "" & !is.na(order)) %>% 
  group_by (order) %>%
  mutate (Colmax = max (mean.prop),
          Normalized.reads = mean.prop / Colmax) %>% 
  dplyr::select(-Colmax, -mean.prop)  # scales index vals to max in their order

# make table with asv cols (and site, mmyy, and reach cols) w/ normalized reads (for capscale)

complete_asv_covariate_tab <- left_join(order_index_insects, asv_reads_annotated, by = c("Sample", "Hash")) %>% 
  select(Site, mmyy, Reach, Hash, Normalized.reads) %>% distinct() %>% 
  pivot_wider(names_from = Hash, values_from = Normalized.reads) %>% 
  mutate_all(~replace(., is.na(.), 0))   

hashes_only <- complete_asv_covariate_tab %>% select(-c(Site, mmyy, Reach))
covariates_only <- complete_asv_covariate_tab %>% select(Site, mmyy, Reach)

# reach (lets make sure only 2 levels)
covariates_only[which(covariates_only$Reach == "Up11" | covariates_only$Reach == "Up5"), "Reach"] <- "Up"

###====Site differencess (capscale)======
# get bcds
ord <- vegdist(hashes_only, method= "bray")

# cap analysis
cap1 <- capscale(hashes_only~ Site,data=covariates_only, distance="bray")
sppscores(cap1) <- hashes_only

# plot
cap2 <- cap1[["CCA"]][["wa"]] %>%
  as.data.frame() %>%
  bind_cols(covariates_only) %>% 
  ggplot(aes(x = CAP1,
             y = CAP2)) +
    geom_point(size = 1.5) +
  geom_point(aes(color = Site), size = 3) +
 scale_color_manual(values=c("darkblue" , wes_palette(n=3, name="FantasticFox1"), "orangered3")) + 
  theme_bw() + 
geom_text(x = -0.42, y = 0.1, label="Site", size = 6)

ggsave(file = here("Figures", "CAPplot_Site.png"), width = 5, height = 4)

# taxa assigments

ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
rownames(CAPorder) <- NULL

#map on the classification of the ASVs to determine
CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
CAP.ASVid <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1))

Top50 <- CAP.ASVid %>% head(50)
write.csv(Top50, "CAPanalysis_Site_Top50DescASVs.csv")

###====temporal*Site differences (capscale)======

cap1 <- capscale(hashes_only~ Site + mmyy + Site*mmyy,data=covariates_only, distance="bray")
sppscores(cap1) <- hashes_only

# plot
cap2 <- cap1[["CCA"]][["wa"]] %>%
  as.data.frame() %>%
  bind_cols(covariates_only) %>%
  mutate(mmyy = as.factor(mmyy)) %>% 
  ggplot(aes(x = CAP1,
             y = CAP2)) +
  geom_point(size = 1.5) +
  geom_point(aes(color = Site, shape = mmyy), size = 3) +
  scale_color_manual(values=c("darkblue" , wes_palette(n=3, name="FantasticFox1"), "orangered3")) + 
  theme_bw()  + 
  geom_text(x = 0.22, y = 0.1, label="Site + mmyy + Site*mmyy", size = 4.2)
  
ggsave(file = here("Figures", "CAPplot_Site_month.png"), width = 5, height = 4)

###====temporal differences (capscale) WHY ONLY ONE CONSTRAINED AXIS?======

cap1 <- capscale(hashes_only ~ mmyy ,data=covariates_only, distance="bray")
sppscores(cap1) <- hashes_only

# plot
cap2 <- cap1[["CCA"]][["wa"]] %>%
  as.data.frame() %>%
  bind_cols(covariates_only) %>%
  mutate(mmyy = as.factor(mmyy)) %>% 
  ggplot(aes(x = CAP1,
             y = CAP1)) + # CHANGE? 
  geom_point(size = 1.5) +
  geom_point(aes(color = mmyy, shape = mmyy), size = 3) +
  scale_color_manual(values=c("darkblue" , wes_palette(n=3, name="FantasticFox1"), "orangered3")) + 
  theme_bw()  + 
  geom_text(x = 0.22, y = 0.1, label="mmyy", size = 4.2)

ggsave(file = here("Figures", "CAPplot_month.png"), width = 5, height = 4)

# taxa
ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
rownames(CAPorder) <- NULL

#map on the classification of the ASVs to determine
CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
CAP.ASVid <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1))

Top50 <- CAP.ASVid %>% head(50)
write.csv(Top50, "CAPanalysis_Month_Top50DescASVs.csv")

###====reach*site differences (capscale)======

cap1 <- capscale(hashes_only~ Site + Reach + Site*Reach ,data=covariates_only, distance="bray")
sppscores(cap1) <- hashes_only

# plot
cap2 <- cap1[["CCA"]][["wa"]] %>%
  as.data.frame() %>%
  bind_cols(covariates_only) %>%
  ggplot(aes(x = CAP1,
             y = CAP2)) +
  geom_point(size = 1.5) +
  geom_point(aes(color = Site, shape = Reach), size = 3) +
  scale_color_manual(values=c("darkblue" , wes_palette(n=3, name="FantasticFox1"), "orangered3")) + 
  theme_bw() + 
  geom_text(x = -0.25, y = 0.3, label="Site + Reach + Site*Reach", size = 4)

ggsave(file = here("Figures", "CAPplot_Site_reach.png"), width = 5, height = 4)

###====Comparison across Creeks using updated eDNA index (permanova)======

adonis(hashes_only~ Site,data=covariates_only, method="bray") 
adonis(hashes_only~ Site + mmyy,data=covariates_only, method="bray") 
adonis(hashes_only~ Site + mmyy + Site*mmyy,data=covariates_only, method="bray") 

adonis(hashes_only ~ Site + Reach , data=covariates_only, method="bray")
adonis(hashes_only ~ Site + Reach + Site*Reach, data=covariates_only, method="bray")

adonis(hashes_only ~ Reach*mmyy, data=covariates_only, method="bray")

adonis(hashes_only~ mmyy,data=covariates_only, method="bray") 






