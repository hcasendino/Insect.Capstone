###=== Characterizing Creek Communities among Bellingham Creeks======
# Written by Helen Casendino (hcas1024@uw.edu) & Ezza 
# Created: 28 Jan 2022   Modified: 30 Jan 2022

asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
asv_reads_annotated <- asv_reads_annotated %>% filter(mmyy != "521" & mmyy != "621" & mmyy != "721") # FOR NOW, we'll remove may june and july because messed up sequencing runs

# summary stats
length(which(asv_reads_annotated$class == "Insecta")) / nrow(asv_reads_annotated) # Classified Insects make up 0.17% of asv instances

asv_reads_annotated %>% group_by(class) %>% mutate(ClassReadSum = sum(nReads)) %>% ungroup() %>% 
  mutate(totalSum = sum(nReads)) %>% mutate(propReads = ClassReadSum/totalSum) %>% 
  group_by(class) %>% summarise(proportion = unique(propReads)) %>% filter(class == "Insecta") # Classified Insects make up 0.48% of total reads
  

###====Dependencies====
library(tidyverse)
library(here)
library(vegan)
library(ggpubr)
library(gt)

###====Gross Insecta Richness (asv) across Creeks======

insect_richness <-  asv_reads_annotated %>% filter(class == "Insecta") %>%  
  group_by(Sample, Biological.replicate) %>% 
  mutate(richness = length(unique(Hash))) %>% ungroup()

# violin plot
ggplot(insect_richness, aes(x=Site, y=richness)) + 
  geom_violin(aes(fill = Site))  + 
  scale_fill_brewer(palette="RdBu") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Insect Richness", y = "ASV Richness Per Bottle") +
  theme_classic()  +  theme(plot.title = element_text(hjust = 0.5))

ggsave(file = here("Figures", "total_insect_rich_asv.png"), width = 5, height = 4)

###====Insecta Richness (asv) for IBI orders across Creeks======

# IBI orders = ephemeroptera, trichoptera, plecoptera

ephem_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & !is.na(genus) & order == "Ephemeroptera") %>% 
  group_by(Sample, Biological.replicate) %>% 
  mutate(e_richness = length(unique(genus))) %>% ungroup() 

pleco_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & !is.na(genus) & order == "Plecoptera") %>% 
  group_by(Sample, Biological.replicate) %>% 
  mutate(p_richness = length(unique(genus))) %>% ungroup() 

trich_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & !is.na(genus) & order == "Trichoptera") %>% 
  group_by(Sample, Biological.replicate) %>% 
  mutate(t_richness = length(unique(genus))) %>% ungroup()

# plots
v1 <- ggplot(ephem_richness, aes(x=Site, y=e_richness)) + 
  geom_violin(aes(fill = Site))  + 
  scale_fill_brewer(palette="Blues") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Genus Richness - Ephemeroptera", y = "Richness Per Bottle") +
  theme_classic()  +  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

v2 <- ggplot(pleco_richness, aes(x=Site, y=p_richness)) + 
  geom_violin(aes(fill = Site))  + 
  scale_fill_brewer(palette="Blues") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Genus Richness - Plecoptera", y = "Richness Per Bottle") +
  theme_classic()  +  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

v3 <- ggplot(trich_richness, aes(x=Site, y=t_richness)) + 
  geom_violin(aes(fill = Site))  + 
  scale_fill_brewer(palette="Blues") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Genus Richness - Trichoptera", y = "Richness Per Bottle") +
  theme_classic()  +  theme(plot.title = element_text(hjust = 0.5))

ggarrange(v1, v2, v3, ncol = 3, nrow = 1)
ggsave(file = here("Figures", "IBI_insect_rich_genus.png"), width = 15, height = 4)










# IBI orders = ephemeroptera, trichoptera, plecoptera

ephem_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & order == "Ephemeroptera") %>% 
  group_by(Sample, Biological.replicate) %>% 
  mutate(e_richness = length(unique(Hash))) %>% ungroup() %>% #to make plotting 0s easier bc portage has no ephem, I'll add a row
  add_row(Site = "1Prt", e_richness = 0)

pleco_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & order == "Plecoptera") %>% 
  group_by(Sample, Biological.replicate) %>% 
  mutate(p_richness = length(unique(Hash))) %>% ungroup() 

trich_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & order == "Trichoptera") %>% 
  group_by(Sample, Biological.replicate) %>% 
  mutate(t_richness = length(unique(Hash))) %>% ungroup()

# plots
v1 <- ggplot(ephem_richness, aes(x=Site, y=e_richness)) + 
  geom_violin(aes(fill = Site))  + ylim(c(0,8)) + 
  scale_fill_brewer(palette="Blues") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Ephemeroptera", y = "Richness Per Bottle") +
  theme_classic()  +  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

v2 <- ggplot(pleco_richness, aes(x=Site, y=p_richness)) + 
  geom_violin(aes(fill = Site))  + ylim(c(0,8)) + 
  scale_fill_brewer(palette="Blues") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Plecoptera", y = "Richness Per Bottle") +
  theme_classic()  +  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

v3 <- ggplot(trich_richness, aes(x=Site, y=t_richness)) + 
  geom_violin(aes(fill = Site))  + ylim(c(0,8)) + 
  scale_fill_brewer(palette="Blues") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Trichoptera", y = "Richness Per Bottle") +
  theme_classic()  +  theme(plot.title = element_text(hjust = 0.5))

ggarrange(v1, v2, v3, ncol = 3, nrow = 1)
ggsave(file = here("Figures", "IBI_insect_rich_asv.png"), width = 15, height = 4)


###====Insecta Characterization across Creeks using updated eDNA index (capscale)======

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

###====temporal differences across Creeks using updated eDNA index (capscale)======

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

###====reach differences across Creeks using updated eDNA index (capscale)======

cap1 <- capscale(hashes_only~ Site + Reach,data=covariates_only, distance="bray")
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
  geom_text(x = -0.25, y = 0.3, label="Reach", size = 6)

ggsave(file = here("Figures", "CAPplot_Site_reach.png"), width = 5, height = 4)

###====Comparison across Creeks using updated eDNA index (permanova)======

adonis(hashes_only~ Site,data=covariates_only, method="bray") 
adonis(hashes_only~ Site + mmyy,data=covariates_only, method="bray") 
adonis(hashes_only~ Site + mmyy + Site*mmyy,data=covariates_only, method="bray") 

# reach (lets make sure only 2 levels)
covariates_only[which(covariates_only$Reach == "Up11" | covariates_only$Reach == "Up5"), "Reach"] <- "Up"

adonis(hashes_only ~ Site + Reach , data=covariates_only, method="bray")
adonis(hashes_only ~ Site + Reach + Site*Reach, data=covariates_only, method="bray")

adonis(hashes_only ~ Reach*mmyy, data=covariates_only, method="bray")



###====Taxa Richness (not using asv) (OLD)===============
 
# Gross Insecta Richness (species & genus) across Creeks !!!! 

insect_sp_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & species != "" & !is.na(species)) %>%  # 0.1% of asv instances are insects and sp classified
          group_by(Sample, Biological.replicate) %>% 
          mutate(sp_richness = length(unique(species))) %>% ungroup()

insect_g_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & genus != "" & !is.na(genus)) %>%  
  group_by(Sample, Biological.replicate) %>% 
  mutate(g_richness = length(unique(genus))) %>% ungroup()

# violin plots
v1 <- ggplot(insect_sp_richness, aes(x=Site, y=sp_richness)) + 
  geom_violin(aes(fill = Site))  + 
  scale_fill_brewer(palette="RdBu") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Insect Species Richness", y = "Richness Per Bottle") +
    theme_classic()  +  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

v2 <- ggplot(insect_g_richness, aes(x=Site, y=g_richness)) + 
  geom_violin(aes(fill = Site))  + 
  scale_fill_brewer(palette="RdBu") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Insect Genus Richness", y = "Richness Per Bottle") +
  theme_classic()  +  theme(plot.title = element_text(hjust = 0.5))

ggarrange(v1, v2,  ncol = 2, nrow = 1)
ggsave(file = here("Figures", "total_insect_rich_sp_genus.png"), width = 10, height = 4)

# Insecta Richness (genus) for IBI orders across Creeks !!!! 

# IBI orders = ephemeroptera, trichoptera, plecoptera

ephem_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & genus != "" & !is.na(genus) & order == "Ephemeroptera") %>% 
  group_by(Sample, Biological.replicate) %>% 
  mutate(e_richness = length(unique(genus))) %>% ungroup() 

pleco_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & genus != "" & !is.na(genus) & order == "Plecoptera") %>% 
  group_by(Sample, Biological.replicate) %>% 
  mutate(p_richness = length(unique(genus))) %>% ungroup() 

trich_richness <-  asv_reads_annotated %>% filter(class == "Insecta" & genus != "" & !is.na(genus) & order == "Trichoptera") %>% 
  group_by(Sample, Biological.replicate) %>% 
  mutate(t_richness = length(unique(genus))) %>% ungroup()

# plots
v1 <- ggplot(ephem_richness, aes(x=Site, y=e_richness)) + 
  geom_violin(aes(fill = Site))  + ylim(c(0,8)) + 
  scale_fill_brewer(palette="Blues") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Genus Richness - Ephemeroptera", y = "Richness Per Bottle") +
  theme_classic()  +  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

v2 <- ggplot(pleco_richness, aes(x=Site, y=p_richness)) + 
  geom_violin(aes(fill = Site))  + ylim(c(0,8)) + 
  scale_fill_brewer(palette="Blues") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Genus Richness - Plecoptera", y = "Richness Per Bottle") +
  theme_classic()  +  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

v3 <- ggplot(trich_richness, aes(x=Site, y=t_richness)) + 
  geom_violin(aes(fill = Site))  + ylim(c(0,8)) + 
  scale_fill_brewer(palette="Blues") +  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  labs(title = "Genus Richness - Trichoptera", y = "Richness Per Bottle") +
  theme_classic()  +  theme(plot.title = element_text(hjust = 0.5))

ggarrange(v1, v2, v3, ncol = 3, nrow = 1)
ggsave(file = here("Figures", "IBI_insect_rich_genus.png"), width = 15, height = 4)



