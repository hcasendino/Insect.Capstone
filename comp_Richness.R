###=== Characterizing Creek Communities among Bellingham Creeks - Richness======
# Written by Helen Casendino (hcas1024@uw.edu) & Ezza 
# Created: 28 Jan 2022   Modified: 30 Jan 2022

asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
asv_reads_annotated <- asv_reads_annotated %>% filter(mmyy != "521" & mmyy != "621" & mmyy != "721") # FOR NOW, we'll remove may june and july because messed up sequencing runs


###====Dependencies====
library(tidyverse)
library(here)
library(vegan)
library(ggpubr)

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



