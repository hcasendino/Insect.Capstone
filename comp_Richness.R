###=== Characterizing Creek Communities among Bellingham Creeks - Richness======
# Written by Helen Casendino (hcas1024@uw.edu) & Ezza 
# Created: 28 Jan 2022   Modified: 1 Mar 2022


# Dependencies
library(tidyverse)
library(here)
library(vegan)
library(ggpubr)
library(gridExtra)

# Read In Data 
asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
asv_reads_annotated <- asv_reads_annotated %>% filter(mmyy != "521" & mmyy != "621" & mmyy != "721") # FOR NOW, we'll remove may june and july because messed up sequencing runs

# summary stats
length(which(asv_reads_annotated$class == "Insecta")) / nrow(asv_reads_annotated) # Classified Insects make up 0.17% of asv instances

asv_reads_annotated %>% group_by(class) %>% mutate(ClassReadSum = sum(nReads)) %>% ungroup() %>% 
  mutate(totalSum = sum(nReads)) %>% mutate(propReads = ClassReadSum/totalSum) %>% 
  group_by(class) %>% summarise(proportion = unique(propReads)) %>% filter(class == "Insecta") # Classified Insects make up 0.48% of total reads


###====Fig. 1a, 1b: Gross Insecta Richness (species) across Creeks, and by reach (PLUS FIGURE FACETED BY SITE)======
# richness is collapsed by bottle, but still separated by reach

asv_reads_annotated[which(asv_reads_annotated$Reach == "Up11" | asv_reads_annotated$Reach == "Up5"), "Reach"] <- "Up" # for richness, group padden up sites

insect_richness <-  asv_reads_annotated %>% filter(class == "Insecta") %>%
                   filter(!is.na(genus)) %>% 
                  group_by(Sample, Reach) %>% 
                  mutate(genus_richness = length(unique(genus))) %>% 
                  filter(!is.na(species) & species != "") %>%  
                  mutate(sp_richness = length(unique(species))) %>% ungroup()

siteplot <- ggplot(insect_richness, aes(x=Site, y=sp_richness)) + 
  geom_boxplot(aes(fill=Site)) + 
  scale_fill_viridis_d(option = "viridis", begin = 0.2, end = 0.9)  +
    theme_bw() +
  theme(strip.background =element_rect(fill="white")) + 
  labs(title = "Insecta Richness by Site", y = "Species Richness", x = "") + 
  theme(legend.position="none") + 
  scale_x_discrete( labels= c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"))

reachplot <- ggplot(insect_richness, aes(x=Reach, y=sp_richness)) + 
  geom_boxplot(aes(fill=Reach)) + 
  #  facet_wrap( ~ Site) + 
  scale_fill_viridis_d(option = "viridis", begin = 0.4, end = 1)  +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) + 
  labs(title = "Insecta Richness by Reach and Site", y = "Species Richness", x = "") + 
  theme(legend.position="none") + 
  scale_x_discrete( labels= c("Downstream", "Upstream"))

ggarrange(siteplot, reachplot, ncol=2)
ggsave(file = here("Figures", "total_insect_rich_species_multipanel.png"), width = 8, height = 4)

###====Fig. 0a, 0b: Gross Insecta Richness (genus) across Creeks, and by reach======

siteplot <- ggplot(insect_richness, aes(x=Site, y=genus_richness)) + 
  geom_boxplot(aes(fill=Site)) + 
  scale_fill_viridis_d(option = "viridis", begin = 0.2, end = 0.9)  +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) + 
  labs(title = "Insecta Richness by Site", y = "Genus Richness", x = "") + 
  theme(legend.position="none") + 
  scale_x_discrete( labels= c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"))

reachplot <- ggplot(insect_richness, aes(x=Reach, y=genus_richness)) + 
  geom_boxplot(aes(fill=Reach)) + 
  scale_fill_viridis_d(option = "viridis", begin = 0.4, end = 1)  +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) + 
  labs(title = "Insecta Richness by Reach", y = "Genus Richness", x = "") + 
  theme(legend.position="none") + 
  scale_x_discrete( labels= c("Downstream", "Upstream"))

ggarrange(siteplot, reachplot, ncol=2)
ggsave(file = here("Figures", "total_insect_rich_genus_multipanel.png"), width = 8, height = 4)
###====Fig. 2a, 2b: Relative Insecta Richness up/downstream, by month and creek ======

reach_diff <- insect_richness %>%
  select(Sample, Reach, sp_richness) %>% 
  distinct() %>%
  pivot_wider(names_from = Reach, values_from = sp_richness) %>% 
  mutate(richness_diffs = Up - Dn) %>% select(Sample, richness_diffs)

richness_reach_diff_df <- left_join(insect_richness, reach_diff)
toplot <- richness_reach_diff_df %>% group_by(mmyy, Site) %>% 
              summarise(richness_diffs = mean(richness_diffs))  # just getting values, not averaging over any different values (bc there's only one difference / month&site)
        
toplot$Site <- factor(toplot$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
toplot$mmyy <- factor(toplot$mmyy, labels = c("March", "April", "August"))
              
  # fix this hard code (just need to show on plot a thin band when the difference b/w up and down is zero because there is a nonzero amt on both sides, not just zero)
toplot[ which(toplot$richness_diffs == 0 & toplot$Site != "5Sqm") , "richness_diffs"] <- 0.01

diffsplot <- ggplot(toplot, aes(x= Site, y=richness_diffs)) + 
        geom_bar(aes(fill = Site), stat="identity") + 
        facet_wrap( ~ mmyy) + 
   scale_fill_viridis_d(option = "viridis", begin = 0.4, end = 1)  +
  theme_bw() +  theme( strip.background =element_rect(fill="white")) + 
  labs(title = "Reach Differences in Insecta Richness", y = "Relative Species Richness (Up - Down)", x = "") + 
  scale_x_discrete(labels= rep("", 5))

ggsave(diffsplot, file = here("Figures", "relative_insect_rich_sp_month_site.png"), width = 7, height = 4)

###====Insecta Richness (asv) for IBI orders across Creeks & months======

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


