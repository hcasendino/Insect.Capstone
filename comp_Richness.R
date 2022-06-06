###=== Characterizing Creek Communities among Bellingham Creeks - Richness======
# Written by Helen Casendino (hcas1024@uw.edu) & Ezza 
# Created: 28 Jan 2022   Modified: 5 June 2022

# Dependencies
library(tidyverse)
library(here)
library(vegan)
library(ggpubr)
library(gridExtra)

# Read In Data 
asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
asv_reads_annotated[which(asv_reads_annotated$Reach == "Up11" | asv_reads_annotated$Reach == "Up5"), "Reach"] <- "Up" # for richness, group padden up sites
asv_reads_annotated <- asv_reads_annotated %>% 
                          filter(mmyy != "Kangaroo" & Site != "Run")
# summary stats
length(which(asv_reads_annotated$class == "Insecta")) / nrow(asv_reads_annotated) # Classified Insects make up 0.17% of asv instances

asv_reads_annotated %>% group_by(class) %>% mutate(ClassReadSum = sum(nReads)) %>% ungroup() %>% 
  mutate(totalSum = sum(nReads)) %>% mutate(propReads = ClassReadSum/totalSum) %>% 
  group_by(class) %>% summarise(proportion = unique(propReads)) %>% filter(class == "Insecta") # Classified Insects make up 0.48% of total reads

###====Fig. 5: Total Insecta Richness among IBI groups by Site (5/29/22)=====

ibi_insect_richness <-  asv_reads_annotated %>% filter(class == "Insecta") %>%
  filter(species != "") 

richness_summary_pt1 <-  ibi_insect_richness %>% filter(order =="Plecoptera" |order =="Trichoptera"| order == "Odonata") %>% 
  group_by(Sample, Reach, Biological.replicate) %>% 
  mutate(ibi_richness = length(unique(species))) %>%  # unique spp per IBI order
  ungroup() %>% 
  select(mmyy, Site, Reach, Biological.replicate, order, ibi_richness) %>% 
  rename("taxon" ="order") %>% 
  mutate("group" = "Trichoptera + Plecoptera + Odonata ")

richness_summary_pt2 <-  ibi_insect_richness %>% filter(family == "Chironomidae" |order =="Ephemeroptera") %>% 
  group_by(Sample, Reach, Biological.replicate) %>% 
  mutate(ibi_richness = length(unique(species))) %>%  # unique spp per IBI order
  ungroup() %>% 
  select(mmyy, Site, Reach, Biological.replicate, family, ibi_richness) %>% 
  rename("taxon" ="family") %>% 
  mutate("group" = "Chironomidae + Baetidae")

richness_summary_pt3 <-  ibi_insect_richness %>%
  group_by(Sample, Reach, Biological.replicate) %>% 
  mutate(ibi_richness = length(unique(species))) %>%  # unique spp per order
  ungroup() %>% 
  select(mmyy, Site, Reach, Biological.replicate, family, ibi_richness) %>% 
  rename("taxon" ="family") %>% 
  mutate("group" = "All Insects")

richness_summary <- rbind(richness_summary_pt1, richness_summary_pt2, richness_summary_pt3)
  
# plot
richplot<- ggplot(richness_summary, aes(x=Site, y=ibi_richness)) + 
  geom_violin(aes(fill=Site)) + 
  facet_wrap(~group)+ 
  scale_fill_brewer(palette = "Greens", name = "group")  +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) + 
  labs(y = "Species Richness", x = "") + 
  theme(legend.position="none") + 
  scale_x_discrete( labels= c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.3,
               colour = "gray6")


ggsave(richplot, file = here("Figures", "SPP_3IBIgroups.png"), width = 9, height = 4)

####======REPEAT BUT WITH ASVs, NOT SPECIES==============

ibi_insect_richness <-  asv_reads_annotated %>% filter(class == "Insecta") 

richness_summary_pt1 <-  ibi_insect_richness %>% filter(order =="Plecoptera" |order =="Trichoptera" | order == "Odonata") %>% 
  group_by(Sample, Reach, Biological.replicate) %>% 
  mutate(ibi_richness = length(unique(Hash))) %>%  # unique spp per IBI order
  ungroup() %>% 
  select(mmyy, Site, Reach, Biological.replicate, order, ibi_richness) %>% 
  rename("taxon" ="order") %>% 
  mutate("group" = "Trichoptera + Plecoptera + Odonata")

richness_summary_pt2 <-  ibi_insect_richness %>% filter(family == "Chironomidae" |order =="Ephemeroptera") %>% 
  group_by(Sample, Reach, Biological.replicate) %>% 
  mutate(ibi_richness = length(unique(Hash))) %>%  # unique spp per IBI order
  ungroup() %>% 
  select(mmyy, Site, Reach, Biological.replicate, family, ibi_richness) %>% 
  rename("taxon" ="family") %>% 
  mutate("group" = "Chironomidae + Baetidae")

richness_summary_pt3 <-  ibi_insect_richness %>%
  group_by(Sample, Reach, Biological.replicate) %>% 
  mutate(ibi_richness = length(unique(Hash))) %>%  # unique spp per order
  ungroup() %>% 
  select(mmyy, Site, Reach, Biological.replicate, family, ibi_richness) %>% 
  rename("taxon" ="family") %>% 
  mutate("group" = "All Insects")

richness_summary <- rbind(richness_summary_pt1, richness_summary_pt2, richness_summary_pt3)

# plot
richplot<- ggplot(richness_summary, aes(x=Site, y=ibi_richness)) + 
  geom_violin(aes(fill=Site)) + 
  facet_wrap(~group)+ 
  scale_fill_brewer(palette = "Blues", name = "group")  +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) + 
  labs(y = "ASV Richness", x = "") + 
  theme(legend.position="none") + 
  scale_x_discrete( labels= c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.3,
               colour = "gray6")


ggsave(richplot, file = here("Figures", "ASV_3IBIgroups.png"), width = 9, height = 4)