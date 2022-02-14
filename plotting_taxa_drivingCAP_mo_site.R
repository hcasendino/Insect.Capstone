###=== Making Plots of CAP results ======
# Make 2 10-panel plots of the top 10 taxa driving diffs among sites, and among months
# Written by Helen Casendino (hcas1024@uw.edu) 
# Created: 13 Feb 2022   Modified: 13 Feb 2022

# Dependencies 
library(tidyverse)
library(here)

siteTaxa <- read.csv("CAPanalysis_Site_Top50DescASVs.csv")
monthTaxa <- read.csv("CAPanalysis_Month_Top50DescASVs.csv")
asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
asv_reads_annotated <- asv_reads_annotated %>% filter(mmyy != "521" & mmyy != "621" & mmyy != "721") # FOR NOW, we'll remove may june and july because messed up sequencing runs

# first, index by max order proportion across bottles

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

# site plot 

top8 <- siteTaxa[1:8,"Hash"]
plotdata<-  order_index_insects %>% filter(Hash %in% top8) %>% 
  separate(Sample, into = c("mmyy", "Site"), sep="_") %>% 
  ungroup() %>% 
  add_row(Site = "1Prt") %>% 
  pivot_wider(names_from = Site, values_from = Normalized.reads, values_fill = 0) %>%
  na.omit() %>% 
  pivot_longer(!c(mmyy,Hash,class,order), names_to = "Site", values_to ="Normalized.reads" )

plotdata$Hash <- factor(plotdata$Hash, labels = c("Ephemeroptera", siteTaxa[2:8,"species"]))

ggplot(plotdata, aes(x = Site, y = Normalized.reads)) +
  geom_point(aes(color = Hash, shape = Site), size = 3)  +
  facet_wrap( ~ Hash)+
  theme_bw() + 
  labs(col = "Taxon", x = "Site", y = "eDNA Index") +   
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white")) +
  scale_colour_brewer(palette = "RdBu")   

ggsave(file = here("Figures", "8DiffTaxa_Site.png"), width = 8, height = 4)


# month plot 

top8 <- monthTaxa[1:8,"Hash"]
plotdata<-  order_index_insects %>% filter(Hash %in% top8) %>% 
  separate(Sample, into = c("mmyy", "Site"), sep="_") %>% 
  pivot_wider(names_from = mmyy, values_from = Normalized.reads, values_fill = 0) %>%
  pivot_longer(!c(Site,Hash,class,order), names_to = "mmyy", values_to ="Normalized.reads" )

plotdata$Hash <- factor(plotdata$Hash, labels = c("Ephemeroptera", monthTaxa[2:8,"species"]))

ggplot(plotdata, aes(x = mmyy, y = Normalized.reads)) +
  geom_point(aes(color = Hash, shape = mmyy), size = 3)  +
  facet_wrap( ~ Hash)+
  theme_bw() + 
  labs(col = "Taxon", x = "Month", y = "eDNA Index") +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white")) +
  scale_colour_brewer(palette = "RdYlBu")   

ggsave(file = here("Figures", "8DiffTaxa_Month.png"), width = 8, height = 4)



