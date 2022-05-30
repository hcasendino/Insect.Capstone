###=== Characterizing Creek Communities among Bellingham Creeks - Richness======
# Written by Helen Casendino (hcas1024@uw.edu) & Ezza 
# Created: 28 Jan 2022   Modified: 29 May 2022


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


###====Fig. 1a, 1b: Gross Insecta Richness (species) across Creeks, and by reach======
# richness is collapsed by bottle, but still separated by reach

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
  labs( title = "Insecta Richness by Site and Reach", y = "Species Richness", x = "") + 
  theme(legend.position="none") + 
  scale_x_discrete( labels= c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"))

reachplot <- ggplot(insect_richness, aes(x=Reach, y=sp_richness)) + 
  geom_boxplot(aes(fill=Reach)) + 
  #  facet_wrap( ~ Site) + 
  scale_fill_viridis_d(option = "viridis", begin = 0.4, end = 1)  +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) + 
  labs( y = "Species Richness", x = "") + 
  theme(legend.position="none") + 
  scale_x_discrete( labels= c("Downstream", "Upstream"))

ggarrange(siteplot, reachplot, ncol=2)
#ggsave(file = here("Figures", "total_insect_rich_species_multipanel.png"), width = 8, height = 4)

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

#ggarrange(siteplot, reachplot, ncol=2)
#ggsave(file = here("Figures", "total_insect_rich_genus_multipanel.png"), width = 8, height = 4)

###====Fig. 2: Relative Insecta Richness up/downstream (sp), by month and creek ======

reach_diff <- insect_richness %>%
  select(Sample, Reach, sp_richness) %>% 
  distinct() %>%
  pivot_wider(names_from = Reach, values_from = sp_richness) %>% 
  mutate(richness_diffs = Up - Dn) %>%
  mutate(new_richness_diffs = case_when(richness_diffs == 0 & Up != 0 ~ 0.01,  #need to show on plot a thin band when the difference b/w up and down is zero because there is a nonzero amt on both sides, not just zero)
                                        TRUE ~ as.numeric(richness_diffs))) %>% ungroup() %>%
  select(Sample, new_richness_diffs)

reach_diff <- reach_diff %>% rename("richness_diffs" = new_richness_diffs )

richness_reach_diff_df <- left_join(insect_richness, reach_diff)
toplot <- richness_reach_diff_df %>% group_by(mmyy, Site) %>% 
              summarise(richness_diffs = mean(richness_diffs))  # just getting values, not averaging over any different values (bc there's only one difference / month&site)
        
toplot$Site <- factor(toplot$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
toplot$mmyy <- factor(toplot$mmyy, labels = c("March", "April", "August"))

diffsplot <- ggplot(toplot, aes(x= Site, y=richness_diffs)) + 
        geom_bar(aes(fill = Site), stat="identity") + 
        facet_wrap( ~ mmyy) + 
   scale_fill_viridis_d(option = "viridis", begin = 0.4, end = 1)  +
  theme_bw() +  theme( strip.background =element_rect(fill="white")) + 
  labs(title = "Reach Differences in Insecta Richness", y = "Relative Species Richness (Up - Down)", x = "") + 
  scale_x_discrete(labels= rep("", 5))

#ggsave(diffsplot, file = here("Figures", "relative_insect_rich_sp_month_site.png"), width = 7, height = 4)

###====Fig. 3: Insecta Richness (genus) for 3 IBI orders across Creeks & months======

ibi_insect_richness <-  asv_reads_annotated %>% filter(class == "Insecta") %>%
                filter(!is.na(genus))  #%>% 
              # filter(!is.na(species) & species != "")

ibi_insect_richness <-  ibi_insect_richness %>% filter(order == "Ephemeroptera" | order =="Plecoptera" |order =="Trichoptera") %>% 
      group_by(Sample, Reach, order) %>% 
      mutate(ibi_richness = length(unique(genus))) %>%  # unique genuses per IBI order
      ungroup() %>% 
      select(mmyy, Site, Reach, order, ibi_richness) %>%
      mutate(Site = case_when(Site == "2Brn" | Site == "3Chk" | Site == "5Sqm" ~ "unrestored",
                              Site == "4Pad" ~ "Padden",
                              Site == "1Prt" ~ "Portage")) %>% distinct()
ibi_insect_richness$mmyy <- factor(ibi_insect_richness$mmyy, labels = c("March", "April", "May", "June", "July", "August"))
ibi_insect_richness$order <- factor(ibi_insect_richness$order, levels = c("Ephemeroptera", "Trichoptera", "Plecoptera"),  labels = c("Ephemeroptera", "Trichoptera", "Plecoptera") )

# plot 
ibiplot <- ibi_insect_richness %>%
        ggplot(aes(x= Site, y=ibi_richness)) + 
        geom_boxplot(aes(fill = Site)) + 
        facet_wrap( ~ order + mmyy) + 
        scale_fill_viridis_d(option = "viridis", begin = 0.4, end = 1)  +
        theme_minimal() +  theme( strip.background =element_rect(fill="white")) + 
        labs(title = "IBI Genus Richness by Month and Site", y = "Insecta Genus Richness", x = "") 

#ggsave(file = here("Figures", "IBI_insect_rich_genus.png"), width = 8, height = 7)

###====Fig. 4: Relative Insecta Richness up/downstream (genus) for 3 IBI orders across Creeks & months======

ibi_insect_richness <-  asv_reads_annotated %>% filter(class == "Insecta") %>%
        filter(!is.na(genus)) %>% 
         filter(order == "Ephemeroptera" | order =="Plecoptera" |order =="Trichoptera") %>% 
          group_by(Sample, Reach, order) %>% 
         mutate(ibi_richness = length(unique(genus))) %>%
         select(mmyy, Site, Reach, order, ibi_richness) %>% 
           distinct() %>%
         pivot_wider(names_from = Reach, values_from = ibi_richness, values_fill = 0) %>% 
         mutate(richness_diffs = Up - Dn) %>%  # Get difference in richness by order, month, site
        mutate(new_richness_diffs = case_when(richness_diffs == 0 & Up != 0 ~ 0.1,  #need to show on plot a thin band when the difference b/w up and down is zero because there is a nonzero amt on both sides, not just zero)
                                          TRUE ~ as.numeric(richness_diffs)))  %>% ungroup() %>% 
         select(!c(Sample, Dn, Up, richness_diffs)) %>% distinct()

ibi_insect_richness$Site <- factor(ibi_insect_richness$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
ibi_insect_richness$mmyy <- factor(ibi_insect_richness$mmyy, labels = c("March", "April", "May", "June", "July", "August"))

diffsplot2 <- ggplot(ibi_insect_richness, aes(x= Site, y=new_richness_diffs)) + 
  geom_bar(aes(fill = Site), stat="identity") + 
  facet_wrap( ~ order + mmyy) + 
  scale_fill_viridis_d(option = "viridis", begin = 0.4, end = 1)  +
  theme_minimal() +  theme( strip.background =element_rect(fill="white")) + 
  labs(title = "Reach Differences in IBI Richness", x = "Site", y = "Relative Genus Richness (Up - Down)", x = "")

# ggsave(diffsplot2, file = here("Figures", "relative_ibi_rich_genus_month_site.png"), width = 10, height = 6)

diffsplot1<- ggplot(ibi_insect_richness, aes(x= Site, y=new_richness_diffs)) + 
  geom_bar(aes(fill = Site), stat="identity") + 
  facet_wrap( ~ mmyy) + 
  scale_fill_viridis_d(option = "viridis", begin = 0.4, end = 1)  +
  theme_minimal() +  theme( strip.background =element_rect(fill="white")) + 
  labs(title = "Reach Differences in IBI Richness", x = "Site", y = "Relative Genus Richness (Up - Down)", x = "")
#ggsave(diffsplot1, file = here("Figures", "relative_insect_rich_sp_month_site.png"), width = 10, height = 6)




# 3/30 proportion of asvs within order that are classified to genus or sp
step1 <- asv_reads_annotated %>%
  select(Hash, order, family, genus, species) %>% 
  filter(order == "Ephemeroptera") %>% distinct() %>% 
  mutate(classified_as = case_when(genus != "" & species == ""  ~ "genus",
                                   genus != "" & species != ""  ~ "species",
                                   genus == "" & species == ""  ~ "neither")) 

step1 %>% group_by(classified_as) %>% 
  summarise(n = n()) %>% mutate(prop_classified = n/ nrow(step1))

###====Fig. 5: Total Insecta Richness among IBI groups by Site (5/29/22)=====

ibi_insect_richness <-  asv_reads_annotated %>% filter(class == "Insecta") %>%
  filter(species != "") 

richness_summary_pt1 <-  ibi_insect_richness %>% filter(order =="Plecoptera" |order =="Trichoptera") %>% 
  group_by(Sample, Reach, Biological.replicate) %>% 
  mutate(ibi_richness = length(unique(species))) %>%  # unique spp per IBI order
  ungroup() %>% 
  select(mmyy, Site, Reach, Biological.replicate, order, ibi_richness) %>% 
  rename("taxon" ="order") %>% 
  mutate("group" = "Trichoptera + Plecoptera")

richness_summary_pt2 <-  ibi_insect_richness %>% filter(family == "Chironomidae" |order =="Ephemeroptera") %>% 
  group_by(Sample, Reach, Biological.replicate) %>% 
  mutate(ibi_richness = length(unique(species))) %>%  # unique spp per IBI order
  ungroup() %>% 
  select(mmyy, Site, Reach, Biological.replicate, family, ibi_richness) %>% 
  rename("taxon" ="family") %>% 
  mutate("group" = "Chironomidae + Baetidae")

richness_summary <- rbind(richness_summary_pt1, richness_summary_pt2)
  
# plot
richplot<- ggplot(richness_summary, aes(x=Site, y=ibi_richness)) + 
  geom_violin(aes(fill=Site)) + 
  facet_wrap(~group)+ 
  scale_fill_brewer(palette = "Blues", name = "group")  +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) + 
  labs(y = "Species Richness", x = "") + 
  theme(legend.position="none") + 
  scale_x_discrete( labels= c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.3,
               colour = "gray6")+ 
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7))

ggsave(richplot, file = here("Figures", "insect_rich_sp_IBIgroups.png"), width = 6, height = 4)









