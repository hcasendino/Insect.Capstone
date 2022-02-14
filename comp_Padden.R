###=== Examining Padden ======
# Looking at differences in insect communities across month and reach for Padden Specifically (using updated eDNA index)
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 13 Feb 2022  Modified: 13 Feb 2022

# Dependencies 
library(tidyverse)
library(here)
library(vegan)
library(wesanderson)

asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
asv_reads_annotated <- asv_reads_annotated %>% filter(mmyy != "521" & mmyy != "621" & mmyy != "721") # FOR NOW, we'll remove may june and july because messed up sequencing runs


########===== (1) CAP on reach*month & PERMANOVA, plot of top 15 sp driving month diffs==== 

# a) index by max order proportion, not hash 

mean_hash_props <- asv_reads_annotated %>% 
  group_by(Sample,Biological.replicate) %>% 
  mutate (Tot = sum(nReads)) %>% 
  group_by(Sample,Biological.replicate, Hash) %>% 
  mutate(Row.prop = nReads / Tot) %>% # This creates the proportion on each biological replicate
  group_by(Sample, Hash) %>% 
  summarise (mean.prop = mean(Row.prop))  #"Averaging ratios between Biological replicates

order_index_insects_padden <- left_join(mean_hash_props, asv_reads_annotated, by = c("Hash", "Sample")) %>% 
  select(Sample,Hash, mean.prop, class, order) %>% distinct() %>% 
  filter(class == "Insecta" & order != "" & !is.na(order)) %>% 
  group_by (order) %>%
  mutate (Colmax = max (mean.prop),
          Normalized.reads = mean.prop / Colmax) %>% 
  dplyr::select(-Colmax, -mean.prop) %>%  # scales index vals to max in their order 
  separate(Sample, into = c("mmyy", "Site"), sep="_") %>% 
  filter(Site == "4Pad") %>% unite(col = Sample, c("mmyy", "Site") , sep="_" )

complete_asv_covariate_tab <- left_join(order_index_insects_padden, asv_reads_annotated, by = c("Sample", "Hash")) %>% 
  select(Site, mmyy, Reach, Hash, Normalized.reads) %>% distinct() %>% 
  pivot_wider(names_from = Hash, values_from = Normalized.reads) %>% 
  mutate_all(~replace(., is.na(.), 0))   

hashes_only <- complete_asv_covariate_tab %>% select(-c(Site, mmyy, Reach))
covariates_only <- complete_asv_covariate_tab %>% select(Site, mmyy, Reach)

covariates_only[which(covariates_only$Reach == "Up11" | covariates_only$Reach == "Up5"), "Reach"] <- "Up"
covariates_only$mmyy <- as.factor(covariates_only$mmyy )

# b) Month*reach CAP

cap1 <- capscale(hashes_only~ Reach + mmyy*Reach ,data=covariates_only, distance="bray")
sppscores(cap1) <- hashes_only

adonis(hashes_only~ Reach + mmyy*Reach,data=covariates_only, method="bray") 



cap2 <- cap1[["CCA"]][["wa"]] %>%
  as.data.frame() %>%
  bind_cols(covariates_only) %>%
  ggplot(aes(x = CAP1,
             y = CAP2)) +
  geom_point(size = 1.5) +
  geom_point(aes(color = mmyy, shape = Reach), size = 3) +
  scale_color_manual(values=wes_palette(n=3, name="FantasticFox1")) + 
  theme_bw() + 
  geom_text(x = -0.2, y = 0.35, label="Padden: Reach + mmyy*Reach", size = 4)

ggsave(file = here("Figures", "Padden_CAPplot_month_reach.png"), width = 5, height = 4)


cap1 <- capscale(hashes_only~ mmyy ,data=covariates_only, distance="bray")
sppscores(cap1) <- hashes_only

ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
rownames(CAPorder) <- NULL
CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
CAP.ASVid <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1))
Top50 <- CAP.ASVid %>% head(50)

Top15 <- Top50[1:15,"Hash"]
plotdata<-  order_index_insects_padden %>% filter(Hash %in% Top15) %>% 
  separate(Sample, into = c("mmyy", "Site"), sep="_") %>% 
  select(!Site) %>% 
  pivot_wider(names_from = mmyy, values_from = Normalized.reads, values_fill = 0) %>%
  pivot_longer(!c(Hash,class,order), names_to = "mmyy", values_to ="Normalized.reads" )

plotdata$Hash <- factor(plotdata$Hash, 
                        labels = c("Ephemeroptera", 
                                   "Barypeithes pellucidus",
                                   "Baetidae",
                                   "Lepidostoma unicolor", 
                                   "Aquarius remigis",
                                   "Lepidostoma",
                                  " Pteronarcys princeps",
                                  "Eucallipterus tiliae",
                                  "Campaea",
                                  "Neophylax",
                                  "Chironomidae",
                                 "Psychoglypha alascensis",
                                 "Ceratopsyche",
                                 "Rhyacophila vedra",
                                 "Parapsyche almota"))
                                  

plotdata %>% ggplot(aes(x = mmyy, y = Hash)) +
  geom_point(aes(size = Normalized.reads, color = Hash))  +
  theme_bw() +  labs(title = "Padden", y = "Taxon", x = "Month", size = "eDNA Index") + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  guides(colour = "none") 

ggsave(file = here("Figures", "Padden_15DiffTaxa_Month.png"), width = 5, height = 7)



########===== (2)  Fish ! (presence absence ) ==== 
sp <- c("Oncorhynchus kisutch","Cottus asper", "Perca flavescens")        

fishspdf <- asv_reads_annotated %>% filter(species %in% sp) %>% filter(Site == "4Pad") %>% 
  select(!c(Site, Sample,taxID,taxon,rank,score, kingdom  , phylum ,  class,order ,family, genus)) 
fishspdf[which(fishspdf$Reach == "Up11" | fishspdf$Reach == "Up5"), "Reach"] <- "Up"
fishspdf$mmyy <- as.factor(fishspdf$mmyy)

fishspdf <- fishspdf %>% 
  group_by(mmyy, species, Reach) %>% summarise(mean_read_count = mean(nReads)) %>% # not actually using this mean, just for P/A later on cause there are no 0s rn
  pivot_wider(names_from = "species", values_from =  mean_read_count, values_fill = 0) %>% 
  pivot_longer(!c(mmyy,Reach), names_to = "species", values_to = "p_a") %>% 
  mutate(p_a = case_when(p_a == 0 ~ "absence",
                         TRUE ~ "presence"))

fishspdf$mmyy <- factor(fishspdf$mmyy, labels = c("March", "April", "August"))

fishspdf %>% filter(p_a != "absence") %>% 
  ggplot(aes(x = Reach, y = species)) +
  geom_point(aes(color = species, shape = Reach), size =5) + 
  scale_color_brewer(palette = "Paired")+
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) + 
  facet_wrap( ~ mmyy) + 
  labs(title = "Padden", y = "", x = "Reach")  


ggsave(file = here("Figures", "Padden_Fishes_Reach_Month.png"), width = 7 , height = 3)


