###=== Comparing Creek Communities among Bellingham Creeks======
# Looking at differences in insect communities across creek, month, and reach (with CAP and PERMANOVA, using presence/absence)
# Written by Helen Casendino (hcas1024@uw.edu) & Ezza 
# Created: 28 Jan 2022   Modified: 13 Apr 2022

# Dependencies 
library(tidyverse)
library(here)
library(vegan)
library(ggpubr)
library(ggrepel)
library(cowplot)

asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))

###====Data wrangling======

# there were technical replicates done on 05215SqmUp and 06214PadDn, but the different # of reads don't matter since I'm looking at p/a. So just averaging over them
covariates_hashes <- asv_reads_annotated %>% filter(class == "Insecta") %>%
  select(mmyy, Site, Reach, Biological.replicate, nReads, species) %>% 
  group_by(mmyy, Site, Reach, Biological.replicate, species) %>% summarise(nReads = mean(nReads)) %>% 
  filter(species != "") %>% 
  pivot_wider(names_from = species, values_from = nReads, values_fill = 0) %>% # each species has column with nReads info
  ungroup()

covariates <- covariates_hashes %>% select(c(mmyy, Site , Reach ,Biological.replicate)) # month, site, reach, bottle info
covariates$Site <- factor(covariates$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
covariates$mmyy <- factor(covariates$mmyy, labels = c("March", "April", "May", "June", "July", "August"))

hashes <- covariates_hashes %>% select(!c(mmyy, Site , Reach ,Biological.replicate)) # hash cols with binary read info (now they're species cols, but I don't want to change variable names rn)
for(i in 1:ncol(hashes)){
  hashes[,i]<-ifelse(hashes[,i] != 0,1,0)
}

###==== Fig 1: Site differences, sp vectors (capscale)======

cap1 <- capscale(hashes~ Site,data=covariates, distance= "jaccard")
sppscores(cap1) <- hashes

# taxa assigments
ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
rownames(CAPorder) <- NULL

#map on the classification of the ASVs to determine
CAP.ASVid <- CAPorder %>% distinct() %>% arrange(desc(CAP1))
Top50 <- CAP.ASVid %>% head(50) %>% rename("species" = "Hash")
write.csv(Top50, "CAPanalysis_Site_Top50DescASVs.csv")

# plot
cap2 <- cap1[["CCA"]][["wa"]] %>%
  as.data.frame() %>%
  bind_cols(covariates) %>% 
  ggplot(aes(x = CAP1,
             y = CAP2)) +
  geom_point(size = 1.5) +
  geom_point(aes(color = Site), size = 3) +
  geom_segment(aes(x = 0, y = 0,
     xend = CAP1,
     yend = CAP2), data = Top50[1:15,], color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
  geom_label_repel(size = 2, aes(x= CAP1  ,
                      y= CAP2 ,  label = species),
                   data = Top50[1:15,], fill = "orange", alpha = 0.75, max.overlaps = 18) + 
  scale_color_viridis_d(option = "viridis", begin = .1, end = 1) + 
  theme_bw() 

ggsave(file = here("Figures", "CAP_Site.png"), width = 8, height = 7)

###====Fig 2: Month differences, sp vectors(capscale)======

   cap1 <- capscale(hashes ~ mmyy,data=covariates, distance= "jaccard") 
   sppscores(cap1) <- hashes
   
   ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
   CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
   CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
   rownames(CAPorder) <- NULL
   CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
   Top50 <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1)) %>% head(50)
  
 cap2<- cap1[["CCA"]][["wa"]] %>% # plotting
     as.data.frame() %>%
     bind_cols(covariates) %>% 
     ggplot(aes(x = CAP1,
                y = CAP2)) +  geom_point(size = 1.5) +  geom_point(aes(color = mmyy), size = 3) +
     geom_segment(aes(x = 0, y = 0,
                      xend = CAP1,
                      yend = CAP2), data = Top50[1:40,], color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
     geom_label_repel(aes(x= CAP1  ,  y= CAP2 ,  label = species), data = Top50[1:40,], fill = "orange", alpha = 0.75, max.overlaps = 23, nudge_y = 0.025, nudge_x= 0.025 ) + 
     scale_color_viridis_d(option = "plasma", begin = .1, end = 0.9) + 
      labs(color = "Month") +  theme_bw()  +
        theme(legend.title = element_text(size = 15), legend.text = element_text(size =10)) 

ggsave(file = here("Figures", "CAP_Month.png"), width = 9, height = 7)

###====Fig 3: Port and Pad only, sp vectors (capscale)======

# step 1: re-wrangle data (PADDEN)
covariates_hashes_PAD <- asv_reads_annotated %>% filter(class == "Insecta") %>% filter(Site == "4Pad") %>%
  select(mmyy, Reach, Biological.replicate, Hash, nReads) %>% 
  group_by(mmyy, Reach, Biological.replicate, Hash) %>% summarise(nReads = mean(nReads)) %>% 
  pivot_wider(names_from = Hash, values_from = nReads, values_fill = 0) # each hash has column with nReads info

covariates_PAD  <- covariates_hashes_PAD  %>% select(c(mmyy , Reach ,Biological.replicate)) 
covariates_PAD$mmyy <- factor(covariates_PAD$mmyy, labels = c("March", "April", "May", "June", "July", "August"))

hashes_PAD <- covariates_hashes_PAD %>% select(!c(mmyy , Reach ,Biological.replicate)) # hash cols with binary read info 
for(i in 1:ncol(hashes_PAD)){
  hashes_PAD[,i]<-ifelse(hashes_PAD[,i] != 0,1,0)
}

# step 2: get CAP of month sp

cap1 <- capscale(hashes_PAD ~ mmyy,data=covariates_PAD, distance= "jaccard") 
sppscores(cap1) <- hashes_PAD

ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
rownames(CAPorder) <- NULL
CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
Top50 <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1)) %>% head(50)

cap2<- cap1[["CCA"]][["wa"]] %>% # plotting
  as.data.frame() %>%
  bind_cols(covariates_PAD) %>% 
  ggplot(aes(x = CAP1,
             y = CAP2)) +  geom_point(size = 1.5) +  geom_point(aes(color=mmyy, shape = Reach), size = 3) +
  geom_segment(aes(x = 0, y = 0,
                   xend = CAP1,
                   yend = CAP2), data = Top50, color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
  geom_label_repel(size = 3, aes(x= CAP1  ,  y= CAP2 ,  label = species), data = Top50, fill = "orange", alpha = 0.75, max.overlaps = 23, nudge_y = 0.025, nudge_x= 0.025 ) + 
  scale_color_viridis_d(option = "plasma", begin = 0, end = 0.9) + 
  labs(color = "Month") +  theme_bw()  +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size =10)) 

ggsave(file = here("Figures", "CAP_Month_PAD.png"), width = 11, height = 7)


# step 1: re-wrangle data (PORTAGE)
covariates_hashes_PRT <- asv_reads_annotated %>% filter(class == "Insecta") %>% filter(Site == "1Prt") %>%
  select(mmyy, Reach, Biological.replicate, Hash, nReads) %>% 
  group_by(mmyy, Reach, Biological.replicate, Hash) %>% summarise(nReads = mean(nReads)) %>% 
  pivot_wider(names_from = Hash, values_from = nReads, values_fill = 0) # each hash has column with nReads info

covariates_PRT  <- covariates_hashes_PRT  %>% select(c(mmyy , Reach ,Biological.replicate)) 
covariates_PRT$mmyy <- factor(covariates_PRT$mmyy, labels = c("March", "April", "May", "June", "July", "August"))

hashes_PRT <- covariates_hashes_PRT %>% select(!c(mmyy , Reach ,Biological.replicate)) # hash cols with binary read info 
for(i in 1:ncol(hashes_PRT)){
  hashes_PRT[,i]<-ifelse(hashes_PRT[,i] != 0,1,0)
}

# step 2: get CAP of month sp

cap1 <- capscale(hashes_PRT ~ mmyy,data=covariates_PRT, distance= "jaccard") 
sppscores(cap1) <- hashes_PRT

ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
rownames(CAPorder) <- NULL
CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
Top50 <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1)) %>% head(50)

cap2<- cap1[["CCA"]][["wa"]] %>% # plotting
  as.data.frame() %>%
  bind_cols(covariates_PRT) %>% 
  ggplot(aes(x = CAP1,
             y = CAP2)) +  geom_point(size = 1.5) +  geom_point(aes(color=mmyy, shape = Reach), size = 3) +
  geom_segment(aes(x = 0, y = 0,
                   xend = CAP1,
                   yend = CAP2), data = Top50, color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
  geom_label_repel(size = 3, aes(x= CAP1  ,  y= CAP2 ,  label = species), data = Top50, fill = "orange", alpha = 0.75, max.overlaps = 23, nudge_y = 0.025, nudge_x= 0.025 ) + 
  scale_color_viridis_d(option = "plasma", begin = 0, end = 0.9) + 
  labs(color = "Month") +  theme_bw()  +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size =10)) 

ggsave(file = here("Figures", "CAP_Month_PORT.png"), width = 11, height = 7)







###====For Ezza, up11 vs dn, up5 vs down=====

#up11
padrows <- which(covariates$Reach != "Up5" & covariates$Site == "Padden")
padcov <- springcovariates[padrows,]
padhash <- springhashes[padrows,]
cap1 <- capscale(padhash~ Reach,data=padcov, distance= "jaccard")

cap2 <- cap1[["CCA"]][["wa"]] %>%
  as.data.frame() %>%
  bind_cols(padcov) %>% 
  ggplot(aes(x = CAP1,
             y = CAP1)) +
  geom_point(size = 1.5) +
  geom_point(aes(color = Reach), size = 3) +
  theme_bw() +labs(title = "Padden: Dn vs. Up11")

#up5
padrows <- which(covariates$Reach != "Up11" & covariates$Site == "Padden")
padcov <- springcovariates[padrows,]
padhash <- springhashes[padrows,]
cap1 <- capscale(padhash~ Reach,data=padcov, distance= "jaccard")

cap3 <- cap1[["CCA"]][["wa"]] %>%
  as.data.frame() %>%
  bind_cols(padcov) %>% 
  ggplot(aes(x = CAP1,
             y = CAP1)) +
  geom_point(size = 1.5) +
  geom_point(aes(color = Reach), size = 3) +
  theme_bw() +labs(title = "Padden: Dn vs. Up5")


ggarrange(cap2, cap3, col = 2)
ggsave(file = here("Figures", "CAP_Up11_Up5_vs_Dn.png"), width = 6, height = 6)

###====Comparison across Creeks (permanova)======

adonis(hashes~ Site,data=covariates, method="jaccard")  # p = 0.001
adonis(hashes~ Site + mmyy + Site*mmyy,data=covariates, method="jaccard")  # p = 0.001








