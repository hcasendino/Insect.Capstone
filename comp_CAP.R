###=== Comparing Creek Communities among Bellingham Creeks======
# Looking at differences in insect communities across creek, month, and reach (with CAP and PERMANOVA, using presence/absence)
# Written by Helen Casendino (hcas1024@uw.edu) & Ezza 
# Created: 28 Jan 2022   Modified: 4 Mar 2022

# Dependencies 
library(tidyverse)
library(here)
library(vegan)
library(ggpubr)
library(ggrepel)
library(cowplot)

asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
asv_reads_annotated <- asv_reads_annotated %>% filter(mmyy != "521" & mmyy != "621" & mmyy != "721") # FOR NOW, we'll remove may june and july because messed up sequencing runs

###====Data wrangling======

covariates_hashes <- asv_reads_annotated %>% filter(class == "Insecta") %>%
        select(mmyy, Site, Reach, Biological.replicate, Hash, nReads) %>% 
        pivot_wider(names_from = Hash, values_from = nReads, values_fill = 0) # each hash has column with nReads info

covariates <- covariates_hashes %>% select(c(mmyy, Site , Reach ,Biological.replicate)) # month, site, reach, bottle info
covariates$Site <- factor(covariates$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
covariates$mmyy <- factor(covariates$mmyy, labels = c("March", "April", "August"))

hashes <- covariates_hashes %>% select(!c(mmyy, Site , Reach ,Biological.replicate)) # hash cols with binary read info 
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
CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
CAP.ASVid <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1))
Top50 <- CAP.ASVid %>% head(50)
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
     yend = CAP2), data = Top50[1:20,], color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
  geom_label_repel(aes(x= CAP1  ,
                      y= CAP2 ,  label = species),
                   data = Top50[1:20,], fill = "orange", alpha = 0.75, max.overlaps = 13) + 
  scale_color_viridis_d(option = "viridis", begin = .1, end = 1) + 
  theme_bw() 

ggsave(file = here("Figures", "CAP_Site.png"), width = 9, height = 6)

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
                      yend = CAP2), data = Top50, color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
     geom_label_repel(aes(x= CAP1  ,  y= CAP2 ,  label = species), data = Top50, fill = "olivedrab3", alpha = 0.75, max.overlaps = 23, nudge_y = 0.025, nudge_x= 0.025 ) + 
     scale_color_viridis_d(option = "plasma", begin = .4, end = 0.9) + 
      labs(color = "Month") +  theme_bw()  +
        theme(legend.title = element_text(size = 15), legend.text = element_text(size =10)) 

ggsave(file = here("Figures", "CAP_Month.png"), width = 10, height = 8)

###====Fig 3: Padden & Portage month and reach, sp vectors (capscale)======

# PADDEN
padrows <- which(covariates$Site == "Padden")
padcovariates <- covariates[padrows,]
padhashes <- hashes[padrows,]
padcovariates[which(padcovariates$Reach == "Up11" | padcovariates$Reach == "Up5"), "Reach"] <- "Up" # group padden up sites

cap1 <- capscale(padhashes ~ Reach*mmyy,data=padcovariates, distance= "jaccard") # all cap analysis 
sppscores(cap1) <- padhashes
 
ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
rownames(CAPorder) <- NULL
CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
Top50 <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1)) %>% head(50)
  
 cap3 <- cap1[["CCA"]][["wa"]] %>% # plotting
    as.data.frame() %>%
    bind_cols(padcovariates) %>% 
    ggplot(aes(x = CAP1,
               y = CAP2)) +  geom_point(size = 1.5) +  geom_point(aes(shape=mmyy, color = Reach), size = 3) +
    geom_segment(aes(x = 0, y = 0,
                     xend = CAP1,
                     yend = CAP2), data = Top50, color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
    geom_label_repel(aes(x= CAP1  ,  y= CAP2 ,  label = species), data = Top50, fill = "mediumpurple1", alpha = 0.75, max.overlaps = 20, nudge_y = 0.1, nudge_x= 0.1 ) + 
    scale_color_viridis_d(option = "magma", begin = .2, end = .9) + 
    labs(title = "Padden", color = "Reach", shape = "Month") + theme_bw()  + 
   coord_cartesian(xlim = c(-0.25, 0.6), ylim = c(-0.5, 0.5)) + 
   theme(legend.position = "none")
 

# PORTAGE (ignore pad naming)
 padrows <- which(covariates$Site == "Portage")
 padcovariates <- covariates[padrows,]
 padhashes <- hashes[padrows,]
 padcovariates[which(padcovariates$Reach == "Up11" | padcovariates$Reach == "Up5"), "Reach"] <- "Up" # group padden up sites
 
 cap1 <- capscale(padhashes ~ Reach*mmyy,data=padcovariates, distance= "jaccard") # all cap analysis 
 sppscores(cap1) <- padhashes
 
 ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
 CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
 CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
 rownames(CAPorder) <- NULL
 CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
 Top50 <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1)) %>% head(50)
 
 cap4 <- cap1[["CCA"]][["wa"]] %>% # plotting
   as.data.frame() %>%
   bind_cols(padcovariates) %>% 
   ggplot(aes(x = CAP1,
              y = CAP2)) +  geom_point(size = 1.5) +  geom_point(aes(shape=mmyy, color = Reach), size = 3) +
   geom_segment(aes(x = 0, y = 0,
                    xend = CAP1,
                    yend = CAP2), data = Top50, color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
   geom_label_repel(aes(x= CAP1  ,  y= CAP2 ,  label = species), data = Top50, fill = "mediumpurple1", alpha = 0.75, max.overlaps = 24, nudge_y = 0.1, nudge_x= 0.1 ) + 
   scale_color_viridis_d(option = "magma", begin = .2, end = .9) + 
   labs(title = "Portage", shape = "Month", color = "Reach", y ="") + theme_bw()  + 
   coord_cartesian(xlim = c(-0.25, 0.6), ylim = c(-0.5, 0.5)) + 
   theme(legend.position = c(0.7, 0.15), legend.box = "horizontal",
         legend.background = element_rect(fill="white",size = 1, linetype="solid", 
                                          colour ="darkslateblue"))
ggarrange(cap3, cap4)
ggsave(file = here("Figures", "CAP_PortPad_Month_Reach.png"), width = 9, height = 6)

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








