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
     yend = CAP2), data = Top50[1:20,], color = "blue", arrow = arrow(length = unit(0.1,"cm"))) +
  geom_label_repel(aes(x= CAP1  ,
                      y= CAP2 ,  label = species),
                   data = Top50[1:20,], fill = "orange", alpha = 0.75, max.overlaps = 13) + 
  scale_color_viridis_d(option = "viridis", begin = .1, end = 1) + 
  theme_bw() 

ggsave(file = here("Figures", "CAP_Site.png"), width = 9, height = 6)

###====Fig 2: Month differences, faceted by site (capscale)======

plist <- list()
sitevec <- unique(covariates$Site)
for(i in 1:length(sitevec)){

   siterows <- which(covariates$Site == sitevec[i]) # everything corresponding to specified site
   cov <- covariates[siterows,]
   hash <- hashes[siterows,]
  
   cap1 <- capscale(hash ~ mmyy,data=cov, distance= "jaccard") # all cap analysis 
   sppscores(cap1) <- hash
   ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
   CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
   CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
   rownames(CAPorder) <- NULL
   CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
   Top50 <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1)) %>% head(50)
  
   plist[[i]] <- cap1[["CCA"]][["wa"]] %>% # plotting
     as.data.frame() %>%
     bind_cols(cov) %>% 
     ggplot(aes(x = CAP1,
                y = CAP1)) +  geom_point(size = 1.5) +  geom_point(aes(color = mmyy), size = 3) +
     geom_segment(aes(x = 0, y = 0,
                      xend = CAP1,
                      yend = CAP1), data = Top50[1:10,], color = "blue", arrow = arrow(length = unit(0.1,"cm"))) +
     geom_label_repel(aes(x= CAP1  ,  y= CAP1 ,  label = species), data = Top50[1:10,], fill = "orange", alpha = 0.75, max.overlaps = 8, nudge_y = 0.025, nudge_x= 0.025 ) + 
     scale_color_viridis_d(option = "plasma", begin = .4, end = 0.9) + 
      labs(color = "Month") + theme_bw() +  theme(legend.position="none")
   }

ggarrange(plist[[1]] + labs(title = "Portage", x = "") , 
          plist[[2]] + labs(title = "Barnes", x = "", y = "")  , 
          plist[[3]]+ labs(title = "Chuckanut", y = "") + theme(legend.position = "right", legend.title = element_text(size=15), legend.text = element_text(size=12))  ,
          plist[[4]]+ labs(title = "Squalicum") ,
          plist[[5]] + labs(title = "Padden", y = "") ,
          col=3)
ggsave(file = here("Figures", "CAP_Site_Month.png"), width = 12, height = 10)


###====Fig 3: Reach differences, faceted by site, where august is excluded or included (NOT WORKING)(capscale)======

# SPRING

springrows <- which(covariates$mmyy != "821")
springcovariates <- covariates[springrows,]
springhashes <- hashes[springrows,]
springcovariates <- springcovariates %>% mutate(Site = case_when(Site == "Barnes" | Site == "Chuckanut" | Site == "Squalicum" ~ "unrestored",
                                                                 Site == "Padden" ~ "Padden",
                                                                 Site == "Portage" ~ "Portage"))

springcovariates[which(springcovariates$Reach == "Up11" | springcovariates$Reach == "Up5"), "Reach"] <- "Up" # for richness, group padden up sites

plist <- list()
sitevec <- unique(springcovariates$Site)
for(i in 1:length(sitevec)){
  
  siterows <- which(springcovariates$Site == sitevec[i]) # everything corresponding to specified site
  cov <- springcovariates[siterows,]
  hash <- springhashes[siterows,]
  
  cap1 <- capscale(hash ~ Reach,data=cov, distance= "jaccard") # all cap analysis 
  sppscores(cap1) <- hash
  ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
  CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
  CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
  rownames(CAPorder) <- NULL
  CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
  Top50 <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1)) %>% head(50)
  
  plist[[i]] <- cap1[["CCA"]][["wa"]] %>% # plotting
    as.data.frame() %>%
    bind_cols(cov) %>% 
    ggplot(aes(x = CAP1,
               y = CAP1)) +  geom_point(size = 1.5) +  geom_point(aes(color = Reach), size = 3) +
    geom_segment(aes(x = 0, y = 0,
                     xend = CAP1,
                     yend = CAP1), data = Top50[1:10,], color = "blue", arrow = arrow(length = unit(0.1,"cm"))) +
    geom_label_repel(aes(x= CAP1  ,  y= CAP1 ,  label = species), data = Top50[1:10,], fill = "orange", alpha = 0.75, max.overlaps = 25, nudge_y = 0.1, nudge_x= 0.1 ) + 
    scale_color_viridis_d(option = "viridis", begin = 0, end = 1) + 
    labs(color = "Reach") + theme_bw() +  theme(legend.position="none")
}





# FALL (ignore labels of spring)

springrows <- which(covariates$mmyy == "821")
springcovariates <- covariates[springrows,]
springhashes <- hashes[springrows,]
springcovariates <- springcovariates %>% mutate(Site = case_when(Site == "Barnes" | Site == "Chuckanut" | Site == "Squalicum" ~ "unrestored",
                                                                 Site == "Padden" ~ "Padden",
                                                                 Site == "Portage" ~ "Portage"))

springcovariates[which(springcovariates$Reach == "Up11" | springcovariates$Reach == "Up5"), "Reach"] <- "Up" # for richness, group padden up sites

plist2 <- list()
sitevec <- unique(springcovariates$Site)
for(i in 1:length(sitevec)){
  
  siterows <- which(springcovariates$Site == sitevec[i]) # everything corresponding to specified site
  cov <- springcovariates[siterows,]
  hash <- springhashes[siterows,]
  
  cap1 <- capscale(hash ~ Reach,data=cov, distance= "jaccard") # all cap analysis 
  sppscores(cap1) <- hash
  ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
  CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
  CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
  rownames(CAPorder) <- NULL
  CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
  Top50 <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1)) %>% head(50)
  
  plist2[[i]] <- cap1[["CCA"]][["wa"]] %>% # plotting
    as.data.frame() %>%
    bind_cols(cov) %>% 
    ggplot(aes(x = CAP1,
               y = CAP1)) +  geom_point(size = 1.5) +  geom_point(aes(color = Reach), size = 3) +
    geom_segment(aes(x = 0, y = 0,
                     xend = CAP1,
                     yend = CAP1), data = Top50[1:10,], color = "blue", arrow = arrow(length = unit(0.1,"cm"))) +
    geom_label_repel(aes(x= CAP1  ,  y= CAP1 ,  label = species), data = Top50[1:10,], fill = "orange", alpha = 0.75, max.overlaps = 25, nudge_y = 0.1, nudge_x= 0.1 ) + 
    scale_color_viridis_d(option = "viridis", begin = 0, end = 1) + 
    labs(color = "Reach") + theme_bw() +  theme(legend.position="none")
}



# PLOT ALL 
ggarrange(plist[[1]] + labs(title = "Portage - Spring", x = "") , 
          plist[[2]] + labs(title = "Unrestored - Spring", x = "", y = "")  , 
          plist[[3]]+ labs(title = "Padden - Spring", y = "", x= "") + theme(legend.position = "right"),
         # plist2[[1]]+ labs(title = "Portage") ,
         # plist2[[2]] + labs(title = "Unrestored", y = "") ,
        #  plist2[[3]] + labs(title = "Padden", y = ""),
          row=1)
ggsave(file = here("Figures", "CAP_Spring_Site_Reach.png"), width = 10, height = 8)

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








