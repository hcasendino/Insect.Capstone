###=== Comparing Creek Communities among Bellingham Creeks======
# this one reverts to using asvs for dissim matrices, not just ID'd sp
# Created: 28 Jan 2022   Modified: 13 Apr 2022 

# Dependencies 
library(tidyverse)
library(here)
library(vegan)
library(ggpubr)
library(ggrepel)
library(gganimate)
library(gifski)
library(av)

asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))

###====Data wrangling======

# there were technical replicates done on 05215SqmUp and 06214PadDn, but the different # of reads don't matter since I'm looking at p/a. So just averaging over them
covariates_hashes <- asv_reads_annotated %>% filter(class == "Insecta") %>%
  select(mmyy, Site, Reach, Biological.replicate, Hash, nReads) %>% 
  group_by(mmyy, Site, Reach, Biological.replicate, Hash) %>% summarise(nReads = mean(nReads)) %>% 
  pivot_wider(names_from = Hash, values_from = nReads, values_fill = 0) # each hash has column with nReads info


covariates <- covariates_hashes %>% select(c(mmyy, Site , Reach ,Biological.replicate)) # month, site, reach, bottle info
covariates$Site <- factor(covariates$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
covariates$mmyy <- factor(covariates$mmyy, labels = c("March", "April", "May", "June", "July", "August"))

hashes <- covariates_hashes %>% select(!c(mmyy, Site , Reach ,Biological.replicate)) # hash cols with binary read info 
for(i in 1:ncol(hashes)){
  hashes[,i]<-ifelse(hashes[,i] != 0,1,0)
}

###==== Fig 1 & 2: Site differences (across and by month), sp vectors (capscale)======

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
                   yend = CAP2), data = Top50[1:30,], color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
  geom_label_repel(size = 2, aes(x= CAP1  ,
                                 y= CAP2 ,  label = species),
                   data = Top50[1:30,], fill = "orange", alpha = 0.75, max.overlaps = 13) + 
  scale_color_viridis_d(option = "viridis", begin = .1, end = 1) + 
  theme_bw() 

ggsave(file = here("Figures", "CAP_Site.png"), width = 7, height = 5)


#### FIG 2 = SITE DIFFS ACCOUNTING FOR MONTH 
wrangle_by_month <- function(month1, month2){
  return(asv_reads_annotated %>% filter(class == "Insecta" & (mmyy == month1 | mmyy == month2)) %>%
           select(mmyy, Site, Reach, Biological.replicate, nReads, Hash) %>% 
           group_by(mmyy, Site, Reach, Biological.replicate, Hash) %>% summarise(nReads = mean(nReads)) %>% 
           pivot_wider(names_from = Hash, values_from = nReads, values_fill = 0) %>% # each species has column with nReads info
           ungroup())
}
get_CAP_plot <- function(covariates_hashes){ 
  
  covariates <- covariates_hashes %>% select(c(mmyy, Site , Reach ,Biological.replicate)) # month, site, reach, bottle info
  covariates$Site <- factor(covariates$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
  
  
  hashes <- covariates_hashes %>% select(!c(mmyy, Site , Reach ,Biological.replicate)) 
  for(i in 1:ncol(hashes)){
    hashes[,i]<-ifelse(hashes[,i] != 0,1,0)
  }
  
  cap1 <- capscale(hashes ~ Site, data=covariates, distance= "jaccard") 
  sppscores(cap1) <- hashes
  ASVvectors <- cap1[["CCA"]][["v"]] %>% as.data.frame 
  CAPorder <- ASVvectors %>% arrange(desc(CAP1)) 
  CAPorder <- cbind(Hash = rownames(CAPorder), CAPorder)
  rownames(CAPorder) <- NULL
  CAP.ASVid <- merge(CAPorder, asv_reads_annotated[,c("Hash","order", "family", "genus", "species")], by=c("Hash"), all.x=FALSE, all.y=TRUE)
  Top50 <- CAP.ASVid %>% distinct() %>% arrange(desc(CAP1)) %>% head(50)
  
  cap2 <- cap1[["CCA"]][["wa"]] %>%
    as.data.frame() %>%
    bind_cols(covariates) %>% 
    ggplot(aes(x = CAP1,
               y = CAP2)) +
    geom_point(size = 1.5) +
    geom_point(aes(color = Site), size = 2.5) +
    geom_segment(aes(x = 0, y = 0,
                     xend = CAP1,
                     yend = CAP2), data = Top50[1:20,], color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
    geom_label_repel(size = 2, aes(x= CAP1  ,
                                   y= CAP2 ,  label = species),
                     data = Top50[1:20,], fill = "orange", alpha = 0.75, max.overlaps = 18) + 
    scale_color_viridis_d(option = "viridis", begin = .1, end = 1) + 
    theme_bw() + 
    coord_cartesian(ylim = c(-0.5,0.5), xlim= c(-0.4,0.45))
  
  return(cap2)
}

covariates_hashes <- wrangle_by_month("0321","0421")
cap1 <- get_CAP_plot(covariates_hashes)

covariates_hashes <- wrangle_by_month("0521","0621")
cap2 <- get_CAP_plot(covariates_hashes)

covariates_hashes <- wrangle_by_month("0721","0821")
cap3 <- get_CAP_plot(covariates_hashes)

fullplot <- ggarrange(cap1 + ggtitle("March-April"), 
                      cap2 + theme(axis.title.y = element_blank())  +  ggtitle("May-June"), 
                      cap3 + theme(axis.title.y = element_blank()) + ggtitle("July-August"),
                      common.legend = TRUE, 
                      legend = "right", ncol =3)

ggsave(file = here("Figures", "CAP_Site_monthfacet.png"), width = 12, height = 5)

###====Fig 3: Month diffs only, sp vectors =====

covariates_hashes <- asv_reads_annotated %>% filter(class == "Insecta") %>%
  select(mmyy, Site, Reach, Biological.replicate, Hash, nReads) %>% 
  group_by(mmyy, Site, Reach, Biological.replicate, Hash) %>% summarise(nReads = mean(nReads)) %>% 
  pivot_wider(names_from = Hash, values_from = nReads, values_fill = 0) # each hash has column with nReads info

covariates <- covariates_hashes %>% select(c(mmyy, Site , Reach ,Biological.replicate)) # month, site, reach, bottle info
covariates$Site <- factor(covariates$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
covariates$mmyy <- factor(covariates$mmyy, labels = c("March", "April", "May", "June", "July", "August"))

hashes <- covariates_hashes %>% select(!c(mmyy, Site , Reach ,Biological.replicate)) # hash cols with binary read info 
for(i in 1:ncol(hashes)){
  hashes[,i]<-ifelse(hashes[,i] != 0,1,0)
}

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
                   yend = CAP2), data = Top50[1:15,], color = "black", arrow = arrow(length = unit(0.1,"cm"))) +
  geom_label_repel(size = 2, aes(x= CAP1  ,  y= CAP2 ,  label = species), data = Top50[1:15,], fill = "orange", alpha = 0.75, max.overlaps = 23) + 
  scale_color_viridis_d(option = "plasma", begin = .1, end = 0.9) + 
  labs(color = "Month") +  theme_bw()  +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size =10)) 

ggsave(file = here("Figures", "CAP_Month.png"), width = 7, height = 5)

##====Fig 4: Animated Site CAP =====

wrangle_by_month <- function(month1){
  return(asv_reads_annotated %>% filter(class == "Insecta" & (mmyy == month1)) %>%
           select(mmyy, Site, Reach, Biological.replicate, nReads, Hash) %>% 
           group_by(mmyy, Site, Reach, Biological.replicate, Hash) %>% summarise(nReads = mean(nReads)) %>% 
           pivot_wider(names_from = Hash, values_from = nReads, values_fill = 0) %>% # each species has column with nReads info
           ungroup())
}
get_CAP_df <- function( month , covariates_hashes){ 
  
  covariates <- covariates_hashes %>% select(c(mmyy, Site , Reach ,Biological.replicate)) 
  
  if(month != "0821"){
    covariates$Site <- factor(covariates$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"))
  } else {
    covariates$Site <- factor(covariates$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden"))
  }
  
  hashes <- covariates_hashes %>% select(!c(mmyy, Site , Reach ,Biological.replicate)) 
  for(i in 1:ncol(hashes)){
    hashes[,i]<-ifelse(hashes[,i] != 0,1,0)
  }
  
  cap1 <- capscale(hashes~ Site,data=covariates, distance= "jaccard")
  sppscores(cap1) <- hashes
  
  return(cap1[["CCA"]][["wa"]] %>%
           as.data.frame() %>%
           bind_cols(covariates))
}

# For each month, get CAP Site points
covariates_hashes <- wrangle_by_month("0321")
output1 <- get_CAP_df("0321", covariates_hashes)
output1 <- output1 %>% select(CAP1, CAP2, mmyy, Site) 

covariates_hashes <- wrangle_by_month("0421")
output2 <- get_CAP_df("0421", covariates_hashes)
output2 <- output2 %>% select(CAP1, CAP2, mmyy, Site) 

covariates_hashes <- wrangle_by_month("0521")
output3 <- get_CAP_df("0521", covariates_hashes)
output3 <- output3 %>% select(CAP1, CAP2, mmyy, Site) 

covariates_hashes <- wrangle_by_month("0621")
output4 <- get_CAP_df("0621", covariates_hashes)
output4 <- output4 %>% select(CAP1, CAP2, mmyy, Site) 

covariates_hashes <- wrangle_by_month("0721")
output5 <- get_CAP_df("0721", covariates_hashes)
output5 <- output5 %>% select(CAP1, CAP2, mmyy, Site) 

covariates_hashes <- wrangle_by_month("0821")
output6 <- get_CAP_df("0821", covariates_hashes)
output6 <- output6 %>% select(CAP1, CAP2, mmyy, Site) 

cap.coords <- rbind(output1, output2, output3, output4, output5, output6) # all together now

# Convert mmyy to date class (annoying oof)
nMoRows <- cap.coords %>% group_by(mmyy) %>% summarise(n =n()) %>% select(n)
cap.coords$mmyy <- c(rep("03", nMoRows$n[1]),
           rep("04", nMoRows$n[2]),
           rep("05", nMoRows$n[3]),
           rep("06", nMoRows$n[4]),
           rep("07", nMoRows$n[5]),
           rep("08", nMoRows$n[6]))

# Now plot animation 
p1 <- cap.coords %>%
  ggplot(aes(x = CAP1,
             y = CAP2)) +
  geom_point(aes(color = Site), size = 5) +
  scale_color_viridis_d(option = "viridis", begin = .1, end = 1) + 
  theme_bw()  +
  transition_states(as.factor(mmyy), transition_length = 7, state_length = 9) +
  labs(title = 'Month: {(next_state)}') + 
  ease_aes('linear') 

anim <- animate(p1, duration = 100, fps = 20, width = 400, height = 400, renderer = gifski_renderer(loop = F))
anim_save('output.gif')

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

