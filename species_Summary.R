###=== Visualizing species present (summary across levels) ======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 13 Feb 2022  Modified: 15 Apr 2022

# Dependencies 
library(tidyverse)
library(here)
library(formattable)

asv_reads_annotated <- read.csv(here("Input","COI_reads_taxonomy.csv"))
a <- asv_reads_annotated %>% filter(mmyy != "Kangaroo" & Site != "Run" & class == "Insecta") # insects only 


####==== Summary Statistics Table, by order ======

b <- a %>% group_by(order) %>%   mutate(n_unique_sp = length(unique(species))) %>% 
                                  mutate(n_unique_hash = length(unique(Hash)))
c <- b %>% ungroup() %>% 
        group_by(Hash) %>% mutate(sp.ID.logical = case_when(species == "" ~ "no",
                                            TRUE ~ "yes")) %>% ungroup() %>% group_by(order, Hash, sp.ID.logical) %>% 
        summarize(n_responses = n())%>% select(-n_responses) %>% ungroup() %>% group_by(order) %>% 
    summarize(n_IDd = length(which(sp.ID.logical == "yes")))
d <- left_join(b, c, by = "order")

summary_table <- d %>% mutate(IDd.hash.prop = n_IDd / n_unique_hash) %>% # so now I have divided, for each order, the number of Hashes ID'd to species by the # of unique hashes in that order
    ungroup() %>% 
      group_by(Sample, Reach, Biological.replicate) %>% # i also want to get average read proportion for each order 
      mutate(totalReads = sum(nReads)) %>% 
      ungroup() %>%  group_by(Sample, Reach, Biological.replicate, order) %>%
      mutate(orderReads = sum(nReads)) %>% 
      ungroup() %>% mutate(orderPropReads = orderReads/totalReads) %>% 
      group_by(order) %>% mutate(averaged.Order.Prop.Reads = mean(orderPropReads)) %>% 
      select(order, n_unique_sp, IDd.hash.prop, averaged.Order.Prop.Reads) %>% 
      distinct()

summary_table <-  summary_table %>% filter(order != "") %>% 
            rename("Order" = "order", "Unique.Species" = "n_unique_sp", "Prop.Hash.toSpeciesID" = "IDd.hash.prop", 
                   "Mean.Read.Proportion" = "averaged.Order.Prop.Reads" )
formattable::formattable(summary_table)

# Based on this summary table, we could include the following in the sp heat map:
# Trichoptera, Diptera, Plecoptera, Lepidoptera, Ephemeroptera
# also Coleoptera, Hemiptera (64 species, but some of these will only pop up up once or twice, and not in the CAP, so we can cut them )

# Trichoptera, Diptera, Ephemeroptera, then Plecoptera and Lepidoptera are taking up the most reads 

####==== "Heat Map" Plots for each site through months, based on common orders/sp ===== 

# Trichoptera, Diptera, Plecoptera, Lepidoptera, Ephemeroptera
# also Coleoptera, Hemiptera

# get data (a) into (Site mmyy Hash order species Normalized.reads)

# padden code (change)
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








