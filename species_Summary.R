###=== Visualizing species present (summary across levels) ======
# Written by Helen Casendino (hcas1024@uw.edu)
# Created: 13 Feb 2022  Modified: 16 Apr 2022

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

summary_table$Prop.Hash.toSpeciesID <- round(summary_table$Prop.Hash.toSpeciesID, 2)
summary_table$Mean.Read.Proportion <- round(summary_table$Mean.Read.Proportion, 2)

formattable::formattable(summary_table)

# Based on this summary table, we could include the following in the sp heat map:
# Trichoptera, Diptera, Plecoptera, Lepidoptera, Ephemeroptera
# also Coleoptera, Hemiptera (64 species, but some of these will only pop up up once or twice, and not in the CAP, so we can cut them )

# Trichoptera, Diptera, Ephemeroptera, then Plecoptera and Lepidoptera are taking up the most reads 

####==== "Heat Map" Plots for each site through months, based on common orders/sp ===== 

# Trichoptera, Diptera, Plecoptera, Lepidoptera, Ephemeroptera
# also Coleoptera, Hemiptera

# get data (a) into (Site mmyy order species Normalized.reads)
sp_plot_df <- a %>% select(Site, mmyy, order, species, nReads) %>%
        filter(order != "" & species != "") %>% 
        group_by(Site, mmyy, order, species) %>% summarize(nReads = mean(nReads)) %>% 
        mutate(p.a = case_when(nReads > 0 ~ 1,
                               TRUE ~ 0)) %>% select(-nReads) %>%
          unite(col = spord, c(order, species), sep = ".") %>% 
          pivot_wider(names_from = spord, values_from = p.a, values_fill = 0)

sp_plot_df %>% group_by(Site) %>% summarize(n = length(unique(mmyy))) # check that all months/sites are represented...Squalicum missing one month (august), so lets add a row
sp_plot_df[nrow(sp_plot_df) + 1, 1] <- "5Sqm"
sp_plot_df[nrow(sp_plot_df), 2] <- "0821"
sp_plot_df[nrow(sp_plot_df), 3:ncol(sp_plot_df)] <- 0

# pivot back to long 
sp_plot_df <- sp_plot_df %>% pivot_longer(!c(Site, mmyy) ,names_to = "spord", values_to = "p.a")  %>% 
  separate(col = "spord" , into = c("order", "species"), sep = "[.]")

sp_plot_df$Site <- factor(sp_plot_df$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
sp_plot_df$mmyy <- factor(sp_plot_df$mmyy, labels = c("March", "April", "May", "June", "July", "August"))
sp_plot_df <- sp_plot_df%>% filter(order != "Blattodea" & order!= "Odonata")

sp_plot_df %>% 
  mutate(species = fct_reorder(species, desc(order))) %>%
  ggplot(aes(x = mmyy, y = species)) +
  geom_point(shape=20, aes(size = ifelse(p.a==0, NA, p.a), color = order))  +
  facet_wrap( ~ Site, ncol = 5) + 
  theme_bw() +  
  labs(y = "", x = "", color = "Order") + 
  scale_color_brewer(palette = "Dark2")+
  guides(size = "none", color = guide_legend(override.aes = list(size = 9))) + 
  theme(axis.text.y = element_text(size = 4, face="italic"),
        axis.text.x = element_text(size = 8, angle = 45, hjust=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(file = here("Figures", "Species_Site_Month.png"), width = 7, height = 6.5)


####=== multi sp plot FROM WHEN I WANTED TO ONLY USE A SUBSET (old) ####

sp.avg.missed.mo <- sp_plot_df %>% group_by(Site, species, p.a) %>% summarize(n = n()) %>% 
  filter(p.a == "0") %>% # for each site, see for each sp how many months theyre absent from
  group_by(species) %>% 
  summarize(avg.missed.mo = mean(n)) # get avg for each sp 

hist(sp.avg.missed.mo$avg.missed.mo)
quantile(sp.avg.missed.mo$avg.missed.mo) # 5.3 is 25% quantile...lets cut it of at 5.3? 

# lets see what lepidopterans we should def include
lep_vec <- sp_plot_df %>% filter(order == "Lepidoptera") %>% group_by(species) %>% summarize(n = n()) %>% slice(3:14)
lep_vec <- c(lep_vec$species)

# back to making species list to put on plot
sp_to_incl <- sp.avg.missed.mo[sp.avg.missed.mo$avg.missed.mo <= 5.3,1]
sp_to_incl <- sp_to_incl$species
sp_to_incl <- c(sp_to_incl, "Baetis bicaudatus", "Chironomus whitseli", # add ephemeroptera and lepidoptera sp as well, and a few others that show up on vectors
                "Abagrotis baueri","Aquarius remigis",lep_vec)

sp_plot_df_reduced <- sp_plot_df %>% filter(species %in% sp_to_incl) # new plot df

# for plot display change mmyy and site labels
sp_plot_df_reduced$Site <- factor(sp_plot_df_reduced$Site, labels = c("Portage", "Barnes", "Chuckanut", "Padden", "Squalicum"  ))
sp_plot_df_reduced$mmyy <- factor(sp_plot_df_reduced$mmyy, labels = c("March", "April", "May", "June", "July", "August"))

# plot 
sp_plot_df_reduced %>% 
  mutate(species = fct_reorder(species, desc(order))) %>%
  ggplot(aes(x = mmyy, y = species)) +
  geom_point(shape=20, aes(size = ifelse(p.a==0, NA, p.a), color = order))  +
  facet_wrap( ~ Site) + 
  theme_bw() +  
  labs(y = "", x = "", color = "Order") + 
  scale_color_brewer(palette = "Dark2")+
  guides(size = "none", color = guide_legend(override.aes = list(size = 9))) + 
  theme(axis.text.y = element_text(size = 4, face="italic"),
        axis.text.x = element_text(size = 8, angle = 45, hjust=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(file = here("Figures", "Species_Site_Month_REDUCED.png"), width = 5.7, height = 6.5)


