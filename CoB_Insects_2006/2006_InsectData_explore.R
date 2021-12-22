# =======WWU/City of Bellingham 2001-2003==========
# From Joan (co-author of City of Bellingham Urban Streams Final Report 2006)

# september 2002, 2001, 2003 (multiple sampling days/yr) at "Padden", "Squalicum", "Whatcom"...different sites at each watershed
# each species has column, also rep column includes both numbers and letter designations? 
insect_table <- read.csv(here("input/Bellingham_Urban_Streams_2006.csv"))

# read in species code key 
insect_key <- read.csv(here("old input files/Bellingham_sp_key.csv")) 
sp_key <- insect_key[-c(1:8, 156:223),]
colnames(sp_key) <- c("speciesCode", "order", "family", "genus")

# I want to get a sense of most abundant benthic inverts for each watershed, across year/day/site/replicate. Can do deeper diver later

invert_list <- list("Bellingham Watersheds - Benthic Inverts") # will hold ranked species lists
watersheds <- c("Padden", "Squalicum", "Whatcom")

for(i in 1:length(watersheds)){
  
  Grouped <- insect_table %>% 
    filter(watershed == watersheds[i]) %>% 
    select(-c(month,watershed, day, year, site, rep)) %>% 
    pivot_longer(cols = -page, names_to = "speciesCode", values_to = "number") %>% select(-page)
  
  species_rank <- Grouped %>% 
    group_by(speciesCode) %>% summarise(average = mean(number)) %>% 
    arrange(desc(average)) # in order of desceding count for each species 
  
  annotated <- left_join(species_rank, sp_key, by = "speciesCode")
  sp_rank_anno <- annotated %>% select(-speciesCode) %>% filter(average > 0)
  
  invert_list[[i+1]] <- sp_rank_anno
}


