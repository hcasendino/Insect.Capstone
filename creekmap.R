#####==Dependencies=====

library(maps)
library(mapdata)
library(maptools) 
library(scales)
library(rgdal)
library(ggmap)
library(ggsn)
library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(ggplot2)
library(here)

#####==a=====

# close up map bounds
bounds<-c(left=-122.7 , bottom=48.15 , right=-122.1 , top=48.81)

plot1 <- get_stamenmap(bounds, zoom=10, maptype = "toner-lite") %>% ggmap()+
  geom_point(aes(x=-122.1286, y=48.18314), colour=1, size=5) + # portage
  geom_label(aes(label="Portage", colour= "black", fontface = "bold", x=-122.22, y=48.18314), size = 2.5) + 
  geom_point(aes(x=-122.40405, y=48.79990), colour=2, size=5) + # squalicum
  geom_label(aes(label="Squalicum", colour= "black", fontface = "bold", x=-122.3127, y=48.79990), size = 2.5) + 
  geom_point(aes(x=-122.4751, y=48.71053), colour=3, size=5) + # padden
  geom_label(aes(label="Padden", colour= "black", fontface = "bold", x=-122.3837, y=48.71053), size = 2.5) + 
  geom_point(aes(x=-122.4094, y=48.68956), colour=4, size=5) + # chuckanut
  geom_label(aes(label="Chuckanut", colour= "black", fontface = "bold", x=-122.318, y=48.68956), size = 2.5) + 
   geom_point(aes(x=-122.369, y=48.66417), colour=5, size=5) + # barnes
  geom_label(aes(label="Barnes", colour= "black", fontface = "bold", x=-122.2776, y=48.66417), size = 2.5) + 
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  labs(x="Latitude (Deg Dec)", y="Longitude (Deg Dec)") +
  theme(axis.text.y = element_text(size = 8)) + 
  theme(axis.text.x = element_text(size = 8))


# far away map 

bounds<-c(left=-125 , bottom=47 , right=-121 , top=49.5)
plot2 <- get_stamenmap(bounds, zoom=8, maptype = "toner-lite") %>% ggmap()+
  geom_rect(xmin=-122.7 , ymin=48.15 , xmax=-122.1 , ymax=48.81, alpha=0.1, col="red")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 8))

plot_grid(plot1, plot2, ncol = 2)
ggsave(file = here("Figures", "map.png"), width = 7, height = 7)
