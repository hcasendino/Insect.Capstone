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

#####==a=====

# close up map bounds
bounds<-c(left=-122.7 , bottom=48.15 , right=-122.1 , top=48.81)

get_stamenmap(bounds, zoom=10, maptype = "toner-lite") %>% ggmap()+
  geom_point(aes(x=-122.1286, y=48.18314), colour=1, size=5) + # portage
  geom_point(aes(x=-122.40405, y=48.79990), colour=2, size=5) + # squalicum
  geom_point(aes(x=-122.4751, y=48.71053), colour=3, size=5) + # padden
  geom_point(aes(x=-122.4094, y=48.68956), colour=4, size=5) + # chuckanut
  geom_point(aes(x=-122.369, y=48.66417), colour=5, size=5) + # barnes
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  #  labs(x="Latitude (Deg Dec)", y="Longitude (Deg Dec)") +
  theme(axis.text.y = element_text(size = 8)) + 
  theme(axis.text.x = element_text(size = 8))


# far away map 

bounds<-c(left=-125 , bottom=47 , right=-121 , top=49.5)
get_stamenmap(bounds, zoom=8, maptype = "watercolor") %>% ggmap()+
  geom_rect(xmin=-122.7 , ymin=48.15 , xmax=-122.1 , ymax=48.81, alpha=0.1, col="red")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 8))

