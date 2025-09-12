install.packages("sf")
install.packages("raster")
install.packages("spData")
install.packages("leaflet")
install.packages("tmap")
install.packages("geobr")
remotes::install_github("Nowosad/spDataLarge")
devtools::install_github("humaniverse/geographr")

library(geographr)
library(sf)
library(raster)
library(dplyr)
library(spData)
library(spDataLarge)
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(ggplot2) # tidyverse data visualization package

WorldData <- map_data('world')
site= c('Arboretum','Peak District','Attenborough') 
lat = c(52.834139,53.365604,52.904268)  
long= c(-1.253199,-1.697046, -1.233926)
points =data.frame(cbind(site,lat,long))

ggplot() + 
  geom_polygon(data = WorldData, aes(x = long, 
                                    y = lat, 
                                    group = group), 
               fill = 'gray90', 
               color = 'black') + 
  coord_fixed(ratio = 1.3, 
              xlim = c(-10,3), 
              ylim = c(50, 59)) + 
  theme_void()+
  geom_point(data = points, 
             aes(x = as.numeric(long), 
                 y = as.numeric(lat),colour='red', alpha = .7))
