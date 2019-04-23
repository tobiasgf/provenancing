# Analysis SoilTracker: map figure
# Manuscript: Predicting provenance of forensic soil samples: 
# Linking soil to ecological habitats by metabarcoding and supervised classification
# Author: Tobias Guldberg Frøslev
# Date: 23-04-2019

library(ggplot2)
library(rgdal)
library(raster)
library(readxl)
library(here)
library(dplyr)
library(tidyr)

#downlod map data from https://ec.europa.eu/eurostat
#https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/countries/download/ref-countries-2016-01m.shp.zip  

world <- readOGR(dsn = here::here("in_data","CNTR_RG_01M_2016_4326.shp"))

#cut to dk area
dk_basic <- crop(world, extent(8, 13, 54.5, 57.8))
dk_basic = fortify(dk_basic)
dk_basic$id[dk_basic$id != "49"] <- "0"  # making all other countries the same for mapping colour

#Combine data for the points
Biowidesites <- read_xlsx(here::here("in_data","DD130sites.xlsx")) # latlong_coordinated for the Biowide plots.
Cam_data <- read.table(here::here("in_data","RareSpecies_position.txt"))
col_data <- cbind(Biowidesites,Cam_data)
col_data2 <- col_data %>% select(-utm_x, -utm_y, -UTMX, -UTMY, -SiteNr, -Site_nr) %>% gather(key="parameter", value = "number", -c(ddlat,ddlong))

#Make theme for plotting
new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(0, 0, 0, 0), unit = "lines", valid.unit = 3L, class = "unit")

#make plot
plot1 <- ggplot() + 
 geom_polygon(data=dk_basic, aes(x=long,y=lat,group=group), fill="gray80") + 
 geom_point(data=col_data2,aes(ddlong,ddlat, size=number, fill=parameter),pch = 21, stroke = 0.1, alpha = 0.6, position = position_jitter(w = 0.02, h = 0.02)) +  
 scale_fill_manual(name="Distribution", breaks = c("rareSpsSSE","rareSpsE","rareSpsW"), labels=c("S/SE","E","W"), values=c('darkseagreen4','dodgerblue4','brown')) +
 #geom_path(data=dk_basic, aes(x=long,y=lat,group=group), colour='black', size = 0.5) + 
 new_theme_empty + 
 theme(panel.background=element_rect(fill = "white"), legend.position = "right") +
 coord_fixed(ratio = 1.79) + 
 guides(fill = guide_legend(order=1, override.aes = list(size=5)), size = guide_legend(order=2)) + 
 scale_size_continuous(name="Number of plants", range = c(1, 20), breaks =c(1,2,4,10)) + facet_wrap(.~parameter)

ggsave(here::here("plots","FigureS1n.pdf"), plot1, width = 12, height = 4, device = "pdf")