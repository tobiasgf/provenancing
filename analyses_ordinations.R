# Analysis SoilTracker: ordinations
# Manuscript: Predicting provenance of forensic soil samples: 
# Linking soil to ecological habitats by metabarcoding and supervised classification
# Author: Tobias Guldberg Frøslev
# Date: 23-04-2019

library(here)
library(vegan)
library(readr)
library(dplyr)

#### load OTU tables ####
fun_otu_tab <- readRDS(here::here("data","samples_focussed_otu_tab.RDS")) # taken from the Man vs Machine study: https://github.com/tobiasgf/man_vs_machine/blob/master/data/samples_focussed_otu_tab.RDS
pla_otu_tab <- read.table(here::here("data","plant_table_final.txt"), sep="\t", header = T) # see earlier steps in separate script
euk_raw_derep <- readRDS(here::here("data","euk_raw_derep.RDS"))
ins_raw_derep <- readRDS(here::here("data","ins_raw_derep.RDS"))

#### make ordinations ####
pla_ord <- metaMDS(t(pla_otu_tab[,1:130]), k=6, try=200, trymax = 1000)
fun_ord <- metaMDS(t(fun_otu_tab[,1:130]), k=6, try=200, trymax = 1000)
euk_ord <- metaMDS(t(euk_raw_derep[,1:130]), k=6, try=200, trymax = 1000)
ins_ord <- metaMDS(t(ins_raw_derep[,1:130]), k=6, try=200, trymax = 1000)

col_dat <- data.frame(cbind(pla_ord$points[,1:4],fun_ord$points[,1:4],euk_ord$points[,1:4],ins_ord$points[,1:4]))
names(col_dat) <- c("b_plant1","b_plant2","b_plant3","b_plant4","b_fung1","b_fung2","b_fung3","b_fung4","b_eukar1","b_eukar2","b_eukar3","b_eukar4","b_insect1","b_insect2","b_insect3","b_insect4")

write.table(col_dat, here::here("data","ordination_results.txt"), sep="\t",quote=FALSE, col.names = NA)

#### count of named plant OTUs ####
pla_otu_tab <- read.table(here::here("data","plant_table_final.txt"), sep="\t", header = T)
pla_otu_tab %>% filter(!grepl("unmatched",species)) %>% filter(!grepl("ambigu",species)) %>% nrow() # 318 unambiguous plant names
