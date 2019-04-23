# Analysis SoilTracker: preparing mt16S insect OTU data
# Manuscript: Predicting provenance of forensic soil samples: 
# Linking soil to ecological habitats by metabarcoding and supervised classification
# Author: Tobias Guldberg Frøslev
# Date: 23-04-2019

library(here)
source(here::here("R","dereplicate.r"))

#### get data and make initial processing ####
#intermediate data from these steps not included in this GitHub repository

#download sequence data from (data from two MiSeq runs): https://sid.erda.dk/public/archives/65a54fa35a5b4621a834fc58a2ce2619/published-archive.html

#Combine reads from the two MiSeq runs, library wise.

#demultiplex the data using the tag-info files available in the seq_processing sub-directory "mt16S_insect" and the script: https://github.com/tobiasgf/man_vs_machine/blob/master/seq_processing/DADA2_demultiplex.sh

#analyse using DADA2 implemented in the script (seq_processing) : dada2_v3.1.R
#  R<dada2_v3.1.R --no-save &>log.txt #(needs info from file variables.txt in dir seq_processing/16S_insect)

#Rename and proceed with the dada2 table with chimeras removed:
#  mv DADA2_nochim.table ins_DADA2_nochim.table

#### dereplicate (collapse pcr triplicates for each sample) ####
ins_raw_tab <- read.table(here::here("data","ins_DADA2_nochim.table"), sep="\t", header = T,row.names = 1) # see earlier steps in separate script
ins_samples <- read.table(here::here("in_data","Sample_Site_Rep_Info_ins.txt"), sep="\t", header = T, stringsAsFactors = F)
ins_raw_derep <- dereplicate(ins_raw_tab, ins_samples, exclude_sites = c("BW","Bla","neg"))
saveRDS(ins_raw_derep, here::here("data","ins_raw_derep.RDS"))