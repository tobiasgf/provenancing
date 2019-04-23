# Analysis SoilTracker: preparing plant OTU data
# Manuscript: Predicting provenance of forensic soil samples: 
# Linking soil to ecological habitats by metabarcoding and supervised classification
# Author: Tobias Guldberg Frøslev
# Date: 23-04-2019
library(here)
library(lulu)
library(dplyr)
library(readr)
library(stringr)
source(here::here("R","dereplicate.r"))
source(here::here("R","Universal_taxonomic_assignment.r"))
options(ENTREZ_KEY = "ENTER_YOUR_KEY_HERE") # key needed for taxon retrievement from NCBI

#### get data and make initial processing ####
#intermediate data from these steps not included in this GitHub repository  
#download sequence data from: https://datadryad.org/resource/doi:10.5061/dryad.n9077  

#demultiplex the data using the tag-info files availablbe here: https://github.com/tobiasgf/lulu using a slightly improved script "DADA2_demultiplex.sh" from: https://github.com/tobiasgf/man_vs_machine/blob/master/seq_processing/DADA2_demultiplex.sh

#analyse using DADA2 implemented in custom script "dada2_sickle.r" from https://github.com/tobiasgf/man_vs_machine/blob/master/seq_processing/dada2_sickle.r  
#  R<dada2_sickle_v2.r --no-save &>log.t  

#Proceed with the samplewise extracted sequences (in directory DADA2_extracted_samples_nochim)  
#Extract the ITS2 regions with custom script:
#  bash itsx_plants_parallel.sh  

#Annotate the sequences with sample and otu ids   
#  mkdir out  
#  for f in *.fas; do awk '/>/{sub(">","&sample="FILENAME";otu=");sub(/\.fas/,x)}1' $f > out/$f ; done  

#Concatenate all sequences (in the out directory)   
#  cat *fas >> all_otus_concatenated.fas  

#Cluster at 97% level   
#  vsearch --cluster_size all_otus_concatenated.fas --id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --uc all.clustered.uc  --centroids all.otus.fasta --otutabout all_plant.otutab.txt  
#(all_plant.otutab.txt is included in this project)

#Get blasthits using custom script: split_and_blast.sh   
#  bash split_and_blast.sh  all.otus.fasta  
#(the resulting file "all.otus.fasta.blasthits" is included in this project)

#lulu matchlist   
#format sequences to only contain sequence ids   
#  sed 's/>.*otu=/>/' all.otus.fasta | sed 's/;size.*$//' > otu_clean_headers.fasta  

#make blast database  
#  makeblastdb -in otu_clean_headers.fasta -dbtype nucl  
#  blastn -db otu_clean_headers.fasta -num_threads 50 -outfmt '6 qseqid sseqid pident' -out plant_otu.matchlist -qcov_hsp_perc 80 -perc_identity 84 -query otu_clean_headers.fasta  
#(the resulting matchlist "plant_otu.matchlist" is included in this project)

#### lulu processing ####

#read otu tab and matchlist and perform lulu on the data
otutab <- read.csv(here::here("data","all_plant.otutab.txt"),sep='\t',header=T,row.names = 1, as.is=TRUE)
matchlist <-  read.csv(here::here("data","plant_otu.matchlist"),sep='\t',header=F,as.is=TRUE)
lulified_tab <- lulu(otutab,matchlist, minimum_match = 84, minimum_relative_cooccurence = 1)
saveRDS(lulified_tab,here::here("data","lulified_plant_tab.RDS"))

#### dereplicate ####
#collapse the three replicates for each sample
samplelist <- read.csv(here::here("data","Plants_ITS2_sample_sites.csv"), sep=";", header=TRUE, as.is = TRUE)
derep_tab <- dereplicate(otutab, samplelist)
saveRDS(derep_tab,here::here("data","derep_plant_tab.RDS"))

#### taxonomic assignment ####
IDtable <- read.csv(here::here("data","all.otus.fasta.blasthits"),sep='\t',header=F,as.is=TRUE)
names(IDtable) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","qcovs","sgi","sseq","ssciname","staxid")

IDtable <- IDtable %>% filter(qcovs>95, length/qlen > 0.95) # filter out matches that are not in full length.

my_classified_results <- assign_taxonomy(IDtable,upper_margin=0.5,lower_margin=2, remove = c("uncultured","environmental","N/A"), useDB=TRUE, appendDB=TRUE, db_path=here::here("data"), tax_db="storedTaxidsRDS", id_cut=c(98, 90, 85, 80, 75, 70, 50), ambiguity_cutoff=c(50,50,50,50,50,50,50)) # classify, using score 50 as cut off

write.table(my_classified_results$adjusted_classified_table, here::here("data","my_classified_otus_adjusted.txt"), sep = "\t", quote = F, row.names = F)

saveRDS(my_classified_results,here::here("data","adjusted_plant_tax.RDS"))

#### focus tables on plant otus and otus not flagged as errors ####
lulified_tab <- readRDS(here::here("data","lulified_plant_tab.RDS"))
derep_tab <- readRDS(here::here("data","derep_plant_tab.RDS"))
taxonomy <- readRDS(here::here("data","adjusted_plant_tax.RDS"))
tax_tab <- taxonomy$adjusted_classified_table
tax_tab$OTU_ID <- str_split_fixed(tax_tab$qseqid,";",3)[,2]
tax_tab$OTU_ID <- gsub("otu=","",tax_tab$OTU_ID)

pa_tab <- derep_tab
pa_tab[pa_tab>1] <- 1 # reduce to presence/absence table

tax_tab2 <- tax_tab %>% group_by(species) %>% mutate(redundancy = n()) # add count of each species name

lulified_to_keep <- tax_tab2$OTU_ID[tax_tab2$OTU_ID %in% lulified_tab$discarded_otus & tax_tab2$redundancy == 1 & tax_tab2$pident == 100] # which lulu-flagged otus should be kept (high match, non-redundant)
discard <- setdiff(lulified_tab$discarded_otus,lulified_to_keep) # which flagged errors should be removed

tax_tab3 <- tax_tab2[!tax_tab2$OTU_ID %in% discard,c(6,9,10,11,12,13,14,15,23,24,25)] # focus taxonomy on otus to keep

plant_tax <- tax_tab3 %>% filter(phylum == "Streptophyta", !class %in% c("Chlorophyta","Sphagnopsida","Jungermanniopsida","Bryopsida","Polytrichopsida","Klebsormidiophyceae"))

otu_tab <- derep_tab[plant_tax$OTU_ID,] # focus otu_tab on otus to keep
otu_tab$OTU_ID <- row.names(otu_tab)

otu_tab_pa <- pa_tab[plant_tax$OTU_ID,] # focus otu_tab on otus to keep
otu_tab_pa$OTU_ID <- row.names(otu_tab_pa)

edna_tab <- left_join(otu_tab, plant_tax, by = "OTU_ID") # add taxonomy
samples_otu_tab <- edna_tab[,-(131:141)]
controls_otu_tab <- edna_tab[,-(1:130)]

edna_tab_pa <- left_join(otu_tab_pa, plant_tax, by = "OTU_ID") # add taxonomy
samples_otu_tab_pa <- edna_tab_pa[,-(131:141)]
controls_otu_tab_pa <- edna_tab_pa[,-(1:130)]

write.table(edna_tab, here::here("data","plant_table_final.txt"), sep = "\t", quote = F, row.names = F)