# Predicting provenance of forensic soil samples - bioinformatic sequence analyses and ecological/provenancing modelling    
___

This repository (R project) contains all data (and links to data) and scripts necessary to run the initial bioinformatic sequence processing analyses and statistical analyses (procenancing anayses) from the study  **Predicting provenance of forensic soil samples: Linking soil to ecological habitats by metabarcoding and supervised classification**  (in review for publication in PLoS ONE).  
All steps/processes for this study can be carried out on the same computer/platform. But, in practise many of the bioinformatic analyses were carried out on a linux server setup with 64 processors (AMD Opteron(tm) 6380), except R-scripts, which were run on a MacBook Pro (2.9 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3).
All analyses were carried out in one directory and sub-directories of this.

## Sequencing data
The Illumina sequence data used in this study is deposited in the following repositories:   

#### Fungal data - originally published and processed in Frøslev et al. (2019)  
https://sid.erda.dk/public/archives/b0d2b22cb7804ff23d1612f4afdc29ae/published-archive.html

#### Plant data - originally published in Frøslev et al. (2017)  
https://datadryad.org/resource/doi:10.5061/dryad.n9077  

#### Eukaryote data - published here
https://sid.erda.dk/public/archives/3d227ac305559c1651c631d447983ff9/published-archive.html

#### Insect data - published here
https://sid.erda.dk/public/archives/65a54fa35a5b4621a834fc58a2ce2619/published-archive.html

## Bioinformatic tools
### CLI tools were used for this study  

 * VSEARCH v.2.9 (or later) (https://github.com/torognes/vsearch) 
 * Cutadapt v 1.17 (https://cutadapt.readthedocs.io/en/stable/)  
 * blastn v2.4.0+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 
 * sickle v.1.33 (https://github.com/najoshi/sickle)
 
Various R-packages were used for this study (see in the relevant R files).  
Variuous custom bash scripts used, are available in the seq_processing sub directory, and R-function available in the sub directory "R"  

## Description of the sub directories  

 * in_data : contains all initial data not produced as part of the analyses here  
 * data : contains all the data produced as part of the analyses here  
 * seq_processing : contains all the scripts and files necessary to perform the initial sequence processing  
 * R : contains a few functions used in the analyses  
 * plots : output directory for the plots/figures produced in the analyses
 
 
## References
Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Hansen, A. J., Læssøe, T., & Heilmann-Clausen, J. (2019). Man against machine: Do fungal fruitbodies and eDNA give similar biodiversity assessments across broad environmental gradients?. Biological Conservation, 233, 201-212.  
Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature communications, 8(1), 1188.  