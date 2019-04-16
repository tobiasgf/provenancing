start.time <- Sys.time()

require(dada2)
require(readr)
require(dplyr)
require(digest)

pwd <- getwd()
setwd(pwd)
main_path <- pwd

variables <- read.table(file.path(main_path, "variables.txt"), stringsAsFactors = F)
trunc_reads <- as.logical(variables[1,2])
cut_r1 <- as.numeric(variables[2,2])
cut_r2 <- as.numeric(variables[3,2])

appendix <- as.character(paste0(cut_r1,"_",cut_r2))

#filtering of the Sense-reads:
path <- file.path(main_path, "DADA2_SS") 
fns <- list.files(path)
fastqs <- fns[grepl("fastq", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_R1.", fastqs)]
fnRs <- fastqs[grepl("_R2.", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_R1"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

if(trunc_reads){
SSout <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minLen=10, maxN=0, maxEE=c(2,2),truncQ=2,truncLen=c(cut_r1,cut_r2),compress=TRUE, multithread=TRUE,matchID=TRUE)
} else {
 SSout <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minLen=10, maxN=0, maxEE=c(2,2),truncQ=2,compress=TRUE, multithread=TRUE, matchIDs = TRUE)
}
head(SSout)

#filtering of the antiSense-reads:
path <- file.path(main_path, "DADA2_AS") 
fns <- list.files(path)
fastqs <- fns[grepl("fastq", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_R2.", fastqs)] # Reverse direction compared to above for the "sense reads". In practice the reads are here complement reversed to be in the same orientation as the "sense" reads.
fnRs <- fastqs[grepl("_R1.", fastqs)] # See above
sample.names <- sapply(strsplit(fnFs, "_R2"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

if(trunc_reads){
 ASout <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minLen=10, maxN=0, maxEE=c(2,2),truncQ=2,truncLen=c(cut_r2,cut_r1),compress=TRUE, multithread=TRUE,matchID=TRUE)
} else {
 ASout <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minLen=10, maxN=0, maxEE=c(2,2),truncQ=2,compress=TRUE, multithread=TRUE,matchIDs=TRUE)
}

head(ASout)

#Processing the set of files containing the forward primer in the R1 reads (the sense reads):
filt_path <- file.path(main_path, "DADA2_SS/filtered") 
fns <- list.files(filt_path)
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
SSsample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- SSsample.names
names(derepRs) <- SSsample.names
SSdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
SSdadaRs <- dada(derepRs, err=errR, multithread = TRUE)
SSmergers <- mergePairs(SSdadaFs, derepFs, SSdadaRs, derepRs, verbose=TRUE,minOverlap = 5)
seqtab_SS <- makeSequenceTable(SSmergers[names(SSmergers)])
seqtab.nochim_SS <- removeBimeraDenovo(seqtab_SS, verbose=TRUE)

stSS <- file.path(main_path,"seqtab_SS_RDS")
stnsSS <- file.path(main_path,"seqtab.nochim_SS_RDS")
saveRDS(seqtab_SS,stSS)
saveRDS(seqtab.nochim_SS,stnsSS)

#Then DADA2 processing of "the antisense" reads:
filt_path <- file.path(main_path, "DADA2_AS/filtered") 
fns <- list.files(filt_path)
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
ASsample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- ASsample.names
names(derepRs) <- ASsample.names
ASdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
ASdadaRs <- dada(derepRs, err=errR, multithread = TRUE)
ASmergers <- mergePairs(ASdadaFs, derepFs, ASdadaRs, derepRs, verbose=TRUE,minOverlap = 5)
seqtab_AS <- makeSequenceTable(ASmergers[names(ASmergers)])
seqtab.nochim_AS <- removeBimeraDenovo(seqtab_AS, verbose=TRUE)
stAS <- file.path(main_path,"seqtab_AS_RDS")
stnsAS <- file.path(main_path,"seqtab.nochim_AS_RDS")
saveRDS(seqtab_AS,stAS)
saveRDS(seqtab.nochim_AS,stnsAS)

#Define a function for combining two or more tables, collapsing samples with similar names:  
sumSequenceTables <- function(table1, table2, ..., orderBy = "abundance") {
 # Combine passed tables into a list
 tables <- list(table1, table2)
 tables <- c(tables, list(...))
 # Validate tables
 if(!(all(sapply(tables, dada2:::is.sequence.table)))) {
  stop("At least two valid sequence tables, and no invalid objects, are expected.")
 }
 sample.names <- rownames(tables[[1]])
 for(i in seq(2, length(tables))) {
  sample.names <- c(sample.names, rownames(tables[[i]]))
 }
 seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
 sams <- unique(sample.names)
 # Make merged table
 rval <- matrix(0L, nrow=length(sams), ncol=length(seqs))
 rownames(rval) <- sams
 colnames(rval) <- seqs
 for(tab in tables) {
  rval[rownames(tab), colnames(tab)] <- rval[rownames(tab), colnames(tab)] + tab
 }
 # Order columns
 if(!is.null(orderBy)) {
  if(orderBy == "abundance") {
   rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
  } else if(orderBy == "nsamples") {
   rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
  }
 }
 rval
}

stAS <- file.path(main_path,"seqtab_AS_RDS")
stnsAS <- file.path(main_path,"seqtab.nochim_AS_RDS")
stSS <- file.path(main_path,"seqtab_SS_RDS")
stnsSS <- file.path(main_path,"seqtab.nochim_SS_RDS")
seqtab.nochim_AS <- readRDS(stnsAS)
seqtab.nochim_SS <- readRDS(stnsSS)
seqtab_AS <- readRDS(stAS)
seqtab_SS <- readRDS(stSS)
total_sumtable <- sumSequenceTables(seqtab_SS,seqtab_AS)
total_nochim_sumtable <- sumSequenceTables(seqtab.nochim_SS,seqtab.nochim_AS)
stBoth <- file.path(main_path,"seqtab_Both_RDS")
stnsBoth <- file.path(main_path,"seqtab.nochim_Both_RDS")
saveRDS(total_sumtable,stBoth)
saveRDS(total_nochim_sumtable,stnsBoth)

extract_and_simplify <- function(dada2_table, out_name = "DADA2_nochim"){
 out_name <- file.path(getwd(),out_name)
 trans_total_nochim_sumtable <- as.data.frame(t(dada2_table))
 no_seqs <- dim(trans_total_nochim_sumtable)[1]
 no_samples <- dim(trans_total_nochim_sumtable)[2]
 sequences <- row.names(trans_total_nochim_sumtable)
 samples <- names(trans_total_nochim_sumtable)
 sha1_names <- sapply(sequences,digest::sha1)
 
 simple_tab_name <- paste0(out_name,".table")
 simple_table <- trans_total_nochim_sumtable
 rownames(simple_table) <- sha1_names
 {write.table(trans_total_nochim_sumtable,simple_tab_name,sep="\t",col.names = NA, quote=FALSE)}
 
 #Extract OTUs (sequences)
 sinkname <- paste0(out_name,".otus")
 {
  sink(sinkname)
  for (seqX in seq.int(length(sequences))) {
   header <- paste0(">",sha1_names[seqX],"\n")
   cat(header)
   seqq <- paste0(sequences[seqX],"\n")
   cat(seqq)
  }
  sink()
 }
 
 out_path <- paste0(out_name,"_extracted_samples")
 if(!file_test("-d", out_path)) dir.create(out_path)
 for (sampleX in seq(1:no_samples)){
  sinkname <- file.path(out_path, paste0(samples[sampleX],".fas"))
  {
   sink(sinkname)
   for (seqX in seq(1:no_seqs)) {
    if (dada2_table[sampleX,seqX] > 0) {
     header <- paste0(">",sha1_names[seqX],";size=",dada2_table[sampleX,seqX],";","\n")
     cat(header)
     seqq <- paste0(sequences[seqX],"\n")
     cat(seqq)
    }
   }
   sink()
  }
 }
}


extract_and_simplify(total_sumtable,out_name = "DADA2_raw")
extract_and_simplify(total_nochim_sumtable,out_name = "DADA2_nochim")


#Get statistics
getN <- function(x) sum(getUniques(x))
SSout <- as.data.frame(SSout)
SSout$sample <-  gsub("_R.\\.fastq","",row.names(SSout))
#SSout <- as.data.frame(SSout)
r1 <- data.frame(sample=names(SSdadaFs), denoisedF=sapply(SSdadaFs, getN),stringsAsFactors=FALSE)
r2 <- data.frame(sample=names(SSdadaRs), denoisedR=sapply(SSdadaRs, getN),stringsAsFactors=FALSE)
r3 <- data.frame(sample=names(SSmergers), merged=sapply(SSmergers, getN),stringsAsFactors=FALSE)
r4 <- data.frame(sample=row.names(seqtab.nochim_SS), nonchim=rowSums(seqtab.nochim_SS),stringsAsFactors=FALSE)

SStrack <- left_join(SSout, r1, by='sample') %>% left_join(., r2, by='sample') %>% left_join(., r3, by='sample') %>% left_join(., r4, by='sample')
row.names(SStrack) <- SStrack$sample
SStrack <- SStrack[,-3]

#saveRDS(SStrack, "read_track_SS_RDS")
write.table(SStrack,paste0("read_track_SS_",appendix,".txt"),col.names = NA, quote=FALSE,sep="\t")

ASout <- as.data.frame(ASout)
ASout$sample <-  gsub("_R.\\.fastq","",row.names(ASout))
#ASout <- as.data.frame(ASout)
r1 <- data.frame(sample=names(ASdadaFs), denoisedF=sapply(ASdadaFs, getN),stringsAsFactors=FALSE)
r2 <- data.frame(sample=names(ASdadaRs), denoisedR=sapply(ASdadaRs, getN),stringsAsFactors=FALSE)
r3 <- data.frame(sample=names(ASmergers), merged=sapply(ASmergers, getN),stringsAsFactors=FALSE)
r4 <- data.frame(sample=row.names(seqtab.nochim_AS), nonchim=rowSums(seqtab.nochim_AS),stringsAsFactors=FALSE)

AStrack <- left_join(ASout, r1, by='sample') %>% left_join(., r2, by='sample') %>% left_join(., r3, by='sample') %>% left_join(., r4, by='sample')
row.names(AStrack) <- AStrack$sample
AStrack <- AStrack[,-3]

#saveRDS(AStrack, "read_track_AS_RDS")
write.table(AStrack,paste0("read_track_AS_",appendix,".txt"),col.names = NA, quote=FALSE,sep="\t")

total_stat <- data.frame(merged_reads=rowSums(total_sumtable),merged_nochim=rowSums(total_nochim_sumtable))
#saveRDS(total_stat, "total_stat_RDS")
write.table(total_stat,paste0("total_stat_",appendix,".txt"),col.names = NA, quote=FALSE,sep="\t")

otu_stat <- data.frame(ss=ncol(seqtab_SS), ss_nochim=ncol(seqtab.nochim_SS), as=ncol(seqtab_AS), as_nochim=ncol(seqtab.nochim_AS), total=ncol(total_sumtable), total_nochim=ncol(total_nochim_sumtable))
write.table(otu_stat,paste0("otu_stat_",appendix,".txt"),col.names = NA, quote=FALSE,sep="\t")


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
