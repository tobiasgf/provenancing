cut_r1 <- 310
cut_r2 <- 310
appendix <- as.character(paste0(cut_r1,"_",cut_r2))
require(dada2)
require(readr)
require(dplyr)

pwd <- getwd()
setwd(pwd)
main_path <- pwd

start.time <- Sys.time()

#Filtering the "Sense"-reads.
path <- file.path(main_path, "DADA2_SS")
fns <- list.files(path)
fastqs <- fns[grepl("fastq$", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_R1.", fastqs)]
fnRs <- fastqs[grepl("_R2.", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_R1"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
match_path <- file.path(path, "filtered")
if(!file_test("-d", match_path)) dir.create(match_path)
filtFs <- file.path(match_path, paste0(sample.names, "_F_filtered.fastq"))
filtRs <- file.path(match_path, paste0(sample.names, "_R_filtered.fastq"))
singles <- file.path(match_path, paste0(sample.names, "_singles.fastq"))

#commandX <- paste0("/usr/local/home/tobiasgf/bin/sickle-master/sickle pe -l 50 -q 28 -x -f -t sanger ",file.path(fnFs[i])," -r ",fnRs[i]," -o ",filtFs[i]," -p ",filtRs[i])

for(i in seq_along(fnFs)) {
 if (file.info(fnFs[i])$size == 0) {
  print(paste(fnFs[i], "empty.")) } else {
   print(paste("processing", fnFs[i]))
   commandX <- paste0("/usr/local/home/tobiasgf/bin/sickle-master/sickle se -l 100 -q 28 -x -t sanger -f ",fnFs[i]," -o ",filtFs[i])
   commandY <- paste0("/usr/local/home/tobiasgf/bin/sickle-master/sickle se -l 100 -q 28 -x -t sanger -f ",fnRs[i]," -o ",filtRs[i])
   print(commandX)
   system(commandX)
   print(commandY)
   system(commandY)
  }
}


#Filtering the "Anti-Sense"-reads.
path <- file.path(main_path, "DADA2_AS")
fns <- list.files(path)
fastqs <- fns[grepl("fastq$", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_R2.", fastqs)] # Reverse direction compared to above for the "sense reads". In practice the reads are here complement reversed to be in the same orientation as the "sense" reads.
fnRs <- fastqs[grepl("_R1.", fastqs)] # See above
sample.names <- sapply(strsplit(fnFs, "_R2"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
match_path <- file.path(path, "filtered")
if(!file_test("-d", match_path)) dir.create(match_path)
filtFs <- file.path(match_path, paste0(sample.names, "_F_filtered.fastq"))
filtRs <- file.path(match_path, paste0(sample.names, "_R_filtered.fastq"))
singles <- file.path(match_path, paste0(sample.names, "_singles.fastq"))
for(i in seq_along(fnFs)) {
 if (file.info(fnFs[i])$size == 0) {
  print(paste(fnFs[i], "empty.")) } else {
   print(paste("processing", fnFs[i]))
   commandX <- paste0("/usr/local/home/tobiasgf/bin/sickle-master/sickle se -l 100 -q 28 -x -t sanger -f ",fnFs[i]," -o ",filtFs[i])
   commandY <- paste0("/usr/local/home/tobiasgf/bin/sickle-master/sickle se -l 100 -q 28 -x -t sanger -f ",fnRs[i]," -o ",filtRs[i])
   print(commandX)
   system(commandX)
   print(commandY)
   system(commandY)
   }
}


#Matching the "Sense"-reads.
path <- file.path(main_path, "DADA2_SS/filtered")
fns <- list.files(path)
fastqs <- fns[grepl("fastq$", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]
sample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
match_path <- file.path(path, "matched")
if(!file_test("-d", match_path)) dir.create(match_path)
filtFs <- file.path(match_path, paste0(sample.names, "_F_matched.fastq.gz"))
filtRs <- file.path(match_path, paste0(sample.names, "_R_matched.fastq.gz"))
for(i in seq_along(fnFs)) {
 if (file.info(fnFs[i])$size == 0) {
  print(paste(fnFs[i], "empty.")) } else {
   print(paste("processing", fnFs[i]))
   fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                     matchIDs=TRUE)
  }
}

#Matching the "Anti-Sense"-reads.
path <- file.path(main_path, "DADA2_AS/filtered")
fns <- list.files(path)
fastqs <- fns[grepl("fastq$", fns)]
fastqs <- sort(fastqs)
fnFs <- fastqs[grepl("_F_", fastqs)] # Reverse direction compared to above for the "sense reads". In practice the reads are here complement reversed to be in the same orientation as the "sense" reads.
fnRs <- fastqs[grepl("_R_", fastqs)] # See above
sample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
match_path <- file.path(path, "matched")
if(!file_test("-d", match_path)) dir.create(match_path)
filtFs <- file.path(match_path, paste0(sample.names, "_F_matched.fastq.gz"))
filtRs <- file.path(match_path, paste0(sample.names, "_R_matched.fastq.gz"))
for(i in seq_along(fnFs)) {
 if (file.info(fnFs[i])$size == 0) {
  print(paste(fnFs[i], "empty.")) } else {
   print(paste("processing", fnFs[i]))
   fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                     matchIDs=TRUE)
  }
}


#Processing the set of files containing the forward primer in the R1 reads (the sense reads):
filt_path <- file.path(main_path, "DADA2_SS/filtered/matched") 
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
filt_path <- file.path(main_path, "DADA2_AS/filtered/matched") 
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

#Transpose table, assign names, extract sequences and saving table, for further processing:
trans_total_nochim_sumtable <- as.data.frame(t(total_nochim_sumtable))
#Get DNA sequences
sequences <- row.names(trans_total_nochim_sumtable)
#Assign new rownames
row.names(trans_total_nochim_sumtable) <- paste0("seq",seq.int(nrow(trans_total_nochim_sumtable)))
tbname <- file.path(main_path,"DADA2_nochim.table")
{write.table(trans_total_nochim_sumtable,tbname,sep="\t",col.names = NA, quote=FALSE)}
#Extract OTUs (sequences)
sinkname <- file.path(main_path,"DADA2_nochim.otus")
{
 sink(sinkname)
 for (seqX in seq.int(nrow(trans_total_nochim_sumtable))) {
  header <- paste0(">","seq",seqX,"\n")
  cat(header)
  seqq <- paste0(sequences[seqX],"\n")
  cat(seqq)
 }
 sink()
}

#Define function to extract sequences sample-wise
extrSamDADA2 <- function(my_table) {
 out_path <- file.path(main_path, "DADA2_extracted_samples_nochim")
 if(!file_test("-d", out_path)) dir.create(out_path)
 for (sampleX in seq(1:dim(my_table)[1])){
  sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
  {
   sink(sinkname)
   for (seqX in seq(1:dim(my_table)[2])) {
    if (my_table[sampleX,seqX] > 0) {
     header <- paste0(">",rownames(my_table)[sampleX],";size=",my_table[sampleX,seqX],";","\n")
     cat(header)
     seqq <- paste0(colnames(my_table)[seqX],"\n")
     cat(seqq)
    }
   }
   sink()
  }
 }
}

#Extract samplewise sequences from the non-chimera table using the above function:
extrSamDADA2(total_nochim_sumtable)



#Transpose table, assign names, extract sequences and saving table, for further processing:
trans_total_raw_sumtable <- as.data.frame(t(total_sumtable))
#Get DNA sequences
sequences <- row.names(trans_total_raw_sumtable)
#Assign new rownames
row.names(trans_total_raw_sumtable) <- paste0("seq",seq.int(nrow(trans_total_raw_sumtable)))
tbname <- file.path(main_path,"DADA2_raw.table")
{write.table(trans_total_raw_sumtable,tbname,sep="\t",col.names = NA, quote=FALSE)}
#Extract OTUs (sequences)
sinkname <- file.path(main_path,"DADA2_raw.otus")
{
 sink(sinkname)
 for (seqX in seq.int(nrow(trans_total_raw_sumtable))) {
  header <- paste0(">","seq",seqX,"\n")
  cat(header)
  seqq <- paste0(sequences[seqX],"\n")
  cat(seqq)
 }
 sink()
}

#Define function to extract sequences sample-wise
extrSamDADA2 <- function(my_table) {
 out_path <- file.path(main_path, "DADA2_extracted_samples_raw")
 if(!file_test("-d", out_path)) dir.create(out_path)
 for (sampleX in seq(1:dim(my_table)[1])){
  sinkname <- file.path(out_path, paste0(rownames(my_table)[sampleX],".fas"))
  {
   sink(sinkname)
   for (seqX in seq(1:dim(my_table)[2])) {
    if (my_table[sampleX,seqX] > 0) {
     header <- paste0(">",rownames(my_table)[sampleX],";size=",my_table[sampleX,seqX],";","\n")
     cat(header)
     seqq <- paste0(colnames(my_table)[seqX],"\n")
     cat(seqq)
    }
   }
   sink()
  }
 }
}

#Extract samplewise sequences from the non-chimera table using the above function:
extrSamDADA2(total_sumtable)

getN <- function(x) sum(getUniques(x))
#SSout <- as.data.frame(SSout)
#SSout$sample <-  gsub("_F_matched.fastq.gz","",row.names(SSout))
#SSout <- as.data.frame(SSout)
r1 <- data.frame(sample=names(SSdadaFs), denoisedF=sapply(SSdadaFs, getN),stringsAsFactors=FALSE)
r2 <- data.frame(sample=names(SSdadaRs), denoisedR=sapply(SSdadaRs, getN),stringsAsFactors=FALSE)
r3 <- data.frame(sample=names(SSmergers), merged=sapply(SSmergers, getN),stringsAsFactors=FALSE)
r4 <- data.frame(sample=row.names(seqtab.nochim_SS), nonchim=rowSums(seqtab.nochim_SS),stringsAsFactors=FALSE)

SStrack <- left_join(r1, r2, by='sample') %>% left_join(., r3, by='sample') %>% left_join(., r4, by='sample')

#saveRDS(SStrack, "read_track_SS_RDS")
write.table(SStrack,paste0("read_track_SS_",appendix,".txt"),col.names = NA, quote=FALSE,sep="\t")

#ASout <- as.data.frame(ASout)
#ASout$sample <-  gsub("_F_matched.fastq.gz","",row.names(ASout))
#ASout <- as.data.frame(ASout)
r1 <- data.frame(sample=names(ASdadaFs), denoisedF=sapply(ASdadaFs, getN),stringsAsFactors=FALSE)
r2 <- data.frame(sample=names(ASdadaRs), denoisedR=sapply(ASdadaRs, getN),stringsAsFactors=FALSE)
r3 <- data.frame(sample=names(ASmergers), merged=sapply(ASmergers, getN),stringsAsFactors=FALSE)
r4 <- data.frame(sample=row.names(seqtab.nochim_AS), nonchim=rowSums(seqtab.nochim_AS),stringsAsFactors=FALSE)

AStrack <- left_join(r1, r2, by='sample') %>% left_join(., r3, by='sample') %>% left_join(., r4, by='sample')

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


