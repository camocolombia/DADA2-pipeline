##################################################
###Chrystian C. Sosa 2018 Chapter 1              #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##################################################

#Calling libraries
library(R.utils)
#Defining paths to be analyzed
# mainDir <- "E:/Dropbox/Dropbox/Paper_PhD"
# chapter <- "chapter_1
mainDir <- "/home/csosa"
chapter <- "chapter_1"

dat_dir <- paste0(mainDir,"/",chapter,"/","data"); if(!file.exists(dat_dir)){dir.create(dat_dir)}
production_dir <- paste0(dat_dir,"/","production"); if(!file.exists(production_dir)){dir.create(production_dir)} # Copy here the production data
#production <- read.csv(paste0(production_dir,"/","production_data.csv"),header=T,sep="|")
out_dir <- paste0(mainDir,"/",chapter,"/","outcomes"); if(!file.exists(out_dir)){dir.create(out_dir)}
graph_dir <- paste0(out_dir,"/","graphics"); if(!file.exists(graph_dir)){dir.create(graph_dir)}
csv_dir <- paste0(out_dir,"/","csv"); if(!file.exists(csv_dir)){dir.create(csv_dir)}
seq_dir <- paste0(dat_dir,"/","sequences"); if(!file.exists(seq_dir)){dir.create(seq_dir)}
seq_raw_dir <- paste0(seq_dir,"/","raw"); if(!file.exists(seq_raw_dir)){dir.create(seq_raw_dir)} #Copy here raw reads

seq_proc_dir <- paste0(seq_dir,"/","proc"); if(!file.exists(seq_proc_dir)){dir.create(seq_proc_dir)}
seq_proc_dir_p <- paste0(seq_proc_dir,"/","paired"); if(!file.exists(seq_proc_dir_p)){dir.create(seq_proc_dir_p)}
seq_proc_dir_up <- paste0(seq_proc_dir,"/","unpaired"); if(!file.exists(seq_proc_dir_up)){dir.create(seq_proc_dir_up)}
seq_filt_dir <- paste0(seq_dir,"/","filter"); if(!file.exists(seq_filt_dir)){dir.create(seq_filt_dir)}
qual_pl_dir <- paste0(graph_dir,"/","qual_plot"); if(!file.exists(qual_pl_dir)){dir.create(qual_pl_dir)}
program_dir <- paste0(mainDir,"/","software")
script_dir <-  paste0(mainDir,"/","scripts")
#List samples sequenced

sample_list <- list.dirs(seq_raw_dir,full.names = F,recursive = F)
sampleNames <- sapply(strsplit(sample_list, "-"), `[`, 1)
# library("RUnit")
# library("survival")
library("ShortRead")
library("knitr")
library("ggplot2")
library("phyloseq")
library("gridExtra")
library("DECIPHER")
library("phangorn")
library("dada2")

#source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")
source(paste0(script_dir,"/","000_trim_adapters.R"))


set.seed(100)



#Uncompress fastq files in gz format in a new folder before trim the primers and adapters.

listing_reads <- function(TRIMMOMATIC){
  
  if(TRIMMOMATIC==F){ 
    cat("using unzip")
    if(length(seq_proc_dir)==0){
      lapply(1:length(sample_list),function(i){
        
        miseq_path <- paste0(seq_raw_dir,"/",sample_list[[i]])
        m_pw <- paste0(miseq_path,"/",list.files(miseq_path))
        m_pw_out <- paste0(seq_proc_dir,"/",list.files(miseq_path))
        m_pw_out <- sub(".gz","",m_pw_out)
        lapply(1:length(m_pw),function(j){
          x <- m_pw[[j]]
          y <- m_pw_out[[j]]
          # #x <- sub(".gz","",x)
          # x <- sub(miseq_path,"",x)
          # x <- sub("/","",x)
          
          R.utils::gunzip(filename=x,destname=y,remove=F)
        })
        
        cat("                                                   ","\n")
        cat("Uncompressed sample: ",i, " of ",length(sample_list),"\n")
        cat("                                                   ","\n")
        cat(round((i/length(sample_list)*100),2)," %","\n")
        cat("                                                   ","\n")
        cat("  DONE!                                            ","\n")
        cat("                                                   ","\n")
      })  
      
    } else {
      cat("NO UNCOMPRESSING WAS REQUIRED!")
    }
  
  } else {
    cat("NO UNCOMPRESSING IS PERFORMED... LISTING FILES TO BE USED IN +TRIMMOMATIC","\n")
    
      lapply(1:length(sample_list),function(i){
        
        miseq_path <- paste0(seq_raw_dir,"/",sample_list[[i]])
        m_pw <- paste0(miseq_path,"/",list.files(miseq_path,pattern=".gz"))
        m_pw_out <- paste0(seq_proc_dir,"/",list.files(miseq_path))
        m_pw_out <- sub(".gz","",m_pw_out)
        lapply(1:length(m_pw),function(j){
          x <- m_pw[[j]]
          y <- m_pw_out[[j]]
          # #x <- sub(".gz","",x)
          # x <- sub(miseq_path,"",x)
          # x <- sub("/","",x)
          
          return(x)
        })
      })
    
  }
}
  
  

x <- unlist(listing_reads(TRIMMOMATIC=T))




# 
# 
# fnFs <- sort(list.files(seq_proc_dir, pattern="_R1_001.fastq"))
# fnRs <- sort(list.files(seq_proc_dir, pattern="_R2_001.fastq"))


trim_folder <- paste0(program_dir,"/","Trimmomatic-0.38")
trim_version <- "trimmomatic-0.38.jar"



trimmomatic_list <- lapply(1:length(sample_list),function(i){
cat("    ","\n");cat(" RUNNING TRIMMOMATIC    ","\n");cat("    ","\n");

miseq_path <- paste0(seq_raw_dir,"/",sample_list[[i]])
paths <- list.files(miseq_path,full.names = T)
names_seq <- list.files(miseq_path,full.names = F)

x <- trim_adapters_function (code = paste0(seq_proc_dir,"/",sampleNames[[i]],"_trim.txt"),trim_folder,trim_version,miseq_path,names_seq,seq_proc_dir_up,seq_proc_dir_p)
cat("    ","\n");cat(round((i/length(sample_list)),2)," %","\n");cat("    ","\n");
return(x)
})








#Sort Forward and reverse reads

fnFs <- sort(list.files(seq_proc_dir, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(seq_proc_dir, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(seq_proc_dir, fnFs)
fnRs <- file.path(seq_proc_dir, fnRs)


F_dir <- paste0(qual_pl_dir,"/","F");if(!file.exists(F_dir)){dir.create(F_dir)}
R_dir <- paste0(qual_pl_dir,"/","R");if(!file.exists(R_dir)){dir.create(R_dir)}
    

#Gathering quality plots for all the samples
lapply(1:length(fnFs),function(i){
cat("Graph ",i," of ",length(fnFs),"\n")
a <- dada2::plotQualityProfile(fnFs[[i]]);b <- plotQualityProfile(fnRs[[i]])
ggsave(paste0(F_dir,"/",as.character(sampleNames[[i]]),"_F.png"),a,width=30,height =30, units = "cm",dpi=100)
ggsave(paste0(R_dir,"/",as.character(sampleNames[[i]]),"_R.png"),b,width=30,height =30, units = "cm",dpi=100)
cat("   ","\n")
})

#Make new filtered sequences
if(!file_test("-d", seq_filt_dir)) dir.create(seq_filt_dir)
filtFs <- file.path(seq_filt_dir, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(seq_filt_dir, paste0(sampleNames, "_R_filt.fastq.gz"))



out <- dada2::filterAndTrim(fwd = fnFs, filt=filtFs, rev = fnRs, filt.rev = filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FAL


out <- as.data.frame(out)
out$coverage <- NA; out$coverage <- (out$reads.out/out$reads.in)*100
out$id <- NA; out$id <- sampleNames
write.csv(out,paste0(seq_dir,"/","fiterAndTrim.csv"),quote = F)


derepFs <- derepFastq(filtFs, verbose=TRUE);names(derepFs) <- sampleNames
saveRDS(derepFs,paste0(seq_dir,"/","F.RDS"));rm(derepFs);gc()

derepRs <- derepFastq(filtRs, verbose=TRUE);names(derepRs) <- sampleNames
saveRDS(derepRs,paste0(seq_dir,"/","R.RDS"));rm(derepRs);gc()
# Name the derep-class objects by the sample names
filtFs <- readRDS(paste0(seq_dir,"/","F.RDS"))
errF <- learnErrors(filtFs, multithread=TRUE)


# x <- fastqPairedFilter(c(fnFs, fnRs), fout=c(filtFs, filtRs), 
# #fastqPairedFilter(c(fnF1, fnR1), fout=c(filtF1, filtR1), 
#                   trimLeft=10, truncLen=c(240, 200), 
#                   maxN=0, maxEE=2,
#                   compress=TRUE, verbose=TRUE)
# 
# 
# 
# testFastqF = system.file("extdata", "sam1F.fastq.gz", package="dada2")
# testFastqR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
# filtFastqF <- tempfile(fileext=".fastq.gz")
# filtFastqR <- tempfile(fileext=".fastq.gz")
# fastqPairedFilter(c(testFastqF, testFastqR), c(filtFastqF, filtFastqR), maxN=0, maxEE=2)
# fastqPairedFilter(c(testFastqF, testFastqR), c(filtFastqF, filtFastqR), trimLeft=c(10, 20),
#                   truncLen=c(240, 200), maxEE=2, rm.phix=TRUE, verbose=TRUE)

