##################################################
###Chrystian C. Sosa 2018 Chapter 1              #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##Denoising paired ended FASTQ reads             #
##################################################

#Calling libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")
library("ShortRead")
library("knitr")
library("ggplot2")
library("phyloseq")
library("gridExtra")
library("DECIPHER")
library("phangorn")
library("dada2")
library("ggExtra")

#Defining paths to be analyzed
# mainDir <- "E:/Dropbox/Dropbox/Paper_PhD"
# chapter <- "chapter_1
# mainDir <- "/home/csosa"
# chapter <- "chapter_1"
mainDir <- "E:/DADA2"
chapter <- "chapter_1"

#Defining workspace folders
dat_dir <- paste0(mainDir,"/",chapter,"/","data"); if(!file.exists(dat_dir)){dir.create(dat_dir)}
production_dir <- paste0(dat_dir,"/","production"); if(!file.exists(production_dir)){dir.create(production_dir)} # Copy here the production data
taxadb_dir <- paste0(mainDir,"/",chapter,"/","db"); if(!file.exists(taxadb_dir)){dir.create(taxadb_dir)}
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

#Calling sources
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
        
        cat(" ","\n");cat("Uncompressed sample: ",i, " of ",length(sample_list),"\n");cat(" ","\n")
        cat(round((i/length(sample_list)*100),2)," %","\n")
        cat(" ","\n");cat("  DONE!  ","\n");cat(" ","\n")
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

#Defining primers and parameters to run BBDUK
primers <- "GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT"
trimmomatic_list <- lapply(1:length(sample_list),function(i){
cat("    ","\n");cat(" RUNNING TRIMMING STEP    ","\n");cat("    ","\n");
miseq_path <- paste0(seq_raw_dir,"/",sample_list[[i]])
paths <- list.files(miseq_path,full.names = T,pattern = ".gz")
names_seq <- list.files(miseq_path,full.names = F,pattern = ".gz")
#x <- trim_adapters_function (code = paste0(seq_proc_dir,"/",sampleNames[[i]],"_trim.txt"),trim_folder,trim_version,miseq_path,names_seq,seq_proc_dir_up,seq_proc_dir_p)
trim_adapters_function(program="bbduk",
                       code=paste0(seq_proc_dir,"/",sampleNames[[i]],"_trim.sh"),
                       program_dir, miseq_path, names_seq, seq_proc_dir_up,
                       seq_proc_dir_p, primers, trimq = F, tpe = T,
                       default = F, length_seqParameters="minlength=220 maxlength=280"
                       )
cat("    ","\n");cat(round((i/length(sample_list)*100),2)," %","\n");cat("    ","\n");
return(x)
})

#rm(primers)
##Running bbduck according the machine OS
if(as.character(Sys.info()[1])=="Windows"){
  x <- list.files(seq_proc_dir,"sh",full.names = T)
  x <- lapply(1:length(x),function(i){
    x2 <- read.table(x[[i]],header = F,sep = "?")
    return(x2)
  })
  x <- do.call(rbind,x)
  x <- x$V1
 write.table(x,paste0(seq_dir,"/","To_bbduck.sh"),sep = "",quote = F,row.names = F,col.names = F)
 shell(paste0(seq_dir,"/","To_bbduck.sh"));gc()
  } else {
  cat("Linux operative system... TASK DONE","\n")
}

#Sort Forward and reverse reads
fnFs <- sort(list.files(seq_proc_dir_p, pattern="_R1_paired.fastq"));fnFs_names <-fnFs
fnRs <- sort(list.files(seq_proc_dir_p, pattern="_R2_paired.fastq"));fnRs_names <-fnRs
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(seq_proc_dir_p, fnFs)
fnRs <- file.path(seq_proc_dir_p, fnRs)
# Specify the  path to save quality plots
F_dir <- paste0(qual_pl_dir,"/","F");if(!file.exists(F_dir)){dir.create(F_dir)}
R_dir <- paste0(qual_pl_dir,"/","R");if(!file.exists(R_dir)){dir.create(R_dir)}

##################################################################################################    
#Experimental suggested trim nucleotide possition using the median of the Phred score

#Gathering quality plots for all the samples and the median Phred score for sample
#lapply(1:length(fnFs),function(i){
med_Q <-lapply(1:length(fnFs),function(i){
cat("Graph ",i," of ",length(fnFs),"\n")
a <- dada2::plotQualityProfile(fnFs[[i]]);b <- plotQualityProfile(fnRs[[i]])
a <- a + scale_x_continuous(limits = c(0,250),breaks=seq(0,250,10))
b <- b + scale_x_continuous(limits = c(0,250),breaks=seq(0,250,10))
median_r <- cbind(a$layers[[4]]$data[,c(1,4,6,7)])
median_r <- rbind(median_r,b$layers[[4]]$data[,c(1,4,6,7)])
ggsave(paste0(F_dir,"/",as.character(sampleNames[[i]]),"_F.png"),a,width=30,height =30, units = "cm",dpi=100)
ggsave(paste0(R_dir,"/",as.character(sampleNames[[i]]),"_R.png"),b,width=30,height =30, units = "cm",dpi=100)
cat("   ","\n")
return(median_r)
});gc()

med_Q <- do.call(rbind, med_Q);med_Q2 <- med_Q

#suggested Forward trimming nucleotide possition
med_Q_F <- subset(med_Q[,-3],file %in% fnFs_names)
med_Q_F <- med_Q_F[,-3]
#Phred Score graph for Forward read
q_F_graph <- ggplot(med_Q_F,aes(x=Cycle,y=Q50,colour="black"))+
  #geom_jitter(height = 0.05) +
    geom_line(show.legend = F,colour="black",alpha=0)+ 
  #geom_ribbon( data = med_Q_F , aes( x = Cycle , ymin = min , ymax = max ), alpha = 0.2 , fill='skyblue' )+
  stat_summary(geom="line", fun.y="max",colour="red",size=0.5,linetype="dotted") +
  stat_summary(geom="line", fun.y="min",colour="red",size=0.5,linetype="dotted") +
  stat_summary(geom="line", fun.y="median",size=0.8,colour="green") +
  stat_summary(geom="line", fun.y="mean",size=0.8,colour="orange") +
  scale_x_continuous(limits = c(0,250),breaks=seq(0,250,10))+
  scale_y_continuous(limits = c(10,40),breaks=seq(10,40,5))+
  scale_color_manual(values=c("black"))+
  geom_hline(yintercept=30)+
  xlab("Cycle")+
  ylab("Phred score median")+
  ggtitle("Forward read phred score")
#q_F_graph <- ggMarginal(q_F_graph, type="histogram")
q_F_graph
sug_F <- min(F_med_p2$Cycle[which(F_med_p2$Median_frequency<30)])


# #med_Q_F$Cycle <- as.factor(as.character(med_Q_F$Cycle))
 F_med_p <- reshape2::acast(data = med_Q_F,formula = Cycle ~ Q50, margins = c("file", "Cycle"))
 F_med_p <- F_med_p[-nrow(F_med_p),]
 #as.data.frame.matrix(table(med_Q_F))
 F_med_p2 <- as.data.frame(cbind(as.numeric(as.character(row.names(F_med_p))),
       as.numeric(as.character(colnames(F_med_p)[max.col(F_med_p,ties.method="first")]))))
colnames(F_med_p2) <- c("Cycle","Median_frequency")
 
ggplot(F_med_p2,aes(x=Cycle,y=Median_frequency,colour="black"))+
  geom_line(show.legend = F)+ 
  scale_x_continuous(limits = c(0,250),breaks=seq(0,250,10))+
  scale_color_manual(values=c("black"))
sug_F <- min(F_med_p2$Cycle[which(F_med_p2$Median_frequency<30)])

#suggested Reverse trimming
med_Q_R <- subset(med_Q[,-3],file %in% fnRs_names)
med_Q_R <- med_Q_R[,-3]
#Phred Score graph for Reverse read
q_R_graph <- ggplot(med_Q_R,aes(x=Cycle,y=Q50,colour="black"))+
  #geom_jitter(height = 0.05) +
  geom_line(show.legend = F,colour="black",alpha=0)+ 
  #geom_ribbon( data = med_Q_F , aes( x = Cycle , ymin = min , ymax = max ), alpha = 0.2 , fill='skyblue' )+
  stat_summary(geom="line", fun.y="max",colour="red",size=0.5,linetype="dotted") +
  stat_summary(geom="line", fun.y="min",colour="red",size=0.5,linetype="dotted") +
  stat_summary(geom="line", fun.y="median",size=0.8,colour="green") +
  stat_summary(geom="line", fun.y="mean",size=0.8,colour="orange") +
  scale_x_continuous(limits = c(0,250),breaks=seq(0,250,10))+
  scale_y_continuous(limits = c(10,40),breaks=seq(10,40,5))+
  scale_color_manual(values=c("black"))+
  geom_hline(yintercept=30)+
  xlab("Cycle")+
  ylab("Phred score median")+
  ggtitle("Reverse read phred score")

# #med_Q_F$Cycle <- as.factor(as.character(med_Q_F$Cycle))
R_med_p <- reshape2::acast(data = med_Q_R,formula = Cycle ~ Q50, margins = c("file", "Cycle"))
R_med_p <- R_med_p[-nrow(R_med_p),]
R_med_p2 <- as.data.frame(cbind(as.numeric(as.character(row.names(R_med_p))),
                                as.numeric(as.character(colnames(R_med_p)[max.col(R_med_p,ties.method="first")]))))
colnames(R_med_p2) <- c("Cycle","Median_frequency")

ggplot(R_med_p2,aes(x=Cycle,y=Median_frequency,colour="black"))+
  geom_line(show.legend = F)+ 
  scale_x_continuous(limits = c(0,250),breaks=seq(0,250,10))+
  scale_color_manual(values=c("black"))
sug_R <- min(R_med_p2$Cycle[which(R_med_p2$Median_frequency<30)])
#Make new filtered sequences
if(!file_test("-d", seq_filt_dir)) dir.create(seq_filt_dir)
filtFs <- file.path(seq_filt_dir, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(seq_filt_dir, paste0(sampleNames, "_R_filt.fastq.gz"))


#Filter and Trim step
out <- dada2::filterAndTrim(fwd = fnFs, filt=filtFs, rev = fnRs, filt.rev = filtRs, truncLen=c(230,210),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,minLen=175,
                     #minQ = 10,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FAL
saveRDS(out,paste0(seq_dir,"/","fiterAndTrim.RDS"))#;rm(derepFs);gc()

out <- as.data.frame(out)
out$coverage <- NA; out$coverage <- (out$reads.out/out$reads.in)*100
out$id <- NA; out$id <- sampleNames
out <- cbind(out[,c(ncol(out),1:(ncol(out)-1))])
write.csv(out,paste0(seq_dir,"/","fiterAndTrim.csv"),quote = F,row.names=F)
#out <- read.csv(paste0(seq_dir,"/","fiterAndTrim.csv"))


#Deduplicate pair ended FASTQ files
derepFs <- derepFastq(filtFs, verbose=TRUE);names(derepFs) <- sampleNames
saveRDS(derepFs,paste0(seq_dir,"/","F.RDS"))#;rm(derepFs);gc()

derepRs <- derepFastq(filtRs, verbose=TRUE);names(derepRs) <- sampleNames
saveRDS(derepRs,paste0(seq_dir,"/","R.RDS"))#;rm(derepRs);gc()

# Name the derep-class objects by the sample names
# filtFs <- readRDS(paste0(seq_dir,"/","F.RDS"))

#Training error rate to run dada2
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,paste0(seq_dir,"/","errF.RDS"))#;rm(derepRs);gc()
saveRDS(errR,paste0(seq_dir,"/","errR.RDS"))#;rm(derepRs);gc()

# errF <- readRDS(paste0(seq_dir,"/","errF.RDS"))#;rm(derepRs);gc()
# errR <- readRDS(paste0(seq_dir,"/","errR.RDS"))#;rm(derepRs);gc()

#Plot learn errors plots
er_F_G <- plotErrors(errF, nominalQ=TRUE);ggsave(paste0(graph_dir,"/","LearnErrors_F.pdf"),er_F_G,width=30,height =40, units = "cm",dpi=100)
er_R_G <- plotErrors(errR, nominalQ=TRUE);ggsave(paste0(graph_dir,"/","LearnErrors_R.pdf"),er_R_G,width=30,height =40, units = "cm",dpi=100)

#Run dada model Forward
derepFs <- readRDS(paste0(seq_dir,"/","F.RDS"))
dadaFs <- dada2::dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs,paste0(seq_dir,"/","dadaFs.RDS"))#;rm(derepRs);gc()

#Run dada model Reverse
derepRs <- readRDS(paste0(seq_dir,"/","R.RDS"))
dadaRs <- dada2::dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs,paste0(seq_dir,"/","dadaRs.RDS"))#;rm(derepRs);gc()


#Merging pair ended files
derepFs <- readRDS(paste0(seq_dir,"/","F.RDS"))
derepRs <- readRDS(paste0(seq_dir,"/","R.RDS"))
dadaFs <- readRDS(paste0(seq_dir,"/","dadaFs.RDS"))#;rm(derepRs);gc()
dadaRs <- readRDS(paste0(seq_dir,"/","dadaRs.RDS"))#;rm(derepRs);gc()

mergers <- dada2::mergePairs(dadaF = dadaFs, derepF = derepFs, dadaR = dadaRs, derepR = derepRs, 

                                                          verbose=TRUE,trimOverhang=T)
quantile(mergers$nmatch, probs=seq(0,1,0.05)) 
#mergers <-dada2::mergePairs(dadaF = dadaFs[[1]], derepF = derepFs[[1]], dadaR = dadaRs[[1]], derepR = derepRs[[1]], verbose=TRUE,trimOverhang=T)

saveRDS(mergers,paste0(seq_dir,"/","mergers.RDS"))#;rm(derepRs);gc()
#mergers <- readRDS(paste0(seq_dir,"/","mergers.RDS"))#;rm(derepRs);gc()

#Making a summary table with the sequences counts
seqtab <- makeSequenceTable(mergers)
write.csv(seqtab,paste0(seq_dir,"/","seqtab.csv"),quote = F)

#Remove chimers using de Novo approach
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,paste0(seq_dir,"/","seqtab.nochim.RDS"))#;rm(derepRs);gc()
#seqtab.nochim <- readRDS(paste0(seq_dir,"/","seqtab.nochim.RDS"))
#Explore coverage
#sum(seqtab.nochim)/sum(seqtab)*100

#Save unique merge reads
#uniquesToFasta(seqtab.nochim,paste0(seq_dir,"/","unique.fasta"),ids=paste0("Sequence", seq(length(getSequences(seqtab.nochim)))))

#Writing a count file using the sequences information
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("id","input", "filtered","coverage", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
write.table(track,paste0(seq_dir,"/","track_quality.csv"),quote = F,sep="|",row.names = F)

#Giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}


asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0(seq_dir,"/","ASVs.fasta"))
###############################################################################################
rm(dadaRs,dadaFs,dadaFs,dadaRs,errF,errR,derepRs,derepFs,er_F_G,er_R_G,filtFs,filtRs);gc()
