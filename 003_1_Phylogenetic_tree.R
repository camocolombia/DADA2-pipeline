##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##Phylogenetic tree                              #
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
library("muscle")
library("Biostrings")
library(ips)


# source("https://bioconductor.org/biocLite.R")
# biocLite("XVector")

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

production_dir <- paste0(dat_dir,"/","production"); if(!file.exists(production_dir)){dir.create(production_dir)}
production <- read.csv(paste0(production_dir,"/","production_data.csv"),header=T,sep="|")


Gblocks_dir <- "E:/DADA2/software/Gblocks/Gblocks_0.91b/Gblocks.exe"


fasta_file <- paste0(seq_dir,"/","ASVs.fa")
DNA_Seq <- Biostrings::readDNAStringSet(fasta_file,format = "fasta")
DNA_Seq <- DECIPHER::OrientNucleotides(DNA_Seq)
aligned <- DECIPHER::AlignSeqs(DNA_Seq)
#BrowseSeqs(aligned, highlight=0)
writeXStringSet(aligned,
                file=paste0(seq_dir,"/","ASVs_Aligned.fasta"))


dna <- ape::read.dna(paste0(seq_dir,"/","ASVs_Aligned.fasta"), format="fasta")
# counts_Seq <- as.data.frame(DNA_Seq@ranges)
# nrow(counts_Seq)

#Gblocks
blocks <- gblocks(dna,b1=0.5,b2=0.5,b3=10,b4=5,b5="h", exec= Gblocks_dir)
write.dna(blocks,format = "fasta",file = paste0(seq_dir,"/","ASVs_GBlocks.fasta"))

blocks_phyDat <- phyDat(blocks, type = "DNA", levels = NULL)

modelTest_phy <-phangorn::modelTest(blocks_phyDat);gc()
# 
dna_dist <- dist.ml(blocks_phyDat, model=modelTest_phy$Model[1])
treeNJ <- NJ(dna_dist)
phy <- list(blocks_phyDat,dna_dist,treeNJ)
# 
# saveRDS(phy,paste0(seq_dir,"/","ASVs_phy.RDS"))
# fit = pml(treeNJ, data=blocks_phyDat)
# fitJC <- optim.pml(fit, model=modelTest_phy$Model[1])
# logLik(fitJC)
# bs = bootstrap.pml(fitJC, bs=1000, optNni=TRUE, control = pml.control(trace = 0))

#https://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/
