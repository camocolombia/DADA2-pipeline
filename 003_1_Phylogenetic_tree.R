<<<<<<< HEAD
################################################################################
##Chrystian C. Sosa 2018 Chapter 1                                             #
##2018-06-29                                                                   #
##NUIG - Teagasc Athenry                                                       #
##Phylogenetic tree                                                            #
#https://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/#
################################################################################

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
library("ips")
library("phytools")

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
Gblocks_dir <- "E:/DADA2/software/Gblocks/Gblocks_0.91b/Gblocks.exe"

#Assigning fasta file to use in the phylogenetic tree
fasta_file <- paste0(seq_dir,"/","ASVs.fasta")
#Reading ASV fasta file
DNA_Seq <- Biostrings::readDNAStringSet(fasta_file,format = "fasta")
#Orient Nucleotides
DNA_Seq <- DECIPHER::OrientNucleotides(DNA_Seq)
#Aligning sequences using DECIPHER
aligned <- DECIPHER::AlignSeqs(DNA_Seq)
#adjust alignment adjusting Gap Placements
aligned <-DECIPHER::AdjustAlignment(aligned)
#BrowseSeqs(aligned, highlight=0)
#Writing alignment
Biostrings::writeXStringSet(aligned,file=paste0(seq_dir,"/","ASVs_Aligned.fasta"))
#Reading new alignment using ape to be used in ape and Phangorn
dna <- ape::read.dna(paste0(seq_dir,"/","ASVs_Aligned.fasta"), format="fasta")
# counts_Seq <- as.data.frame(DNA_Seq@ranges) #Counting samples
#Deleting low complex regions using gblocks
blocks <- ips::gblocks(dna,b1=0.5,b2=0.5,b3=10,b4=5,b5="h", exec= Gblocks_dir)
#Writing cleaned DNA alignment
write.dna(blocks,format = "fasta",file = paste0(seq_dir,"/","ASVs_GBlocks.fasta"))
#Converting to phyDat object to be used in phangorn
blocks_phyDat <- phyDat(blocks, type = "DNA", levels = NULL)
#Calculating best Nucleotide Substitution model
modelTest_phy <-phangorn::modelTest(blocks_phyDat,multicore=T,mc.cores=2);gc()
saveRDS(modelTest_phy,paste0(seq_dir,"/","modelTest.RDS"))
#Calculating pairwise distancces using the substitution model
dna_dist <- phangorn::dist.ml(blocks_phyDat, model=modelTest_phy$Model[1])
#Calculating pre-existing Neiborghjoining file to optimize
treeNJ <- NJ(dna_dist)
#saving NJ results
phy <- list(blocks_phyDat,dna_dist,treeNJ);saveRDS(phy,paste0(seq_dir,"/","ASVs_phy.RDS"))
#Maximum Likelihood tree
fit = pml(treeNJ, data=blocks_phyDat)
#Optimizing maximum likelihood
fitJC <- optim.pml(fit, model=modelTest_phy$Model[1],pml.control(trace = 0),rearrangement = "NNI")
#logLik(fitJC) Log ML of the tree
#Bootstrap using 100 replicates to validate results
 bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0))
 cnet <- consensusNet(bs, p=0.2)
# plot(cnet, "2D", show.edge.label=TRUE) #Plots
#Joining ML tree with Bootstrap results
 tree <-  plotBS(fit$tree, bs) 
 ss <- list(tree,bs,cnet,fitJC)
#Saving results
 saveRDS(ss,paste0(seq_dir,"/","ASVs_phy_final.RDS"))
 ss <-  readRDS(paste0(seq_dir,"/","ASVs_phy_final.RDS"))
#Saving tree to be used in Phyloseq R package
 write.tree(tree,paste0(seq_dir,"/","pml.tree")) 
=======
################################################################################
##Chrystian C. Sosa 2018 Chapter 1                                             #
##2018-06-29                                                                   #
##NUIG - Teagasc Athenry                                                       #
##Phylogenetic tree                                                            #
#https://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/#
################################################################################

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
library("ips")
library("phytools")

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
Gblocks_dir <- "E:/DADA2/software/Gblocks/Gblocks_0.91b/Gblocks.exe"

#Assigning fasta file to use in the phylogenetic tree
fasta_file <- paste0(seq_dir,"/","ASVs.fasta")
#Reading ASV fasta file
DNA_Seq <- Biostrings::readDNAStringSet(fasta_file,format = "fasta")
#Orient Nucleotides
DNA_Seq <- DECIPHER::OrientNucleotides(DNA_Seq)
#Aligning sequences using DECIPHER
aligned <- DECIPHER::AlignSeqs(DNA_Seq)
#adjust alignment adjusting Gap Placements
aligned <-DECIPHER::AdjustAlignment(aligned)
#BrowseSeqs(aligned, highlight=0)
#Writing alignment
Biostrings::writeXStringSet(aligned,file=paste0(seq_dir,"/","ASVs_Aligned.fasta"))
#Reading new alignment using ape to be used in ape and Phangorn
dna <- ape::read.dna(paste0(seq_dir,"/","ASVs_Aligned.fasta"), format="fasta")
# counts_Seq <- as.data.frame(DNA_Seq@ranges) #Counting samples
#Deleting low complex regions using gblocks
blocks <- ips::gblocks(dna,b1=0.5,b2=0.5,b3=10,b4=5,b5="h", exec= Gblocks_dir)
#Writing cleaned DNA alignment
write.dna(blocks,format = "fasta",file = paste0(seq_dir,"/","ASVs_GBlocks.fasta"))
#Converting to phyDat object to be used in phangorn
blocks_phyDat <- phyDat(blocks, type = "DNA", levels = NULL)
#Calculating best Nucleotide Substitution model
modelTest_phy <-phangorn::modelTest(blocks_phyDat,multicore=T,mc.cores=2);gc()
saveRDS(modelTest_phy,paste0(seq_dir,"/","modelTest.RDS"))
#Calculating pairwise distancces using the substitution model
dna_dist <- phangorn::dist.ml(blocks_phyDat, model=modelTest_phy$Model[1])
#Calculating pre-existing Neiborghjoining file to optimize
treeNJ <- NJ(dna_dist)
#saving NJ results
phy <- list(blocks_phyDat,dna_dist,treeNJ);saveRDS(phy,paste0(seq_dir,"/","ASVs_phy.RDS"))
#Maximum Likelihood tree
fit = pml(treeNJ, data=blocks_phyDat)
#Optimizing maximum likelihood
fitJC <- optim.pml(fit, model=modelTest_phy$Model[1],pml.control(trace = 0),rearrangement = "NNI")
#logLik(fitJC) Log ML of the tree
#Bootstrap using 100 replicates to validate results
 bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0))
 cnet <- consensusNet(bs, p=0.2)
# plot(cnet, "2D", show.edge.label=TRUE) #Plots
#Joining ML tree with Bootstrap results
 tree <-  plotBS(fit$tree, bs) 
 ss <- list(tree,bs,cnet,fitJC)
#Saving results
 saveRDS(ss,paste0(seq_dir,"/","ASVs_phy_final.RDS"))
 ss <-  readRDS(paste0(seq_dir,"/","ASVs_phy_final.RDS"))
#Saving tree to be used in Phyloseq R package
 write.tree(tree,paste0(seq_dir,"/","pml.tree")) 
>>>>>>> 9848fdeb301491061dfcb7d8f57b257d0432bd81
 