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













# giving our seq headers more manageable names (ASV_1, ASV_2...)
#  



tax_test <- cbind(as.data.frame(unname(taxa_silva)),as.data.frame(asv_tab))
colnames(tax_test)[1:6] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")




seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
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

