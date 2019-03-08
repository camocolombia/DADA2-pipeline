
#Calling libraries
library("ShortRead")
library("DECIPHER")
library("Biostrings")
library("seqinr")
#Reading adjusted taxonomy from JGI
tax_table_csv <- read.csv("E:/DADA2/chapter_1/data/CowPI/otu_table.tab",header = T,sep="\t")
#Reading 16S sequences
#DNA_tab <- read.csv("E:/DADA2/chapter_1/hungate_16S.csv",header = T,sep=",",row.names = 1)


DNA_Seq <- Biostrings::readDNAStringSet("E:/DADA2/chapter_1/data/CowPI/ASVs.fasta",format = "fasta")

#DNA_Seq <- Biostrings::readDNAStringSet("E:/DADA2/chapter_1/db/Hungate_fas_new.fasta",format = "fasta",use.names = T)
to_catch <- as.character(tax_table_csv$OTU)


selected_sequences <- DNA_Seq[to_catch]


Biostrings::writeXStringSet(selected_sequences,
                            file="E:/DADA2/chapter_1/data/CowPI/CowPI_Fil222e.fasta")
