#Calling libraries
library("ShortRead")
library("DECIPHER")
library("Biostrings")
#Reading adjusted taxonomy from JGI
tax_table_csv <- read.csv("E:/DADA2/chapter_1/Hungate_taxonomy.csv",header = T,sep="|")
#Reading 16S sequences
DNA_Seq <- Biostrings::readDNAStringSet("E:/DADA2/chapter_1/db/all_rumen_16S_RG.fas",format = "fasta")to_mod <- as.data.frame(names(DNA_Seq));colnames(to_mod) <- "seqID"
#Joining names using the Hungate ids
to_mod2 <- merge(to_mod,tax_table_csv,by.x = "seqID",by.y="Hungate_id_ori",all=F) #try  Hungate_id_ori or Hungate_id
to_mod3 <- paste(to_mod2$seqID,to_mod2$species_strain)
#Writing species assignment fasta dataset
names(DNA_Seq)  <- to_mod3
Biostrings::writeXStringSet(DNA_Seq,
                file="E:/DADA2/chapter_1/db/all_rumen_16S_RG_RF_sp.fas")
#Writing training fasta dataset (Kingdom to Genus)
to_mod3 <- paste(to_mod2$kingdom,to_mod2$phylum,to_mod2$class,to_mod2$order,to_mod2$family,to_mod2$genus,"",sep = ";")
names(DNA_Seq)  <- to_mod3
Biostrings::writeXStringSet(DNA_Seq,
                file="E:/DADA2/chapter_1/db/all_rumen_16S_RG_TrainSet.fas")