##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##ASSIGNING OTUs                                 #
##################################################

#Calling libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")
library("ggplot2")
library("gridExtra")
library("reshape")
library("plyr")
library("phyloseq")
library("microbiome")
library("dada2")
library("tidyverse")
library("stringr")
#Defining paths to be analyzed
# mainDir <- "E:/Dropbox/Dropbox/Paper_PhD"
# chapter <- "chapter_1
# mainDir <- "/home/csosa"
# chapter <- "chapter_1"
mainDir <- "E:/DADA2"
chapter <- "chapter_1"
levels <- c("kingdom","phylum","class","order","family","genus","species")
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
script_dir <-  paste0(mainDir,"/","scripts")#Load de novo cleaned sequences


production_dir <- paste0(dat_dir,"/","production"); if(!file.exists(production_dir)){dir.create(production_dir)}
production <- read.csv(paste0(production_dir,"/","production_data.csv"),header=T,sep="|")

seqtab.nochim <- readRDS(paste0(seq_dir,"/","seqtab.nochim.RDS"))
###########################################################################
#Assigning OTUs using SILVA
if(!file.exists(paste0(seq_dir,"/","taxa_silva.RDS"))){
taxa_silva <- dada2::assignTaxonomy(seqs=seqtab.nochim,
                             refFasta=paste0(taxadb_dir,"/","silva_nr_v132_train_set.fa.gz"), 
                             multithread = T,
                             verbose = FALSE,
                             minBoot=80,
                             tryRC = T)
saveRDS(taxa_silva,paste0(seq_dir,"/","taxa_silva.RDS"))

} else {
  taxa_silva <- readRDS(paste0(seq_dir,"/","taxa_silva.RDS"))  
}

#Assigning Species using SILVA
if(!file.exists(paste0(seq_dir,"/","genus.species_silva.RDS"))){
#genus.species_silva <- assignSpecies(seqtab.nochim,
  genus.species_silva <- dada2::addSpecies(taxa_silva,
                                       paste0(taxadb_dir,"/","silva_species_assignment_v132.fa.gz"),
                                     n=1,
                                     tryRC = T)

saveRDS(genus.species_silva,paste0(seq_dir,"/","genus.species_silva.RDS"))
} else {
  genus.species_silva <- readRDS(paste0(seq_dir,"/","genus.species_silva.RDS"))  
}

###########################################################################
#Assigning OTUs using RDP
if(!file.exists(paste0(seq_dir,"/","taxa_rdp.RDS"))){
taxa_rdp <- dada2::assignTaxonomy(seqs=seqtab.nochim,
                           refFasta=paste0(taxadb_dir,"/","rdp_train_set_16.fa.gz"), 
                           multithread = T,
                           verbose = FALSE,
                           minBoot=80,
                           tryRC = T)
saveRDS(taxa_rdp,paste0(seq_dir,"/","taxa_rdp.RDS"))#;rm(derepRs);gc
} else {
   taxa_rdp <- readRDS(paste0(seq_dir,"/","taxa_rdp.RDS"))
}

#Adding Species using RDP
if(!file.exists(paste0(seq_dir,"/","genus.species_rdp.RDS"))){
#genus.species_rdp <- dada2::assignSpecies(seqtab.nochim,
genus.species_rdp <- dada2::addSpecies(taxa_rdp,
                                   paste0(taxadb_dir,"/","rdp_species_assignment_16.fa.gz"),
                                   n=1,
                                   tryRC = T)
                                   
saveRDS(genus.species_rdp,paste0(seq_dir,"/","genus.species_rdp.RDS"))
} else {
  genus.species_rdp <- readRDS(paste0(seq_dir,"/","genus.species_rdp.RDS"))
}

###########################################################################
#Assigning OTUs using Greengenes
if(!file.exists(paste0(seq_dir,"/","taxa_gg.RDS"))){
taxa_gg <- dada2::assignTaxonomy(seqs=seqtab.nochim,
                          refFasta=paste0(taxadb_dir,"/","gg_13_8_train_set_97.fa.gz"), 
                          multithread = T,
                          verbose = FALSE,
                          minBoot=80,
                          tryRC = T)
saveRDS(taxa_gg,paste0(seq_dir,"/","taxa_gg.RDS"))#;rm(derepRs);gc()
} else {
taxa_gg <- readRDS(paste0(seq_dir,"/","taxa_gg.RDS"))
  
}

###########################################################################
#Assigning OTUs using Hungate1000 Collection
if(!file.exists(paste0(seq_dir,"/","taxa_hungate.RDS"))){
  taxa_hungate <- dada2::assignTaxonomy(seqs=seqtab.nochim,
                                   refFasta=paste0(taxadb_dir,"/","all_rumen_16S_RG_TrainSet.fasta"), 
                                   multithread = T,
                                   verbose = FALSE,
                                   minBoot=80,
                                   tryRC = T)
  saveRDS(taxa_hungate,paste0(seq_dir,"/","taxa_hungate.RDS"))#;rm(derepRs);gc()
} else {
  taxa_hungate <- readRDS(paste0(seq_dir,"/","taxa_hungate.RDS"))
  
}

#Adding Species Hungate1000 Collection
if(!file.exists(paste0(seq_dir,"/","species_hungate.RDS"))){
  #genus.species_rdp <- dada2::assignSpecies(seqtab.nochim,
  genus.species_hungate <- dada2::addSpecies(taxa_hungate,
                                         paste0(taxadb_dir,"/","all_rumen_16S_RG_RF_sp.fasta"),
                                         n=1,
                                         tryRC = T)
  
  saveRDS(genus.species_hungate,paste0(seq_dir,"/","species_hungate.RDS"))
} else {
  genus.species_hungate <- readRDS(paste0(seq_dir,"/","species_hungate.RDS"))
}


########################3333

#Giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}


asv_fasta <- c(rbind(asv_headers, asv_seqs))
# write(asv_fasta, paste0(seq_dir,"/","ASVs.fa"))

#Counting abundances
asv_tab <- t(seqtab.nochim)
sequences_found <- row.names(asv_tab)
row.names(asv_tab) <- sub(">", "", asv_headers)
asv_tab <- as.data.frame(asv_tab)
for(i in 1:ncol(asv_tab)){
  asv_tab[,i] <- as.numeric(as.character(asv_tab[,i])) 
};rm(i)

asv_tab$seq <- NA; asv_tab$seq <- row.names(asv_tab) 
asv_tab <- cbind(asv_tab[,ncol(asv_tab)],sequences_found,asv_tab[,c(1:(ncol(asv_tab)-1))]);colnames(asv_tab)[1:2] <-c("seq_id","sequence")

write.table(asv_tab,  paste0(seq_dir,"/","ASVs_counts.csv"), sep=",", quote=F,row.names = F)

#Getting a SILVA taxonomy
asv_tab_silva <- as.data.frame(unname(genus.species_silva))
asv_tab_silva$sequence <- NA; asv_tab_silva$sequence <- row.names(genus.species_silva)
s <-  as.character(merge(asv_tab_silva,asv_tab,"sequence",all = T,sort = F)$seq_id)
asv_tab_silva$seq_id <- NA; asv_tab_silva$seq_id <- s
asv_tab_silva <- asv_tab_silva[,c(9,1:7)]
colnames(asv_tab_silva) <- c("seq_id",levels)

 write.csv(asv_tab_silva,paste0(out_dir,"/","csv","/","asv_tab_silva.csv"),row.names = F,quote = F,na = "")

#Getting a RDP taxonomy
 asv_tab_rdp <- as.data.frame(unname(genus.species_rdp))
 asv_tab_rdp$sequence <- NA; asv_tab_rdp$sequence <- row.names(genus.species_rdp)
 s <-  as.character(merge(asv_tab_rdp,asv_tab,"sequence",all = T)$seq_id)
 asv_tab_rdp$seq_id <- NA; asv_tab_rdp$seq_id <- s
 asv_tab_rdp <- asv_tab_rdp[,c(9,1:7)]
 colnames(asv_tab_rdp) <- c("seq_id",levels)
 
 write.csv(asv_tab_rdp,paste0(out_dir,"/","csv","/","asv_tab_rdp.csv"),row.names = F,quote = F,na = "")

 #Getting a Hungate taxonomy
 asv_tab_hungate <- as.data.frame(unname(genus.species_hungate))
 asv_tab_hungate$sequence <- NA; asv_tab_hungate$sequence <- row.names(genus.species_hungate)
 s <-  as.character(merge(asv_tab_hungate,asv_tab,"sequence",all = T,sort = F)$seq_id)
 asv_tab_hungate$seq_id <- NA; asv_tab_hungate$seq_id <- s
 asv_tab_hungate <- asv_tab_hungate[,c(9,1:7)]
 colnames(asv_tab_hungate) <- c("seq_id",levels)
 
 write.csv(asv_tab_hungate,paste0(out_dir,"/","csv","/","asv_tab_hungate.csv"),row.names = F,quote = F,na = "")
 
#Getting a Greengenes taxonomy
g <- unname(taxa_gg);colnames(g) <- levels
asv_tab_gg <- cbind(asv_tab,g)
rm(g)

asv_tab_gg$kingdom <- sub("k__","",asv_tab_gg$kingdom)
asv_tab_gg$phylum <- sub("p__","",asv_tab_gg$phylum)
asv_tab_gg$class <- sub("c__","",asv_tab_gg$class)
asv_tab_gg$order <- sub("o__","",asv_tab_gg$order)
asv_tab_gg$family <- sub("f__","",asv_tab_gg$family)
asv_tab_gg$genus <- sub("g__","",asv_tab_gg$genus)
asv_tab_gg$species <- sub("s__","",asv_tab_gg$species)
asv_tab_gg$species <- paste0(asv_tab_gg$genus," ",asv_tab_gg$species)
asv_tab_gg$species[which(asv_tab_gg$species=="NA NA")] <- NA
asv_tab_gg$species <- sub("NA","",asv_tab_gg$species)

for (i in 1:nrow(asv_tab_gg)){
  if(sapply(strsplit(asv_tab_gg$species[[i]], " "), length)==1){
    asv_tab_gg$species[[i]] <- NA
  } else {
    asv_tab_gg$species[[i]] <- asv_tab_gg$species[[i]]
  }
};rm(i)

for (i in 1:nrow(asv_tab_gg)){
  if(sapply(strsplit(asv_tab_gg$species[[i]], " "), length)==1){
    asv_tab_gg$species[[i]] <- NA
  } else {
    asv_tab_gg$species[[i]] <- asv_tab_gg$species[[i]]
  }
};rm(i)

write.csv(asv_tab_gg,paste0(out_dir,"/","csv","/","asv_tab_gg.csv"),row.names = F,quote = F,na = "")

#Loading NCBI results BLAST
blast_match <- read.csv(paste0(seq_dir,"/","16S_BLAST.csv"))
blast_match <- blast_match[,c("QueryID","superkingdom","phylum","class","order","family","genus","species")]
colnames(blast_match)[2] <- "kingdom"
asv_tab_blast <- merge(asv_tab,blast_match,by.x="seq_id",by.y="QueryID",sort = F)
write.csv(asv_tab_blast,paste0(out_dir,"/","csv","/","asv_tab_blast.csv"),row.names = F,quote = F,na = "")

#Merging all taxonomies databases in a unique file
asv_tab_silva2 <- asv_tab_silva[,c("seq_id","kingdom","phylum","class","order","family","genus","species")]
asv_tab_rdp2 <- asv_tab_rdp[,c("seq_id","kingdom","phylum","class","order","family","genus","species")]
asv_tab_gg2 <- asv_tab_gg[,c("seq_id","kingdom","phylum","class","order","family","genus","species")]
asv_tab_blast2 <- asv_tab_blast[,c("seq_id","kingdom","phylum","class","order","family","genus","species")]
asv_tab_hungate2 <- asv_tab_hungate[,c("seq_id","kingdom","phylum","class","order","family","genus","species")]

colnames(asv_tab_silva2) <- c("seq_id","SILVA_kingdom","SILVA_phylum","SILVA_class","SILVA_order","SILVA_family","SILVA_genus","SILVA_species")
colnames(asv_tab_rdp2) <- c("seq_id","RDP_kingdom","RDP_phylum","RDP_class","RDP_order","RDP_family","RDP_genus","RDP_species")
colnames(asv_tab_gg2) <- c("seq_id","GG_kingdom","GG_phylum","GG_class","GG_order","GG_family","GG_genus","GG_species")
colnames(asv_tab_blast2) <- c("seq_id","BLAST_kingdom","BLAST_phylum","BLAST_class","BLAST_order","BLAST_family","BLAST_genus","BLAST_species")
colnames(asv_tab_hungate2) <- c("seq_id","HUNGATE_kingdom","HUNGATE_phylum","HUNGATE_class","HUNGATE_order","HUNGATE_family","HUNGATE_genus","HUNGATE_species")


asv_tab$seq_id <- as.character(asv_tab$seq_id)
asv_tab_silva2$seq_id <- as.character(asv_tab_silva2$seq_id)
asv_tab_gg2$seq_id <- as.character(asv_tab_gg2$seq_id)
asv_tab_blast2$seq_id <- as.character(asv_tab_blast2$seq_id)
asv_tab_hungate2$seq_id <- as.character(asv_tab_hungate2$seq_id)

asv0 <- dplyr::full_join(asv_tab,asv_tab_silva2,by="seq_id")
asv1 <- dplyr::full_join(asv0,asv_tab_rdp2,by="seq_id")
asv2 <- dplyr::full_join(asv1,asv_tab_gg2,by="seq_id")
asv3 <- dplyr::full_join(asv2,asv_tab_blast2,by="seq_id")
asv4 <- dplyr::full_join(asv3,asv_tab_hungate2,by="seq_id")
# asv0 <- merge(asv_tab,asv_tab_silva2,by="seq_id",all=T,sort = F)
# asv1 <- merge(asv0,asv_tab_rdp2,by="seq_id",all=T,sort = F)
# asv2 <- merge(asv1,asv_tab_gg2,by="seq_id",all=T,sort = F)
# asv3 <- merge(asv2,asv_tab_blast2,by="seq_id",all=T,sort = F)
# asv4 <- merge(asv3,asv_tab_hungate2,by="seq_id",all=T,sort = F)

#Detecting taxonomic levels covered by taxonomic db
asv4$SILVA_LEVELS <- NA
asv4$RDP_LEVELS <- NA
asv4$GG_LEVELS <- NA
asv4$BLAST_LEVELS <- NA
asv4$HUNGATE_LEVELS <- NA

tax_labels <- c(colnames(asv_tab_silva2),colnames(asv_tab_rdp2),colnames(asv_tab_gg2),colnames(asv_tab_blast2),colnames(asv_tab_hungate2))
tax_labels <- tax_labels[-grep("seq_id",c(colnames(asv_tab_silva2),colnames(asv_tab_rdp2),colnames(asv_tab_gg2),colnames(asv_tab_blast2),colnames(asv_tab_hungate2)))]
for(i in 1:length(tax_labels)){
  for(j in 1:nrow(asv4)){
    cat("Col: ",tax_labels[[i]]," | row: ",j,"\n")
    if(is.na(asv4[j,tax_labels[[i]]])){
      asv4[j,tax_labels[[i]]] <- NA
    } else if(asv4[j,tax_labels[[i]]]==""){
      asv4[j,tax_labels[[i]]] <- NA
    } else {
      asv4[j,tax_labels[[i]]] <- asv4[j,tax_labels[[i]]]
    }
  };rm(j)
};rm(i)
  
for(i in 1:nrow(asv4)){
asv4$SILVA_LEVELS[[i]] <- sum(!is.na(asv4[i,c("SILVA_kingdom","SILVA_phylum","SILVA_class","SILVA_order","SILVA_family","SILVA_genus","SILVA_species")]))
asv4$RDP_LEVELS[[i]] <- sum(!is.na(asv4[i,c("RDP_kingdom","RDP_phylum","RDP_class","RDP_order","RDP_family","RDP_genus","RDP_species")]))
asv4$GG_LEVELS[[i]] <- sum(!is.na(asv4[i,c("GG_kingdom","GG_phylum","GG_class","GG_order","GG_family","GG_genus","GG_species")]))
asv4$BLAST_LEVELS[[i]] <- sum(!is.na(asv4[i,c("BLAST_kingdom","BLAST_phylum","BLAST_class","BLAST_order","BLAST_family","BLAST_genus","BLAST_species")]))
asv4$HUNGATE_LEVELS[[i]] <- sum(!is.na(asv4[i,c("HUNGATE_kingdom","HUNGATE_phylum","HUNGATE_class","HUNGATE_order","HUNGATE_family","HUNGATE_genus","HUNGATE_species")]))

};rm(i)

write.csv(asv4,paste0(out_dir,"/","csv","/","taxonomy_4DB.csv"),row.names = F,quote = F,na = "")

########################################################
########################################################
########################################################
#i <- 378
#Suggesting taxonomy
suggest_db_list <-lapply(1:nrow(asv4),function(i){
cat(i,"\n")
x <- asv4[i,c("SILVA_LEVELS","RDP_LEVELS","GG_LEVELS","BLAST_LEVELS","HUNGATE_LEVELS")]
x <- x[which(x==max(x))]
lab_x <- colnames(x)
lab_x <- sub("_LEVELS","",lab_x)

#if(sum(x==max(x))>1){
if(any(lab_x=="BLAST" |
   lab_x=="HUNGATE")){
  

  
  if(any(lab_x=="BLAST") &
         any(lab_x=="HUNGATE")
) {   lab_x <- "HUNGATE"
} else if(any(lab_x=="BLAST") &
          !any(lab_x=="HUNGATE")
) {   lab_x <- "BLAST"
} else if(!any(lab_x=="BLAST") &
          any(lab_x=="HUNGATE")
          ) {   lab_x <- "HUNGATE"
}

lab_x <- lab_x
x2  <- asv4[i,grepl(lab_x, names(asv4))]
x2 <- x2[,-c(ncol(x2))]
colnames(x2) <- levels
x2$suggested_db <- NA
x2$suggested_db <- lab_x
x2$LEVELS_COVERED <- NA
x2$LEVELS_COVERED <- sum(!is.na(x2[,levels]))
} else {
  x2  <- asv4[i,grepl(lab_x, names(asv4))]
  x2 <- x2[,-c(ncol(x2))]
  colnames(x2) <- levels
  x2$suggested_db <- NA
  x2$suggested_db <- lab_x[[1]]
  x2$LEVELS_COVERED <- NA
  x2$LEVELS_COVERED <- sum(!is.na(x2[,levels]))
}
return(x2)
})

#Removing taxonomy with kindgom and phylum level
suggest_db_list  <- do.call(rbind,suggest_db_list)
suggest_db_list$suggested_db[which(suggest_db_list$LEVELS_COVERED==0)] <- NA
suggest_db_list$LEVELS_COVERED[which(suggest_db_list$LEVELS_COVERED==0)] <- NA
suggest_db_list$suggested_db[which(suggest_db_list$LEVELS_COVERED==1)] <- NA
suggest_db_list$kingdom[which(suggest_db_list$LEVELS_COVERED==1)] <- NA
suggest_db_list$LEVELS_COVERED[which(suggest_db_list$LEVELS_COVERED==1)] <- NA

#Saving final taxonomy table
asv5 <- cbind(asv4,suggest_db_list)

asv5$final_taxon <- NA
for(i in 1:nrow(asv5)){
  cat(i,"\n")
if(is.na(asv5$LEVELS_COVERED[[i]])){
  asv5$final_taxon[[i]] <- "Unclassified"
} else if(asv5$LEVELS_COVERED[[i]]==1){
  asv5$final_taxon[[i]] <- as.character(asv5$kingdom[[i]])
} else if(asv5$LEVELS_COVERED[[i]]==2){
  asv5$final_taxon[[i]] <- as.character(asv5$phylum[[i]])
} else if(asv5$LEVELS_COVERED[[i]]==3){
  asv5$final_taxon[[i]] <- as.character(asv5$class[[i]])
} else if(asv5$LEVELS_COVERED[[i]]==4){
  asv5$final_taxon[[i]] <- as.character(asv5$order[[i]])
} else if(asv5$LEVELS_COVERED[[i]]==5){
  asv5$final_taxon[[i]] <- as.character(asv5$family[[i]])
} else if(asv5$LEVELS_COVERED[[i]]==6){
  asv5$final_taxon[[i]] <- as.character(asv5$genus[[i]])
} else if(asv5$LEVELS_COVERED[[i]]==7){
  asv5$final_taxon[[i]] <- as.character(asv5$species[[i]])
  } 
};rm(i)

#asv5$LEVELS_COVERED[which(is.na(asv5$final_taxon))] <-7
#asv5$genus[which(is.na(asv5$final_taxon) & asv5$LEVELS_COVERED>1)] <- "Eubacterium"
asv5$final_taxon[which(is.na(asv5$final_taxon) & asv5$LEVELS_COVERED>1)] <-   
asv5$species[which(is.na(asv5$final_taxon) & asv5$LEVELS_COVERED>1)]


asv5$final_taxon[which(asv5$final_taxon=="Unclassified")] <-
paste0(asv5$final_taxon[which(asv5$final_taxon=="Unclassified")],"_",1:length(asv5$final_taxon[which(asv5$final_taxon=="Unclassified")]))

asv5$LEVELS_COVERED[which(asv5$final_taxon=="Mitochondria")] <- 4
asv5$family[which(asv5$final_taxon=="Mitochondria")] <-NA
asv5$final_taxon[which(asv5$final_taxon=="Mitochondria")] <- "Rickettsiales"

asv5$LEVELS_COVERED[which(asv5$final_taxon=="mitochondria")] <- 4
asv5$family[which(asv5$final_taxon=="mitochondria")] <-NA
asv5$final_taxon[which(asv5$final_taxon=="mitochondria")] <- "Rickettsiales"



asv5$final_taxon[which(asv5$final_taxon=="[Eubacterium] cylindroides")] <- "Eubacterium cylindroides"
asv5$final_taxon[which(asv5$final_taxon=="[Eubacterium] eligens")] <- "Eubacterium eligens"
asv5$final_taxon[which(asv5$final_taxon=="[Mogibacteriaceae]")] <- "Mogibacteriaceae"
asv5$final_taxon[which(asv5$final_taxon=="[Paraprevotellaceae]")] <- "Paraprevotellaceae"
asv5$final_taxon[which(asv5$final_taxon=="[Prevotella]")] <- "Prevotella"

asv5$sequence_length <- NA
for(i in 1:nrow(asv5)){
  asv5$sequence_length[[i]] <-nchar(as.character(asv5$sequence[[i]]))
};rm(i)

asv5$species <- as.character(asv5$species)
for(i in 1:nrow(asv5)){
  cat(i,"\n")
  
  if(is.na(asv5$LEVELS_COVERED[[i]])){
    asv5$final_taxon[[i]] <- asv5$final_taxon[[i]]
  } else if (asv5$LEVELS_COVERED[[i]]!=7){
    asv5$final_taxon[[i]] <- asv5$final_taxon[[i]]
  } else if(asv5$LEVELS_COVERED[[i]]==7 & asv5$suggested_db[[i]]=="HUNGATE"){
    asv5$final_taxon[[i]] <- paste(asv5$genus[[i]],asv5$species[[i]],sep = " ")
    asv5$species[[i]] <- as.character(asv5$final_taxon[[i]]) 
  } else if(asv5$LEVELS_COVERED[[i]]==7 & asv5$suggested_db[[i]]=="RDP"){
    asv5$final_taxon[[i]] <- paste(asv5$genus[[i]],asv5$species[[i]],sep = " ")
    asv5$species[[i]] <-  as.character(asv5$final_taxon[[i]]) 
    } else if(asv5$LEVELS_COVERED[[i]]==7 & asv5$suggested_db[[i]]=="SILVA"){
    asv5$final_taxon[[i]] <- paste(asv5$genus[[i]],asv5$species[[i]],sep = " ")
    asv5$species[[i]] <-  as.character(asv5$final_taxon[[i]])  
  } else if(asv5$LEVELS_COVERED[[i]]==7 & asv5$suggested_db[[i]]=="BLAST"){
    asv5$final_taxon[[i]] <- asv5$final_taxon[[i]]
  } else if(asv5$LEVELS_COVERED[[i]]==7 & asv5$suggested_db[[i]]=="GG"){
    asv5$final_taxon[[i]] <- asv5$final_taxon[[i]]}
};rm(i) 

  write.csv(asv5,paste0(out_dir,"/","csv","/","taxonomy_final.csv"),row.names = F,quote = F,na = "")

  ####Taxonomic coverage
saveRDS(asv5,paste0(out_dir,"/","csv","/","taxonomy_final.RDS"))
asv_summary <- as.data.frame(rowSums(asv_tab[,-c(1,2)]))
asv_summary$ID <- row.names(asv_summary)
asv_summary <- asv_summary[,c(2,1)]
colnames(asv_summary) <- c("ID","COUNT")
asv_summary$REL <- NA; asv_summary$REL <- (asv_summary$COUNT/sum(asv_summary$COUNT,na.rm = T)*100)
q1 <- quantile(x =asv_summary$COUNT,c(0.25,0.5,0.75))
asv_summary$QUA <- NA;
asv_summary$QUA[which(asv_summary$COUNT< as.numeric(q1[1]))] <- "OUT_BOTTOM"
asv_summary$QUA[which(asv_summary$COUNT>= as.numeric(q1[1]) & asv_summary$COUNT< as.numeric(q1[2]))] <- "Q1"
asv_summary$QUA[which(asv_summary$COUNT== as.numeric(q1[2]))] <- "Q2"
asv_summary$QUA[which(asv_summary$COUNT> as.numeric(q1[2]) & asv_summary$COUNT< as.numeric(q1[3]))] <- "Q3"
asv_summary$QUA[which(asv_summary$COUNT>= as.numeric(q1[3]))] <- "OUT_TOP"

########################################################
###Calculating coverage per db
taxonomic_summary <- as.data.frame(matrix(ncol=35,nrow=1))
lev_col <- asv4[,(((ncol(asv4)-4)-35):(ncol(asv4)-5))]
colnames(taxonomic_summary) <- colnames(lev_col)

for(i in 1:ncol(lev_col)){
  taxonomic_summary[1,i] <- length((na.omit(lev_col[,i])))/nrow(lev_col)*100
};rm(i,lev_col)

taxonomic_summary$kingdom <- NA; taxonomic_summary$kingdom <-median(taxonomic_summary$SILVA_kingdom,taxonomic_summary$RDP_kingdom,taxonomic_summary$GG_kingdom,taxonomic_summary$BLAST_kingdom,taxonomic_summary$HUNGATE_kingdom,na.rm=T)
taxonomic_summary$phylum <- NA; taxonomic_summary$phylum <-median(taxonomic_summary$SILVA_phylum,taxonomic_summary$RDP_phylum,taxonomic_summary$GG_phylum,taxonomic_summary$BLAST_phylum,taxonomic_summary$HUNGATE_phylum,na.rm=T)
taxonomic_summary$class <- NA; taxonomic_summary$class <-median(taxonomic_summary$SILVA_class,taxonomic_summary$RDP_class,taxonomic_summary$GG_class,taxonomic_summary$BLAST_class,taxonomic_summary$HUNGATE_class,na.rm=T)
taxonomic_summary$order <- NA; taxonomic_summary$order <-median(taxonomic_summary$SILVA_order,taxonomic_summary$RDP_order,taxonomic_summary$GG_order,taxonomic_summary$BLAST_order,taxonomic_summary$HUNGATE_order,na.rm=T)
taxonomic_summary$family <- NA; taxonomic_summary$family <-median(taxonomic_summary$SILVA_family,taxonomic_summary$RDP_family,taxonomic_summary$GG_family,taxonomic_summary$BLAST_family,taxonomic_summary$HUNGATE_family,na.rm=T)
taxonomic_summary$genus <- NA; taxonomic_summary$genus <-median(taxonomic_summary$SILVA_genus,taxonomic_summary$RDP_genus,taxonomic_summary$GG_genus,taxonomic_summary$BLAST_genus,taxonomic_summary$HUNGATE_genus,na.rm=T)
taxonomic_summary$species <- NA; taxonomic_summary$species <-median(taxonomic_summary$SILVA_species,taxonomic_summary$RDP_species,taxonomic_summary$GG_species,taxonomic_summary$BLAST_species,taxonomic_summary$HUNGATE_species,na.rm=T)

taxonomic_summary_coverage <- as.data.frame(matrix(ncol=6,nrow=7))
#kingdom
taxonomic_summary_coverage[1,1] <- taxonomic_summary$SILVA_kingdom
taxonomic_summary_coverage[1,2] <- taxonomic_summary$RDP_kingdom
taxonomic_summary_coverage[1,3] <- taxonomic_summary$GG_kingdom
taxonomic_summary_coverage[1,4] <- taxonomic_summary$BLAST_kingdom
taxonomic_summary_coverage[1,5] <- taxonomic_summary$HUNGATE_kingdom
taxonomic_summary_coverage[1,6] <- taxonomic_summary$kingdom
#phylum
taxonomic_summary_coverage[2,1] <- taxonomic_summary$SILVA_phylum
taxonomic_summary_coverage[2,2] <- taxonomic_summary$RDP_phylum
taxonomic_summary_coverage[2,3] <- taxonomic_summary$GG_phylum
taxonomic_summary_coverage[2,4] <- taxonomic_summary$BLAST_phylum
taxonomic_summary_coverage[2,5] <- taxonomic_summary$HUNGATE_phylum
taxonomic_summary_coverage[2,6] <- taxonomic_summary$phylum
#class
taxonomic_summary_coverage[3,1] <- taxonomic_summary$SILVA_class
taxonomic_summary_coverage[3,2] <- taxonomic_summary$RDP_class
taxonomic_summary_coverage[3,3] <- taxonomic_summary$GG_class
taxonomic_summary_coverage[3,4] <- taxonomic_summary$BLAST_class
taxonomic_summary_coverage[3,5] <- taxonomic_summary$HUNGATE_class
taxonomic_summary_coverage[3,6] <- taxonomic_summary$class
#order
taxonomic_summary_coverage[4,1] <- taxonomic_summary$SILVA_order
taxonomic_summary_coverage[4,2] <- taxonomic_summary$RDP_order
taxonomic_summary_coverage[4,3] <- taxonomic_summary$GG_order
taxonomic_summary_coverage[4,4] <- taxonomic_summary$BLAST_order
taxonomic_summary_coverage[4,5] <- taxonomic_summary$HUNGATE_order
taxonomic_summary_coverage[4,6] <- taxonomic_summary$order
#family
taxonomic_summary_coverage[5,1] <- taxonomic_summary$SILVA_family
taxonomic_summary_coverage[5,2] <- taxonomic_summary$RDP_family
taxonomic_summary_coverage[5,3] <- taxonomic_summary$GG_family
taxonomic_summary_coverage[5,4] <- taxonomic_summary$BLAST_family
taxonomic_summary_coverage[5,5] <- taxonomic_summary$HUNGATE_family
taxonomic_summary_coverage[5,6] <- taxonomic_summary$family
#genus
taxonomic_summary_coverage[6,1] <- taxonomic_summary$SILVA_genus
taxonomic_summary_coverage[6,2] <- taxonomic_summary$RDP_genus
taxonomic_summary_coverage[6,3] <- taxonomic_summary$GG_genus
taxonomic_summary_coverage[6,4] <- taxonomic_summary$BLAST_genus
taxonomic_summary_coverage[6,5] <- taxonomic_summary$HUNGATE_genus
taxonomic_summary_coverage[6,6] <- taxonomic_summary$genus
#species
taxonomic_summary_coverage[7,1] <- taxonomic_summary$SILVA_species
taxonomic_summary_coverage[7,2] <- taxonomic_summary$RDP_species
taxonomic_summary_coverage[7,3] <- taxonomic_summary$GG_species
taxonomic_summary_coverage[7,4] <- taxonomic_summary$BLAST_species
taxonomic_summary_coverage[7,5] <- taxonomic_summary$HUNGATE_species
taxonomic_summary_coverage[7,6] <- taxonomic_summary$species
#
taxonomic_summary_coverage[,7] <- levels
taxonomic_summary_coverage <- cbind(taxonomic_summary_coverage[,ncol(taxonomic_summary_coverage)],taxonomic_summary_coverage[,c(1:(ncol(taxonomic_summary_coverage)-1))])
colnames(taxonomic_summary_coverage) <- c("LEVEL","SILVA","RDP","GREENGENES","BLAST","HUNGATE","MEDIAN")

write.csv(taxonomic_summary_coverage,paste0(out_dir,"/","csv","/","taxonomic_summary_coverage.csv"),row.names = F,quote = F,na = "")

########################################################
#Plotting coverage per taxonomic level
 lev_to_plot <- melt(taxonomic_summary_coverage[,-7])#melt(taxonomic_summary_coverage[-7,-c(5,7)])
 lev_to_plot$LEVEL <- factor(lev_to_plot$LEVEL,levels = levels)
 lev_to_plot<-lev_to_plot[order(
   lev_to_plot$LEVEL == "kingdom",
   lev_to_plot$LEVEL == "phylum",
   lev_to_plot$LEVEL == "class",
   lev_to_plot$LEVEL == "order",
   lev_to_plot$LEVEL == "family",
   lev_to_plot$LEVEL == "genus",
   lev_to_plot$LEVEL == "species",
   decreasing = T),]
 colnames(lev_to_plot) <- c("Level","Database","value")
 coverage_tax_plot <- ggplot(data=lev_to_plot,aes(x=Level,y=value,fill=Database))+
   geom_bar(stat="identity",show.legend = T,position = "dodge")+
   xlab("")+
   ylab("Coverage (%)")+
   # ylim(0,0.1)+
#geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
   theme(text=element_text(size=60),
         legend.text=element_text(size=60),
         axis.text.x  = element_text(size=60,colour="black"),
         axis.text.y  = element_text(size=60,colour="black"))
#scale_y_continuous(limits = c(0,0.8),breaks=seq(0,0.8,0.05))

 ggsave(paste0(graph_dir,"/","Taxonomy_coverage",".pdf"),coverage_tax_plot,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)
 ##################
##Current taxonomic level
 asv5 <-readRDS(paste0(out_dir,"/","csv","/","taxonomy_final.RDS"))
 
 
tax_cov_curr <- as.data.frame.matrix(table(asv5$suggested_db,asv5$LEVELS_COVERED))
colnames(tax_cov_curr) <- levels[-1]
tax_cov_curr$db <- NA; tax_cov_curr$db <- row.names(tax_cov_curr)
tax_cov_curr <- tax_cov_curr[,c(7,1:6)]
tax_cov_curr <- melt(tax_cov_curr)
tax_cov_curr$db <- sub("GG","GREENGENES",tax_cov_curr$db)
  
tax_cov_curr<-tax_cov_curr[order(
  tax_cov_curr$variable == "phylum",
  tax_cov_curr$variable == "class",
  tax_cov_curr$variable == "order",
  tax_cov_curr$variable == "family",
  tax_cov_curr$variable == "genus",
  tax_cov_curr$variable == "species",
  decreasing = T),]

tax_cov_curr<-tax_cov_curr[order(
  tax_cov_curr$db == "SILVA",
  tax_cov_curr$db == "RDP",
  tax_cov_curr$db == "GREENGENES",
  tax_cov_curr$db == "BLAST",
  tax_cov_curr$db == "HUNGATE",
  decreasing = T),]
tax_cov_curr$db <- factor( tax_cov_curr$db,levels = c("SILVA","RDP","GREENGENES","BLAST","HUNGATE"))

colnames(tax_cov_curr) <- c("Database","Level","value")
coverage_tax_curr_plot <- ggplot(data=tax_cov_curr,aes(x=Level,y=value,fill=Database))+
  geom_bar(stat="identity",show.legend = T,position = "dodge")+
  xlab("")+
  ylab("Count")+
  # ylim(0,0.1)+
  #geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black"),
        axis.text.y  = element_text(size=60,colour="black"))

ggsave(paste0(graph_dir,"/","Taxonomy_coverage_current",".pdf"),coverage_tax_curr_plot,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)

##Percentage
tax_cov_curr$value <- tax_cov_curr$value/sum(tax_cov_curr$value)*100
coverage_tax_curr_plot3 <- ggplot(data=tax_cov_curr,aes(x=Database,y=value,fill=Level))+
  geom_bar(stat="identity",show.legend = T,position = "stack")+
  xlab("")+
  ylab("Coverage (%)")+
  # ylim(0,0.1)+
  #geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black"),
        axis.text.y  = element_text(size=60,colour="black"))+
scale_y_continuous(limits = c(0,45),breaks=seq(0,45,5))

ggsave(paste0(graph_dir,"/","Taxonomy_coverage_current_per_db",".pdf"),coverage_tax_curr_plot3,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)
##pre


coverage_tax_curr_plot2 <- ggplot(data=tax_cov_curr,aes(x=Level,y=value,fill=Database))+
  geom_bar(stat="identity",show.legend = T,position = "dodge")+
  xlab("")+
  ylab("Coverage (%)")+
  # ylim(0,0.1)+
  #geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black"),
        axis.text.y  = element_text(size=60,colour="black"))+
  scale_y_continuous(limits = c(0,45),breaks=seq(0,45,5))

ggsave(paste0(graph_dir,"/","Taxonomy_coverage_current_perc",".pdf"),coverage_tax_curr_plot2,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)
