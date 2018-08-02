<<<<<<< HEAD
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
s <-  as.character(merge(asv_tab_silva,asv_tab,"sequence",all = T)$seq_id)
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
 s <-  as.character(merge(asv_tab_hungate,asv_tab,"sequence",all = T)$seq_id)
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
asv_tab_blast <- merge(asv_tab,blast_match,by.x="seq_id",by.y="QueryID")
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


asv0 <- merge(asv_tab,asv_tab_silva2,by="seq_id",all=T)
asv1 <- merge(asv0,asv_tab_rdp2,by="seq_id",all=T)
asv2 <- merge(asv1,asv_tab_gg2,by="seq_id",all=T)
asv3 <- merge(asv2,asv_tab_blast2,by="seq_id",all=T)
asv4 <- merge(asv3,asv_tab_hungate2,by="seq_id",all=T)

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

#Suggesting taxonomy
suggest_db_list <-lapply(1:nrow(asv4),function(i){
cat(i,"\n")
x <- asv4[i,c("SILVA_LEVELS","RDP_LEVELS","GG_LEVELS","BLAST_LEVELS","HUNGATE_LEVELS")]

lab_x <- colnames(x)[apply(x,1,which.max)]
lab_x <- sub("_LEVELS","",lab_x)
x2  <- asv4[i,grepl(lab_x, names(asv4))]
x2 <- x2[,-c(ncol(x2))]
colnames(x2) <- levels
x2$suggested_db <- NA
x2$suggested_db <- lab_x
x2$LEVELS_COVERED <- NA
x2$LEVELS_COVERED <- sum(!is.na(x2[,levels]))

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

asv5$LEVELS_COVERED[which(is.na(asv5$final_taxon))] <-7
asv5$genus[which(is.na(asv5$final_taxon) & asv5$LEVELS_COVERED>1)] <- "Eubacterium"
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
########################################################
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

taxonomic_summary_coverage[1,1] <- taxonomic_summary$SILVA_kingdom
taxonomic_summary_coverage[1,2] <- taxonomic_summary$RDP_kingdom
taxonomic_summary_coverage[1,3] <- taxonomic_summary$GG_kingdom
taxonomic_summary_coverage[1,4] <- taxonomic_summary$BLAST_kingdom
taxonomic_summary_coverage[1,5] <- taxonomic_summary$HUNGATE_kingdom
taxonomic_summary_coverage[1,6] <- taxonomic_summary$kingdom
#
taxonomic_summary_coverage[2,1] <- taxonomic_summary$SILVA_phylum
taxonomic_summary_coverage[2,2] <- taxonomic_summary$RDP_phylum
taxonomic_summary_coverage[2,3] <- taxonomic_summary$GG_phylum
taxonomic_summary_coverage[2,4] <- taxonomic_summary$BLAST_phylum
taxonomic_summary_coverage[2,5] <- taxonomic_summary$HUNGATE_phylum
taxonomic_summary_coverage[2,6] <- taxonomic_summary$phylum
#
taxonomic_summary_coverage[3,1] <- taxonomic_summary$SILVA_class
taxonomic_summary_coverage[3,2] <- taxonomic_summary$RDP_class
taxonomic_summary_coverage[3,3] <- taxonomic_summary$GG_class
taxonomic_summary_coverage[3,4] <- taxonomic_summary$BLAST_class
taxonomic_summary_coverage[3,5] <- taxonomic_summary$HUNGATE_class
taxonomic_summary_coverage[3,6] <- taxonomic_summary$class
#
taxonomic_summary_coverage[4,1] <- taxonomic_summary$SILVA_order
taxonomic_summary_coverage[4,2] <- taxonomic_summary$RDP_order
taxonomic_summary_coverage[4,3] <- taxonomic_summary$GG_order
taxonomic_summary_coverage[4,4] <- taxonomic_summary$BLAST_order
taxonomic_summary_coverage[4,5] <- taxonomic_summary$HUNGATE_order
taxonomic_summary_coverage[4,6] <- taxonomic_summary$order
#
taxonomic_summary_coverage[5,1] <- taxonomic_summary$SILVA_family
taxonomic_summary_coverage[5,2] <- taxonomic_summary$RDP_family
taxonomic_summary_coverage[5,3] <- taxonomic_summary$GG_family
taxonomic_summary_coverage[5,4] <- taxonomic_summary$BLAST_family
taxonomic_summary_coverage[5,5] <- taxonomic_summary$HUNGATE_family
taxonomic_summary_coverage[5,6] <- taxonomic_summary$family
#
taxonomic_summary_coverage[6,1] <- taxonomic_summary$SILVA_genus
taxonomic_summary_coverage[6,2] <- taxonomic_summary$RDP_genus
taxonomic_summary_coverage[6,3] <- taxonomic_summary$GG_genus
taxonomic_summary_coverage[6,4] <- taxonomic_summary$BLAST_genus
taxonomic_summary_coverage[6,5] <- taxonomic_summary$HUNGATE_genus
taxonomic_summary_coverage[6,6] <- taxonomic_summary$genus
#
taxonomic_summary_coverage[7,1] <- taxonomic_summary$SILVA_species
taxonomic_summary_coverage[7,2] <- taxonomic_summary$RDP_species
taxonomic_summary_coverage[7,3] <- taxonomic_summary$GG_species
taxonomic_summary_coverage[7,4] <- taxonomic_summary$BLAST_species
taxonomic_summary_coverage[7,5] <- taxonomic_summary$HUNGATE_species
taxonomic_summary_coverage[7,6] <- taxonomic_summary$species
#
taxonomic_summary_coverage[,7] <- levels
#
taxonomic_summary_coverage <- cbind(taxonomic_summary_coverage[,ncol(taxonomic_summary_coverage)],taxonomic_summary_coverage[,c(1:(ncol(taxonomic_summary_coverage)-1))])
colnames(taxonomic_summary_coverage) <- c("LEVEL","SILVA","RDP","GREENGENES","BLAST","HUNGATE","MEDIAN")

write.csv(taxonomic_summary_coverage,paste0(out_dir,"/","csv","/","taxonomic_summary_coverage.csv"),row.names = F,quote = F,na = "")

########################################################
########################################################
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
 

 ########################################################
 ########################################################
 ########################################################











# 
# 
# plot_richness(ps, x="Rumen", measures=c("Shannon", "Simpson"), color="Phase")
# 

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")


plot_ordination(ps.prop, ord.nmds.bray, color="Phase", title="Bray NMDS")

 production$Treatment <- gsub("bottom","Low feed efficiency",production$Treatment); production$Treatment <- gsub("top","High feed efficiency",production$Treatment)
 
 samples.out <- rownames(seqtab.nochim)
 #samples.out <- as.numeric(as.character(samples.out))
 subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
 subject <- as.numeric(as.character(subject));subject2 <- as.data.frame(subject);colnames(subject2) <- "Reference.number"
 prod <- join(subject2,production,by="Reference.number",match="all")
 
 samdf <- data.frame(Subject=subject, Treatment=prod$Treatment,Phase=prod$Phase,Rumen=prod$Rumen.pH.,ADG=prod$ADG,FCE=prod$FCE,ADI=prod$av.dailyintake,TWG=prod$total.wt.Gain)
rownames(samdf) <- samples.out

taxa <- cbind(as.character(asv5$kingdom),
              as.character(asv5$phylum),
              as.character(asv5$class),
              as.character(asv5$order),
              as.character(asv5$family),
              as.character(asv5$genus),
              as.character(asv5$species)
)
# header1 <- matrix(ncol = 3,nrow=nrow(asv5))
# for(i in 1:nrow(asv5)){
#   header1[i,1] <- as.character(paste(taxa[i,],collapse = ";"))
#   header1[i,2] <- as.character(paste(asv5$sequence[i],collapse = ";"))
#   header1[i,3] <- paste(header1[i,1],header1[i,2],sep ='"',collapse = ";")
# };rm(i)
# 
# row.names(taxa) <- as.character(asv5$sequence)
# colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
# 
# ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
#                sample_data(samdf), 
#                tax_table(taxa))
# 
# 
# plot_richness(ps, x="Rumen", measures=c("Shannon", "Simpson"), color="Phase")
# 

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
ord.nmds.CCA <- ordinate(ps, method="RDA")

po_CCA <- plot_ordination(physeq = ps,ordination = ord.nmds.CCA,
                      type="biplot", 
                      color="Treatment",title="Bray NMDS",shape = "Phase")


po <- plot_ordination(physeq = ps.prop,ordination = ord.nmds.bray,
                      #type="samples", 
                      color="Treatment",title="Bray NMDS",shape = "Phase")
po <- po +   scale_shape_manual(values=c(16, 17))
po <- po + scale_size_manual(values=c(9,9))
po <- po + scale_fill_manual("Treatment", values = c("red","blue"))
 po <- po +  theme(panel.background = element_rect(fill = "gray95"),
                   text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
                   axis.text.y  = element_text(size=60,colour="black"))#,

#po <- po + scale_size_manual(values=60)
#po + facet_wrap(~"Treatment", 3)


ggsave(paste0(graph_dir,"/","ordination",".pdf"),po,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)

#https://benjjneb.github.io/dada2/tutorial.html


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Phase", fill="Family") + facet_wrap(~Treatment, scales="free_x")



summarize_phyloseq(ps)
taxonomy <- tax_table(ps)
otu.absolute <- abundances(ps)
otu.relative <- abundances(ps, "compositional")
df <- psmelt(ps)
kable(head(df))

ps.compositional <- transform(ps, "compositional")
g <- global(ps);
g$Reference.number <- NA;g$Reference.number <- row.names(g)
g2 <- join(prod,g,by="Reference.number",match="all")


in_bp <- ggplot(data=g2,aes(x=Treatment,y=evenness_simpson,fill=Phase))+
  geom_boxplot(position = "dodge",outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Simpson values")+ xlab("")+
  ggtitle("")+
  scale_fill_manual("legend", values = c("red","blue"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black")#,
        #legend.title=element_blank(),legend.position="none")
  )
ggsave(paste0(graph_dir,"/","simpson",".pdf"),in_bp,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)


boxplot(g2$evenness_simpson ~ g2$Treatment,fill=g2$Phase)
# Estimate Shannon diversity and add it to the phyloseq object
sample_data(ps)$diversity <- global(ps, index = "shannon")[,1]
plot_regression(diversity ~ Treatment)

# p <- prevalence(ps, detection = 0, sort = TRUE)
 p.core <- core(ps, 
                      detection = .2/100, prevalence = 50/100)
# p.core
# 
# p <- plot_core(transform(p.core, "compositional"), 
#                plot.type = "heatmap", 
#                colours = gray(seq(0,1,length=5)),
#                prevalences = seq(.05, 1, .05), 
#                detections = 10^seq(log10(1e-3), log10(.2), length = 10), 
#                horizontal = TRUE) +
#   xlab("Detection Threshold (Relative Abundance (%))") 
# print(p)  


mp <- plot_composition(p.core, plot.type = "heatmap", transform = "Z", 
                       mar = c(6, 13, 1, 1), sample.sort = "Treatment")
plot_composition(transform(p.core, "compositional"), 
                 plot.type = "barplot", sample.sort = "neatmap")



GPUF <- (ps)


# ig <- make_network(ps, max.dist = 1)
# (p <- plot_network(ig, ps, color = "Treatment", shape = "Phase", 
#                    line_weight = 0.4, label = NULL))


GP <- prune_taxa(taxa_sums(ps) > 0, ps)
human <- get_variable(GP, "Treatment")

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson","evenness_pielou","dominance_gini")
(p <- plot_richness(GP, "Phase", "Treatment", measures=alpha_meas))
p + geom_boxplot(data=p$data, aes(x=Phase, y=value, color=NULL,fill=Treatment), alpha=0.1)



####

ps.compositional <- transform(ps, "compositional")
ps.core <- core(ps.compositional, 
                      detection = .2/100, prevalence = 50/100)

x <- ps.core
quiet(x.ord <- ordinate(x, method = "NMDS", distance = "bray"))
# Pick the projected data (first two columns + metadata)
quiet(proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE))
# Rename the projection axes
names(proj)[1:2] <- paste("Comp", 1:2, sep=".")

p <- plot_landscape(proj[, 1:2], col = proj$Phase,legend = TRUE)
print(p)



log10(ps@tax_tableps)

x <- as.data.frame(ps@otu_table);colnames(x) <- asv5$seq_id
x <- sqrt(x)
y <- ps@sam_data[,-c(1:3)]


correlation.table <- associate(x, y, 
                               method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
kable(head(correlation.table))

correlation.table


heat(correlation.table, "X1", "X2", fill = "Correlation", 
     star = "p.adj", p.adj.threshold = 0.05) 




plot_heatmap(ps)
plot_heatmap(ps, "NMDS", "bray", "Treatment", "Family", low="#000033", high="#FF3300", na.value="white")
heatmap(otu_table(ps))
=======
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
s <-  as.character(merge(asv_tab_silva,asv_tab,"sequence",all = T)$seq_id)
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
 s <-  as.character(merge(asv_tab_hungate,asv_tab,"sequence",all = T)$seq_id)
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
asv_tab_blast <- merge(asv_tab,blast_match,by.x="seq_id",by.y="QueryID")
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


asv0 <- merge(asv_tab,asv_tab_silva2,by="seq_id",all=T)
asv1 <- merge(asv0,asv_tab_rdp2,by="seq_id",all=T)
asv2 <- merge(asv1,asv_tab_gg2,by="seq_id",all=T)
asv3 <- merge(asv2,asv_tab_blast2,by="seq_id",all=T)
asv4 <- merge(asv3,asv_tab_hungate2,by="seq_id",all=T)

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

#Suggesting taxonomy
suggest_db_list <-lapply(1:nrow(asv4),function(i){
cat(i,"\n")
x <- asv4[i,c("SILVA_LEVELS","RDP_LEVELS","GG_LEVELS","BLAST_LEVELS","HUNGATE_LEVELS")]

lab_x <- colnames(x)[apply(x,1,which.max)]
lab_x <- sub("_LEVELS","",lab_x)
x2  <- asv4[i,grepl(lab_x, names(asv4))]
x2 <- x2[,-c(ncol(x2))]
colnames(x2) <- levels
x2$suggested_db <- NA
x2$suggested_db <- lab_x
x2$LEVELS_COVERED <- NA
x2$LEVELS_COVERED <- sum(!is.na(x2[,levels]))

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

asv5$LEVELS_COVERED[which(is.na(asv5$final_taxon))] <-7
asv5$genus[which(is.na(asv5$final_taxon) & asv5$LEVELS_COVERED>1)] <- "Eubacterium"
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
########################################################
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

taxonomic_summary_coverage[1,1] <- taxonomic_summary$SILVA_kingdom
taxonomic_summary_coverage[1,2] <- taxonomic_summary$RDP_kingdom
taxonomic_summary_coverage[1,3] <- taxonomic_summary$GG_kingdom
taxonomic_summary_coverage[1,4] <- taxonomic_summary$BLAST_kingdom
taxonomic_summary_coverage[1,5] <- taxonomic_summary$HUNGATE_kingdom
taxonomic_summary_coverage[1,6] <- taxonomic_summary$kingdom
#
taxonomic_summary_coverage[2,1] <- taxonomic_summary$SILVA_phylum
taxonomic_summary_coverage[2,2] <- taxonomic_summary$RDP_phylum
taxonomic_summary_coverage[2,3] <- taxonomic_summary$GG_phylum
taxonomic_summary_coverage[2,4] <- taxonomic_summary$BLAST_phylum
taxonomic_summary_coverage[2,5] <- taxonomic_summary$HUNGATE_phylum
taxonomic_summary_coverage[2,6] <- taxonomic_summary$phylum
#
taxonomic_summary_coverage[3,1] <- taxonomic_summary$SILVA_class
taxonomic_summary_coverage[3,2] <- taxonomic_summary$RDP_class
taxonomic_summary_coverage[3,3] <- taxonomic_summary$GG_class
taxonomic_summary_coverage[3,4] <- taxonomic_summary$BLAST_class
taxonomic_summary_coverage[3,5] <- taxonomic_summary$HUNGATE_class
taxonomic_summary_coverage[3,6] <- taxonomic_summary$class
#
taxonomic_summary_coverage[4,1] <- taxonomic_summary$SILVA_order
taxonomic_summary_coverage[4,2] <- taxonomic_summary$RDP_order
taxonomic_summary_coverage[4,3] <- taxonomic_summary$GG_order
taxonomic_summary_coverage[4,4] <- taxonomic_summary$BLAST_order
taxonomic_summary_coverage[4,5] <- taxonomic_summary$HUNGATE_order
taxonomic_summary_coverage[4,6] <- taxonomic_summary$order
#
taxonomic_summary_coverage[5,1] <- taxonomic_summary$SILVA_family
taxonomic_summary_coverage[5,2] <- taxonomic_summary$RDP_family
taxonomic_summary_coverage[5,3] <- taxonomic_summary$GG_family
taxonomic_summary_coverage[5,4] <- taxonomic_summary$BLAST_family
taxonomic_summary_coverage[5,5] <- taxonomic_summary$HUNGATE_family
taxonomic_summary_coverage[5,6] <- taxonomic_summary$family
#
taxonomic_summary_coverage[6,1] <- taxonomic_summary$SILVA_genus
taxonomic_summary_coverage[6,2] <- taxonomic_summary$RDP_genus
taxonomic_summary_coverage[6,3] <- taxonomic_summary$GG_genus
taxonomic_summary_coverage[6,4] <- taxonomic_summary$BLAST_genus
taxonomic_summary_coverage[6,5] <- taxonomic_summary$HUNGATE_genus
taxonomic_summary_coverage[6,6] <- taxonomic_summary$genus
#
taxonomic_summary_coverage[7,1] <- taxonomic_summary$SILVA_species
taxonomic_summary_coverage[7,2] <- taxonomic_summary$RDP_species
taxonomic_summary_coverage[7,3] <- taxonomic_summary$GG_species
taxonomic_summary_coverage[7,4] <- taxonomic_summary$BLAST_species
taxonomic_summary_coverage[7,5] <- taxonomic_summary$HUNGATE_species
taxonomic_summary_coverage[7,6] <- taxonomic_summary$species
#
taxonomic_summary_coverage[,7] <- levels
#
taxonomic_summary_coverage <- cbind(taxonomic_summary_coverage[,ncol(taxonomic_summary_coverage)],taxonomic_summary_coverage[,c(1:(ncol(taxonomic_summary_coverage)-1))])
colnames(taxonomic_summary_coverage) <- c("LEVEL","SILVA","RDP","GREENGENES","BLAST","HUNGATE","MEDIAN")

write.csv(taxonomic_summary_coverage,paste0(out_dir,"/","csv","/","taxonomic_summary_coverage.csv"),row.names = F,quote = F,na = "")

########################################################
########################################################
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
 

 ########################################################
 ########################################################
 ########################################################











# 
# 
# plot_richness(ps, x="Rumen", measures=c("Shannon", "Simpson"), color="Phase")
# 

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")


plot_ordination(ps.prop, ord.nmds.bray, color="Phase", title="Bray NMDS")

 production$Treatment <- gsub("bottom","Low feed efficiency",production$Treatment); production$Treatment <- gsub("top","High feed efficiency",production$Treatment)
 
 samples.out <- rownames(seqtab.nochim)
 #samples.out <- as.numeric(as.character(samples.out))
 subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
 subject <- as.numeric(as.character(subject));subject2 <- as.data.frame(subject);colnames(subject2) <- "Reference.number"
 prod <- join(subject2,production,by="Reference.number",match="all")
 
 samdf <- data.frame(Subject=subject, Treatment=prod$Treatment,Phase=prod$Phase,Rumen=prod$Rumen.pH.,ADG=prod$ADG,FCE=prod$FCE,ADI=prod$av.dailyintake,TWG=prod$total.wt.Gain)
rownames(samdf) <- samples.out

taxa <- cbind(as.character(asv5$kingdom),
              as.character(asv5$phylum),
              as.character(asv5$class),
              as.character(asv5$order),
              as.character(asv5$family),
              as.character(asv5$genus),
              as.character(asv5$species)
)
# header1 <- matrix(ncol = 3,nrow=nrow(asv5))
# for(i in 1:nrow(asv5)){
#   header1[i,1] <- as.character(paste(taxa[i,],collapse = ";"))
#   header1[i,2] <- as.character(paste(asv5$sequence[i],collapse = ";"))
#   header1[i,3] <- paste(header1[i,1],header1[i,2],sep ='"',collapse = ";")
# };rm(i)
# 
# row.names(taxa) <- as.character(asv5$sequence)
# colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
# 
# ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
#                sample_data(samdf), 
#                tax_table(taxa))
# 
# 
# plot_richness(ps, x="Rumen", measures=c("Shannon", "Simpson"), color="Phase")
# 

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
ord.nmds.CCA <- ordinate(ps, method="RDA")

po_CCA <- plot_ordination(physeq = ps,ordination = ord.nmds.CCA,
                      type="biplot", 
                      color="Treatment",title="Bray NMDS",shape = "Phase")


po <- plot_ordination(physeq = ps.prop,ordination = ord.nmds.bray,
                      #type="samples", 
                      color="Treatment",title="Bray NMDS",shape = "Phase")
po <- po +   scale_shape_manual(values=c(16, 17))
po <- po + scale_size_manual(values=c(9,9))
po <- po + scale_fill_manual("Treatment", values = c("red","blue"))
 po <- po +  theme(panel.background = element_rect(fill = "gray95"),
                   text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
                   axis.text.y  = element_text(size=60,colour="black"))#,

#po <- po + scale_size_manual(values=60)
#po + facet_wrap(~"Treatment", 3)


ggsave(paste0(graph_dir,"/","ordination",".pdf"),po,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)

#https://benjjneb.github.io/dada2/tutorial.html


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Phase", fill="Family") + facet_wrap(~Treatment, scales="free_x")



summarize_phyloseq(ps)
taxonomy <- tax_table(ps)
otu.absolute <- abundances(ps)
otu.relative <- abundances(ps, "compositional")
df <- psmelt(ps)
kable(head(df))

ps.compositional <- transform(ps, "compositional")
g <- global(ps);
g$Reference.number <- NA;g$Reference.number <- row.names(g)
g2 <- join(prod,g,by="Reference.number",match="all")


in_bp <- ggplot(data=g2,aes(x=Treatment,y=evenness_simpson,fill=Phase))+
  geom_boxplot(position = "dodge",outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Simpson values")+ xlab("")+
  ggtitle("")+
  scale_fill_manual("legend", values = c("red","blue"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black")#,
        #legend.title=element_blank(),legend.position="none")
  )
ggsave(paste0(graph_dir,"/","simpson",".pdf"),in_bp,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)


boxplot(g2$evenness_simpson ~ g2$Treatment,fill=g2$Phase)
# Estimate Shannon diversity and add it to the phyloseq object
sample_data(ps)$diversity <- global(ps, index = "shannon")[,1]
plot_regression(diversity ~ Treatment)

# p <- prevalence(ps, detection = 0, sort = TRUE)
 p.core <- core(ps, 
                      detection = .2/100, prevalence = 50/100)
# p.core
# 
# p <- plot_core(transform(p.core, "compositional"), 
#                plot.type = "heatmap", 
#                colours = gray(seq(0,1,length=5)),
#                prevalences = seq(.05, 1, .05), 
#                detections = 10^seq(log10(1e-3), log10(.2), length = 10), 
#                horizontal = TRUE) +
#   xlab("Detection Threshold (Relative Abundance (%))") 
# print(p)  


mp <- plot_composition(p.core, plot.type = "heatmap", transform = "Z", 
                       mar = c(6, 13, 1, 1), sample.sort = "Treatment")
plot_composition(transform(p.core, "compositional"), 
                 plot.type = "barplot", sample.sort = "neatmap")



GPUF <- (ps)


# ig <- make_network(ps, max.dist = 1)
# (p <- plot_network(ig, ps, color = "Treatment", shape = "Phase", 
#                    line_weight = 0.4, label = NULL))


GP <- prune_taxa(taxa_sums(ps) > 0, ps)
human <- get_variable(GP, "Treatment")

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson","evenness_pielou","dominance_gini")
(p <- plot_richness(GP, "Phase", "Treatment", measures=alpha_meas))
p + geom_boxplot(data=p$data, aes(x=Phase, y=value, color=NULL,fill=Treatment), alpha=0.1)



####

ps.compositional <- transform(ps, "compositional")
ps.core <- core(ps.compositional, 
                      detection = .2/100, prevalence = 50/100)

x <- ps.core
quiet(x.ord <- ordinate(x, method = "NMDS", distance = "bray"))
# Pick the projected data (first two columns + metadata)
quiet(proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE))
# Rename the projection axes
names(proj)[1:2] <- paste("Comp", 1:2, sep=".")

p <- plot_landscape(proj[, 1:2], col = proj$Phase,legend = TRUE)
print(p)



log10(ps@tax_tableps)

x <- as.data.frame(ps@otu_table);colnames(x) <- asv5$seq_id
x <- sqrt(x)
y <- ps@sam_data[,-c(1:3)]


correlation.table <- associate(x, y, 
                               method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
kable(head(correlation.table))

correlation.table


heat(correlation.table, "X1", "X2", fill = "Correlation", 
     star = "p.adj", p.adj.threshold = 0.05) 




plot_heatmap(ps)
plot_heatmap(ps, "NMDS", "bray", "Treatment", "Family", low="#000033", high="#FF3300", na.value="white")
heatmap(otu_table(ps))
>>>>>>> 9848fdeb301491061dfcb7d8f57b257d0432bd81
