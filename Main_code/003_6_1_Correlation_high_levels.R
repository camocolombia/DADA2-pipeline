


##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##Correlate variables using high taxonomic levels#
##################################################

##################################################
#Calling libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("microbiomeSeq")
library("ggplot2")
library("gridExtra")
library("reshape")
library("plyr")
library("phyloseq")
library("microbiome")
#library("microbiomeSeq")
library("vegan")
library("ggvegan")
library("dplyr")
library("ggrepel")
library("ggtree")
library("data.table")
library("DESeq2")
library("BSDA")
library("adegenet")
#install.packages("devtools")
library("devtools")
library("phytools")
library("massMap")
library("structSSI")
library("microbiomeViz");
library("stringr")
library("StructFDR")
library("gtools")
library("agricolae")
library("car")
library("ggdendro")
library("rgl")

#install_github("JiyuanHu/massMap")

##################################################
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



ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))
#ps@sam_data$Treatment <- factor(ps@sam_data$Treatment,levels = c("High feed efficiency","Low feed efficiency"))
ps3 <-  readRDS(paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))
metadata <- as(sample_data(ps3), "data.frame")


ps_object <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list",".RDS"))
ps_object2 <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list_trams",".RDS"))

taxa_df <- data.frame(seqId=row.names(ps3@tax_table),ps3@tax_table)
taxa_df$seqId <- as.character(taxa_df$seqId)


correlation.table_taxa <- lapply(1:(length(levels)-1),function(i){
  
    cat("Proccesing: ",levels[[i]],"\n")

# HFE_Sub <- phyloseq::subset_samples(ps_object2[[i]],Treatment="High feed efficiency")
# LFE_Sub <- phyloseq::subset_samples(ps_object2[[i]],Treatment="Low feed efficiency")
# LI_Sub <- phyloseq::subset_samples(ps_object2[[i]],Phase="liquid")
# SO_Sub <- phyloseq::subset_samples(ps_object2[[i]],Phase="solid")


  HFE_Sub <- as.data.frame((otu_table(subset_samples(ps_object2[[i]], Treatment=="High feed efficiency"))))
  LFE_Sub <- as.data.frame((otu_table(subset_samples(ps_object2[[i]], Treatment=="High feed efficiency"))))
  LI_Sub <- as.data.frame((otu_table(subset_samples(ps_object2[[i]], Phase=="liquid"))))
  SO_Sub <- as.data.frame((otu_table(subset_samples(ps_object2[[i]], Phase=="solid"))))

correlation.table <- microbiome::associate(ps_object2[[i]]@otu_table, meta(ps_object2[[i]])[,-c(1:4,10)], 
                                           method = "spearman",
                                           mode = "table", p.adj.threshold = 0.05, n.signif = 1,
                                           p.adj.method = "fdr")
colnames(correlation.table)[1] <- "ASV"

tax_table_i <- as.data.frame(as(tax_table(ps_object2[[i]]), "matrix"))#as.as.data.frame(matrix(ps_object2[[i]]@tax_table))
tax_table_i$ASV <- row.names(ps_object2[[i]]@tax_table)
tax_table_i <- tax_table_i[,c("ASV",Hmisc::capitalize(levels[i]))]
tax_table_i$mean <- colMeans(ps_object2[[i]]@otu_table)
tax_table_i$sd <- matrixStats::colSds(ps_object2[[i]]@otu_table)
tax_table_i$HFE_mean <- colMeans(HFE_Sub)
tax_table_i$HFE_sd <- matrixStats::colSds(as.matrix(HFE_Sub))
tax_table_i$LFE_mean <- colMeans(LFE_Sub)
tax_table_i$LFE_sd <- matrixStats::colSds(as.matrix(LFE_Sub))
tax_table_i$SO_mean <- colMeans(SO_Sub)
tax_table_i$SO_sd <- matrixStats::colSds(as.matrix(SO_Sub))
tax_table_i$LI_mean <- colMeans(LI_Sub)
tax_table_i$LI_sd <- matrixStats::colSds(as.matrix(LI_Sub))
tax_table_i$status_FE <- NA;#tax_table_i$status_FE <- tax_table_i$HFE_mean > tax_table_i$LFE_mean
tax_table_i$status_PHASE <- NA;#tax_table_i$status_PHASE <- tax_table_i$LI_mean > tax_table_i$SO_mean

for(j in 1:nrow(tax_table_i)){
 if(tax_table_i$HFE_mean[[j]] > tax_table_i$LFE_mean[[j]]){
   tax_table_i$status_FE[[j]] <- "HE" 
   } else if(tax_table_i$HFE_mean[[i]] == tax_table_i$LFE_mean[[j]]) {
   tax_table_i$status_FE[[j]] <- "BOTH"    
   } else {
   tax_table_i$status_FE[[j]] <- "LE" 
   }
  #########
  if(tax_table_i$LI_mean[[j]] > tax_table_i$SO_mean[[j]]){
    tax_table_i$status_PHASE[[j]] <- "LI" 
  } else if(tax_table_i$LI_mean[[j]] == tax_table_i$SO_mean[[j]]) {
    tax_table_i$status_PHASE[[j]] <- "BOTH"    
  } else {
    tax_table_i$status_PHASE[[j]] <- "SO" 
  }
  
};rm(j)

tax_table_i$status_FE[which(tax_table_i$status_FE==T)] <- "HE"
tax_table_i$status_FE[which(tax_table_i$status_FE==F)] <- "LE"

tax_table_i$status_PHASE[which(tax_table_i$status_PHASE==T)] <- "LI"
tax_table_i$status_PHASE[which(tax_table_i$status_PHASE==F)] <- "SO"

tax_table_i$taxon <- levels[[i]]
#tax_table_i <- tax_table_i[which(tax_table_i$p_sign=="*"),]


#tax_table_i <- tax_table_i[complete.cases(tax_table_i[,i]),]

#colnames(tax_table_i)[1] <-    colnames(ps_object2[[i]]@tax_table)[i]

correlation.table <- dplyr::left_join(correlation.table,tax_table_i,"ASV")
correlation.table$p_sign <- NA
correlation.table$p_sign[which(correlation.table$p.adj<0.05)] <- "*"
correlation.table <- correlation.table[which(correlation.table$p_sign=="*"),]
#correlation.table <- correlation.table[,-5]
colnames(correlation.table)[5] <- "Name"
return(correlation.table)

cat("Done: ",levels[[i]],"\n")

})

correlation.table_taxa <- base::do.call(rbind,correlation.table_taxa)
colnames(correlation.table_taxa)[2] <- "Variable"

correlation.table_taxa <- correlation.table_taxa[,c(1:4,19,18,5,6:17)]
write.csv(correlation.table_taxa,paste0(out_dir,"/","csv/","correlation.table_high_taxomomy_full.csv"),quote = F,row.names = F)

#ps_object$ps_genus





# correlation.table_fam <- microbiome::associate(ps_object2$ps_family2@otu_table, ps_object2$ps_family2@otu_table, 
#                                            method = "spearman",
#                                            mode = "table", p.adj.threshold = 0.05, n.signif = 1,
#                                            p.adj.method = "fdr") #correlation.table <- correlation.table[which(correlation.table$p.adj<=0.05),]
# hm_plot_fam <- microbiome::heat(correlation.table_genus, "X1", "X2", fill = "Correlation", 
#                             star = "p.adj", p.adj.threshold = 0.05,
#                             colours = c("darkred","red","white","blue","darkblue")
#                             
# ) 
