


##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##ASSIGNING OTUs                                 #
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
ps_object3 <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list_melt",".RDS"))


correlation.table <- microbiome::associate(ps3@otu_table, meta(ps3)[,-c(1:4,10)], 
                                           method = "spearman",
                                           mode = "table", p.adj.threshold = 0.05, n.signif = 1,
                                           p.adj.method = "fdr")
#correlation.table <- correlation.table[which(correlation.table$p.adj<=0.05),]
hm_plot <- microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", 
                 star = "p.adj", p.adj.threshold = 0.05,
                 colours = c("darkred","red","white","blue","darkblue")
                 
) 
# hm_plot <-  hm_plot  +
#   theme(text=element_text(size=60),
#                              legend.text=element_text(size=60),
#                              axis.text.x  = element_text(size=60,colour="black",angle=90),
#                              axis.text.y  = element_text(size=60,colour="black"))
correlation.table
ggsave(paste0(graph_dir,"/","heat_correlation_abund",".pdf"),hm_plot,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)



write.csv(correlation.table,paste0(out_dir,"/","csv/","correlation.table.csv"),quote = F,row.names = F)


ps_object$ps_genus

correlation.table_fam <- microbiome::associate(ps_object2$ps_family2@otu_table, ps_object2$ps_family2@otu_table, 
                                           method = "spearman",
                                           mode = "table", p.adj.threshold = 0.05, n.signif = 1,
                                           p.adj.method = "fdr") #correlation.table <- correlation.table[which(correlation.table$p.adj<=0.05),]
hm_plot_fam <- microbiome::heat(correlation.table_genus, "X1", "X2", fill = "Correlation", 
                            star = "p.adj", p.adj.threshold = 0.05,
                            colours = c("darkred","red","white","blue","darkblue")
                            
) 