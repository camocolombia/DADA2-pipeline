##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##ASSIGNING OTUs                                 #
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
library("picante")
library("phytools")
#Defining paths to be analyzed
# mainDir <- "E:/Dropbox/Dropbox/Paper_PhD"
# chapter <- "chapter_1
# mainDir <- "/home/csosa"
# chapter <- "chapter_1"
mainDir <- "E:/DADA2"
chapter <- "chapter_1"
levels <- c("kingdom","phylum","class","order","family","genus","species")
ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))

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

#Calling production data
production <- read.csv(paste0(production_dir,"/","production_data.csv"),header=T,sep="|")
production$Treatment <- gsub("bottom","Low feed efficiency",production$Treatment); production$Treatment <- gsub("top","High feed efficiency",production$Treatment)

#################
#Calling functions
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Treatment",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
################
psOrd = subset_taxa(ps, !is.na(Order))
plot_abundance(psOrd,"")


#Loading Phyloseq object
#ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object",".RDS"))
levels(ps@sam_data$Treatment) <- c("Low feed efficiency","High feed efficiency")
sample_data(ps)$UsableReads<-sample_sums(ps)


alpha.1<-as.data.frame(estimate_richness(ps,measures=c("Observed","Shannon","InvSimpson")))
alpha.1.sd<-as(sample_data(ps),'data.frame')
alpha.2<-cbind(alpha.1.sd,alpha.1)
alpha.2$E2<-alpha.2$InvSimpson/alpha.2$Observed
alpha.2$log.reads<-log(alpha.2$UsableReads)




InvSimp_plot <- ggplot(data=alpha.2,aes(x=Treatment,y=InvSimpson))+
  geom_boxplot(aes(fill=Treatment)) +
  stat_boxplot(geom ='errorbar') +
  stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  stat_summary(fun.y=mean, geom="point")+
  guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Simpson´s D")+ xlab("Group")+
  ggtitle("")+
  scale_fill_manual("legend", values = c("blue","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=70),axis.text.x  = element_text(size=70,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=70,colour="black"),
        legend.title=element_blank(),legend.position="none")

ggsave(paste0(graph_dir,"/","SIMPSON_BP.pdf"),InvSimp_plot,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)

Shannon_plot <- ggplot(data=alpha.2,aes(x=Treatment,y=Shannon))+
  geom_boxplot(aes(fill=Treatment)) +
  stat_boxplot(geom ='errorbar') +
  stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  stat_summary(fun.y=mean, geom="point")+
  guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Shannon")+ xlab("Group")+
  ggtitle("")+
  scale_fill_manual("legend", values = c("blue","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=70),axis.text.x  = element_text(size=70,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=70,colour="black"),
        legend.title=element_blank(),legend.position="none")

ggsave(paste0(graph_dir,"/","SHANNON_BP.pdf"),Shannon_plot,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)


###Phylogenetic signal
merge(ps@phy_tree$tip.label,ps@otu_table)

