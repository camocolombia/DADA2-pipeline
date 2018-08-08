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

#Calling production data
production <- read.csv(paste0(production_dir,"/","production_data.csv"),header=T,sep="|")
production$Treatment <- gsub("bottom","Low feed efficiency",production$Treatment); production$Treatment <- gsub("top","High feed efficiency",production$Treatment)

#Calling table sequence
seqtab.nochim <- readRDS(paste0(seq_dir,"/","seqtab.nochim.RDS"))
asv5 <- readRDS(paste0(out_dir,"/","csv","/","taxonomy_final.RDS"))
seqid<-cbind(as.character(asv5$seq_id),as.character(asv5$final_taxon))
colnames(seqid) <- c("seqid","final_taxon")
samples.out <- rownames(seqtab.nochim)
#Importing tree
tree_imp <- phyloseq::read_tree(paste0(seq_dir,"/","pml.tree"),errorIfNULL = T);tree2 <- tree_imp
labels_T <- as.data.frame(matrix(ncol=1,nrow=length(tree_imp$tip.label)));colnames(labels_T) <- c("seqid")
labels_T$seqid <- tree_imp$tip.label
labels_T <- merge(labels_T,seqid,by="seqid",all=T);labels_T$final_taxon <- as.character(labels_T$final_taxon)
tree_imp$tip.label <- labels_T$final_taxon
tree_imp$tip.label <- gsub("_\\w+", "", tree_imp$tip.label)

#groupInfo <- split(tree_imp$tip.label, gsub("_\\w+", "", tree_imp$tip.label))
#OTU_t <- groupOTU(tree_imp, groupInfo)
tree_plot <- ggtree(tree_imp, layout='rectangular') + geom_tiplab(size=0.5)#, aes(angle=angle))
#write.tree(tree_imp,paste0(seq_dir,"/","pml_ft.tree")) 
ggsave(paste0(graph_dir,"/","tree_plot",".pdf"),tree_plot,dpi=300,width =90,height=80,units = "cm",scale=0.2,limitsize = FALSE)

###################################################
#Taxonomic coverage per rank

tc1bp <- ggplot(data=asv5[which(!is.na(asv5$LEVELS_COVERED)),],aes(x=factor(LEVELS_COVERED),y=sequence_length,fill=suggested_db))+
  geom_boxplot(position = "dodge",na.rm=T)+
               #,outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Sequence length")+ xlab("Taxonomic rank")+
  ggtitle("")+
  #scale_fill_manual("legend", values = c("red","blue"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black")#,
        #legend.title=element_blank(),legend.position="none"
        )+
  scale_y_continuous(limits = c(220,420),breaks=seq(220,420,20))

ggsave(paste0(graph_dir,"/","taxomomic_depth_reads",".pdf"),tc1bp,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)

####################################################
#samples.out <- as.numeric(as.character(samples.out))
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
subject <- as.numeric(as.character(subject));subject2 <- as.data.frame(subject);colnames(subject2) <- "Reference.number"
prod <- join(subject2,production,by="Reference.number",match="all")


#(my_data) <- c("ADG","FCR","Average daily intake","pH","Total weight gain","Treatment")
samdf <- data.frame(Subject=subject,
                    sampleID=prod$Sample_ID,
                    Treatment=prod$Treatment,
                    Phase=prod$Phase,
                    pH=prod$Rumen.pH.,
                    ADG=prod$ADG,
                    FCR=prod$FCE,
                    Average_daily_intake=prod$av.dailyintake,
                    Total_weight_gain=prod$total.wt.Gain)
rownames(samdf) <- samples.out

taxa <- cbind(as.character(asv5$kingdom),
              as.character(asv5$phylum),
              as.character(asv5$class),
              as.character(asv5$order),
              as.character(asv5$family),
              as.character(asv5$genus),
              as.character(asv5$species)
)

row.names(taxa) <- as.character(asv5$seq_id);#
#row.names(taxa) <- as.character(asv5$sequence)
colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
colnames(seqtab.nochim) <- as.character(asv5$seq_id)

##################################################
#Formatting to phyloseq library
ps <- phyloseq::phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa),#,
               phy_tree(tree2)
)

saveRDS(ps,paste0(csv_dir,"/","Phyloseq_object",".RDS"))
##################################################
#http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html
sdt = data.table(as(sample_data(ps), "data.frame"),
                 TotalReads = sample_sums(ps), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")

#Plotting reads coverage
read_plot <- ggplot(sdt, aes(TotalReads)) + 
  geom_histogram(color = "black", fill = "gray", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank(),
        panel.background = element_rect(fill = "gray95"),
        text=element_text(size=90),
        axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1))

read_plot <- read_plot +   facet_wrap(Phase~Treatment)


ggsave(paste0(graph_dir,"/","seq_reads",".pdf"),read_plot,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)


tdt = data.table(tax_table(ps),
                 TotalCounts = taxa_sums(ps),
                 OTU = taxa_names(ps))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts")

tdt[(TotalCounts <= 0), .N]
tdt[(TotalCounts <= 2), .N]
#Taxa cumulative sum

taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
xlab("Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("") +
  scale_y_continuous(limits = c(0,1510),breaks=seq(0,1510,50)) +
  scale_x_continuous(limits = c(0,700000),breaks=seq(0,700000,10000)) +
  theme(text=element_text(size=20),
        axis.text.x  = element_text(angle=90,size=20,colour="black"),
        axis.text.y  = element_text(size=20,colour="black"))
ggsave(paste0(graph_dir,"/","OTU_counts",".pdf"),pCumSum,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)

##################################################


#Plotting tree
ps_tree2 <- plot_tree(ps,nodelabf=nodeplotboot(),# nodelabf=nodeplotboot(80,0,3),
                      color="Treatment", 
                      size = "Abundance",
                      #label.tips="taxa_names",
                      justify = "yes please", 
                      ladderize="left",
                      plot.margin=0.4#,
                      #shape="Phase"
                      )+
  scale_size_continuous(range = c(1, 3)) +
  scale_color_manual("Treatment", values = c("blue","red"))+
  theme(text=element_text(size=60),
         legend.text=element_text(size=60)) +
          coord_polar(theta="y")    
        

ggsave(paste0(graph_dir,"/","ps_tree",".pdf"),ps_tree2,dpi=300,width =100,height=150,units = "cm",scale=0.8,limitsize = FALSE)

