
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
library("vegan")
library("ggvegan")
library("dplyr")
library("ggrepel")
library("ggtree")
library("data.table")
library("DESeq2")
library("BSDA")
#library("adegenet")
library("gtools")
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


##################################################
##################################################
#Defining functions to be use in the FDR approaches
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

##################################################
##################################################
#Loading Phyloseq object
ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object",".RDS"))

levels(ps@sam_data$Treatment) <- c("Low feed efficiency","High feed efficiency")
sample_data(ps)$UsableReads<-sample_sums(ps)


##################################################
##################################################
#Exploring differences in libraries per group

HL <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="High feed efficiency" & Phase=="liquid"))))
HS <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="High feed efficiency" & Phase=="solid"))))
LL <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="Low feed efficiency" & Phase=="liquid"))))
LS <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="Low feed efficiency" & Phase=="solid"))))


##Summarizing libraries size


lib_df <- data.frame(samples = c(names(HL),names(HS),names(LL),names(LS)),
                        sampleType = c(rep("liquid",length(HL)),rep("solid",length(HS)),rep("liquid",length(LL)),rep("solid",length(LS))),
                        Treatment = c(rep("High feed efficiency",length(HL)),rep("High feed efficiency",length(HS)),rep("Low feed efficiency",length(LL)),rep("Low feed efficiency",length(LS))),  
                          libSize = as.numeric(c(colSums(HL),colSums(HS),colSums(LL),colSums(LS)))
                        )
lib_df$libSize_rel <- NA; lib_df$libSize_rel <- lib_df$libSize/sum(lib_df$libSize)

#lib__df_df2 <- lib_df[which(lib_df$libSize>=100000),]

####


#ggsave(paste0(graph_dir,"/","Production_BP_scaled",".pdf"),p1bp,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)
#p1bp
##################################################
##################################################
#Prevalence cleaning up approach



#ps = 
ps <- subset_taxa(ps, Phylum!=is.na(Phylum))

prev0 = apply(X = otu_table(ps),
              MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)]
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
ps1 = prune_taxa((prev0 > prevalenceThreshold), ps)
ps1
# Filter entries with unidentified Phylum.
ps2 = subset_taxa(ps1, Phylum %in% names(keepPhyla))
saveRDS(ps2,paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))
ps <- ps2
####prevalence graph

#Prevalence
prev_plot <- ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = "black")) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7,show.legend = F) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum)+ 
  scale_color_manual("Phylum", values = c("black"))

ggsave(paste0(graph_dir,"/","ps_prev",".pdf"),prev_plot,dpi=300,width =100,height=150,units = "cm",scale=0.8,limitsize = FALSE)


############after prevalence treatment


HL <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="High feed efficiency" & Phase=="liquid"))))
HS <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="High feed efficiency" & Phase=="solid"))))
LL <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="Low feed efficiency" & Phase=="liquid"))))
LS <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="Low feed efficiency" & Phase=="solid"))))

lib_df$libSize_T <- NA;lib_df$libSize_T <- as.numeric(c(colSums(HL),colSums(HS),colSums(LL),colSums(LS)))
lib_df$libSize_T_rel <- NA; lib_df$libSize_T_rel <- lib_df$libSize_T/sum(lib_df$libSize_T)

rm(HL,HS,LL,LS)
lib_df2 <- melt(lib_df)
lib_df2 <- lib_df2[which(lib_df2$variable=="libSize" | lib_df2$variable=="libSize_T"),]
lib_df2$variable <- gsub("libSize","Before filter",lib_df2$variable)
lib_df2$variable <- gsub("Before filter_T","After filter",lib_df2$variable)

lib_df2$variable[which(lib_df2$variable=="libSize")] <- "Before filter"
lib_df2$variable[which(lib_df2$variable=="libSize_T")] <- "After filter"
lib_df2$Type <- NA;lib_df2$Type <- paste0(lib_df2$Treatment,"-",lib_df2$sampleType)

lib_size <- ggplot(data=lib_df2,aes(x=Type,y=value,fill=variable))+
  geom_boxplot(position = "dodge",outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Scaled values")+ xlab("")+
  ggtitle("")+
    scale_fill_manual("", values = c("blue","lightblue","orange","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black",angle = 45, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black")+
          scale_y_continuous(limits = c(50000,140000),breaks=seq(50000,140000,5000))
        
          #,
        #legend.title=element_blank(),legend.position="none"
  )#+ 
  #coord_flip()
ggsave(paste0(graph_dir,"/","lib_size_plot",".pdf"),lib_size,dpi=300,width =100,height=150,units = "cm",scale=0.8,limitsize = FALSE)

#size difference
#lib_df$libSize - lib_df$libSize_T
#lib_size

ps3 <- transform_sample_counts(ps,function(x) {x/sum(x)} )
saveRDS(ps3,paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))

#Aggregating taxa per taxonomic level
ps_kingdom <- aggregate_taxa(ps, 'Kingdom');ps_kingdom2 <- transform_sample_counts(ps_kingdom,function(x) {x/sum(x)} ) 
ps_phylum <- aggregate_taxa(ps, 'Phylum');ps_phylum2 <- transform_sample_counts(ps_phylum,function(x) {x/sum(x)} ) 
ps_class <- aggregate_taxa(ps, 'Class');ps_class2 <- transform_sample_counts(ps_class,function(x) {x/sum(x)} ) 
ps_order <- aggregate_taxa(ps, 'Order');ps_order2 <- transform_sample_counts(ps_order,function(x) {x/sum(x)} ) 
ps_family <- aggregate_taxa(ps, 'Family');ps_family2 <- transform_sample_counts(ps_family,function(x) {x/sum(x)} ) 
ps_genus <- aggregate_taxa(ps, 'Genus');ps_genus2 <- transform_sample_counts(ps_genus,function(x) {x/sum(x)} ) 


ps_object <- list(ps_kingdom = ps_kingdom,
                  ps_phylum = ps_phylum,
                  ps_class = ps_class,
                  ps_order = ps_order,
                  ps_family = ps_family,
                  ps_genus = ps_genus
                  )

saveRDS(ps_object,paste0(csv_dir,"/","Phyloseq_taxa_list",".RDS"))

ps_object2 <- list(ps_kingdom2 = ps_kingdom2,
                  ps_phylum2 = ps_phylum2,
                  ps_class2 = ps_class2,
                  ps_order2 = ps_order2,
                  ps_family2 = ps_family2,
                  ps_genus2 = ps_genus2
)
saveRDS(ps_object2,paste0(csv_dir,"/","Phyloseq_taxa_list_trams",".RDS"))s

ps_object3 <- list(
pserie_kingdom = psmelt(ps_kingdom2) %>%   arrange(Kingdom),
pserie_phylum = psmelt(ps_phylum2) %>%   arrange(Phylum),
pserie_class = psmelt(ps_class2) %>%   arrange(Class),
pserie_order = psmelt(ps_order2) %>%   arrange(Order),
pserie_family = psmelt(ps_family2) %>%   arrange(Family),
pserie_genus = psmelt(ps_genus2) %>%   arrange(Genus)
)
ps_object3$pserie_kingdom$tax_level <- "Kingdom"
ps_object3$pserie_phylum$tax_level <- "Phylum"
ps_object3$pserie_class$tax_level <- "Class"
ps_object3$pserie_order$tax_level <- "Order"
ps_object3$pserie_family$tax_level <- "Family"
ps_object3$pserie_genus$tax_level <- "Genus"

saveRDS(ps_object3,paste0(csv_dir,"/","Phyloseq_taxa_list_melt",".RDS"))

rm(ps_kingdom,ps_kingdom2,ps_phylum,ps_phylum2,ps_class,ps_class2,
   ps_order,ps_order2,ps_family,ps_family2,ps_genus,ps_genus2,prev0,prevalenceThreshold)
##################################################
##################################################
#####DESEQ2

#Converting to DESeq2object
diagdds = phyloseq_to_deseq2(ps, ~ Treatment)
#Applying geometric means
geoMeans = apply(counts(diagdds), 1, gm_mean)
#Estimating size factors
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds <- estimateDispersions(diagdds)
#varianceStabilizingTransformation(diagdds, blind = TRUE, fitType = "parametric")
abundds2 <- getVarianceStabilizedData(diagdds)
diagdds2 <- DESeq2::DESeq(diagdds,fitType="local") 
res2 <- as.data.frame(results(diagdds2)); 
res2 <- res2[,c("baseMean", "log2FoldChange", "lfcSE", "padj")]

sigtab = res2[(res2$padj <= 0.05), ]; sigtab <- sigtab[complete.cases(sigtab),]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
#head(sigtab)

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
sigtabgen = subset(sigtab, !is.na(Genus))
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



DESeq_list <- list( abundds2 = abundds2,
                 diagdds2= diagdds2,
                 sigtab = sigtab,
                 posigtab = posigtab,
                 sigtabgen = sigtabgen
)
saveRDS(DESeq_list,paste0(csv_dir,"/","DESeq_list",".RDS"))

 ###Check DESeq2 VST Transformation
# 
# abund_sums <- rbind(data.frame(sum = colSums(abundds2),
#                                sample = colnames(abundds2),
#                                type = "DESeq2"),
#                     data.frame(sum = rowSums(otu_table(ps3)),
#                                sample = rownames(otu_table(ps3)),
#                                type = "log(1 + x)"))
# ggplot(abund_sums) +
#   geom_histogram(aes(x = sum), binwidth = 50) +
#   facet_grid(type ~ .) +
#   xlab("Total abundance within sample")


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
        legend.text=element_text(size=60)) #+
  #coord_polar(theta="y")    


ggsave(paste0(graph_dir,"/","ps_tree",".pdf"),ps_tree2,dpi=300,width =100,height=150,units = "cm",scale=0.8,limitsize = FALSE)

rm(geoMeans,keepPhyla,abundds2,res2,sigtab,sigtabgen,posigtab,ps1,ps2,ps3,ps,x)
