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
boxplot(asv5$sequence_length ~ asv5$LEVELS_COVERED,ylab="Sequence length",xlab="Taxonomic rank")

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

#Formatting to phyloseq library
ps <- phyloseq::phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa),#,
               phy_tree(tree2)
)

saveRDS(ps,paste0(csv_dir,"/","Phyloseq_object",".RDS"))

sample_sum_df <- data.frame(sum = sample_sums(ps))

#Plotting reads coverage
read_plot <- ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "gray", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank(),
        panel.background = element_rect(fill = "gray95"),
              text=element_text(size=90),
              axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1)
        )#,


ggsave(paste0(graph_dir,"/","seq_reads",".pdf"),read_plot,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)


treat_phase <- as.data.frame(cbind(samdf$Subject,paste(samdf$Treatment,samdf$Phase,sep = "_")))
colnames(treat_phase) <- c("Subject","Treat_Phase")
treat_phase$Subject <- as.character(treat_phase$Subject)
ps_otu_hclust$labels <-merge(ps_otu_hclust$labels,treat_phase,by.x="x",by.y="Subject")[,2]
plot(ps_otu_hclust)
# #Stacked barplots
# 
# # Set colors for plotting
# phylum_colors <- c(
#   "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
#   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
#   "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861","red","green","orange","lightgreen","darkgrey","blue","violet"
# )

#Subsetting phyloseq objects per taxonomic range
#Phylum
pserie_phylum <- ps %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Phylum)     
#Family
pserie_family <- ps %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort d
#Genus
pserie_genus <- ps %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort d
#Species
pserie_sp <- ps %>%
  tax_glom(taxrank = "Species") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Species)  

# Plot Phylum
phylum_graph <- ggplot(pserie_sp, aes(x = Treatment, y = Abundance, fill = reorder(Phylum, -Abundance))) + 
  facet_grid(Kingdom~.) +
  #facet_grid(Treatment~.) +
  geom_bar(stat = "identity",position = "dodge") +
 # scale_fill_manual(values = phylum_colors) +
  # scale_x_discrete(
  #   breaks = c("7/8", "8/4", "9/2", "10/6"),
  #   labels = c("Jul", "Aug", "Sep", "Oct"), 
  #   drop = FALSE
  # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        text=element_text(size=40),
        legend.text=element_text(size=40),
        axis.text.y  = element_text(size=49,colour="black")) +  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 5%) \n") +
  ggtitle("Phylum Composition") +
  labs(fill="Phylum")

ggsave(paste0(graph_dir,"/","phylum_graph",".pdf"),phylum_graph,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)

# Plot family
family_graph <- ggplot(pserie_sp, aes(x = Treatment, y = Abundance, fill = reorder(Family, -Abundance))) + 
  facet_grid(Kingdom~.) +
  geom_bar(stat = "identity",position = "dodge") +
  # scale_fill_manual(values = phylum_colors) +
  # scale_x_discrete(
  #   breaks = c("7/8", "8/4", "9/2", "10/6"),
  #   labels = c("Jul", "Aug", "Sep", "Oct"), 
  #   drop = FALSE
  # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        text=element_text(size=40),
        legend.text=element_text(size=40),
        axis.text.y  = element_text(size=49,colour="black")) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Family > 5%) \n") +
  ggtitle("Family Composition") +
  labs(fill="Family")

ggsave(paste0(graph_dir,"/","family_graph",".pdf"),family_graph,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)

# Plot genus
genus_graph <- ggplot(pserie_sp, aes(x = Treatment, y = Abundance, fill = reorder(Genus, -Abundance))) + 
  facet_grid(Kingdom~.) +
  geom_bar(stat = "identity",position = "dodge") +
  # scale_fill_manual(values = phylum_colors) +
  # scale_x_discrete(
  #   breaks = c("7/8", "8/4", "9/2", "10/6"),
  #   labels = c("Jul", "Aug", "Sep", "Oct"), 
  #   drop = FALSE
  # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        text=element_text(size=40),
        legend.text=element_text(size=40),
        axis.text.y  = element_text(size=49,colour="black")) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 5%) \n") +
  ggtitle("Genus Composition") +
  labs(fill="Genus")

ggsave(paste0(graph_dir,"/","genus_graph",".pdf"),genus_graph,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)

# Plot species
sp_graph <- ggplot(pserie_sp, aes(x = Treatment, y = Abundance, fill = reorder(Species, -Abundance))) + 
  facet_grid(Kingdom~.) +
  #facet_grid(Treatment~.) +
  geom_bar(stat = "identity",position = "dodge") +
  # scale_fill_manual(values = phylum_colors) +
  # scale_x_discrete(
  #   breaks = c("7/8", "8/4", "9/2", "10/6"),
  #   labels = c("Jul", "Aug", "Sep", "Oct"), 
  #   drop = FALSE
  # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        text=element_text(size=40),
        legend.text=element_text(size=40),
        axis.text.y  = element_text(size=49,colour="black")) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Species > 5%) \n") +
  ggtitle("Species Composition") +
  labs(fill="Species")

ggsave(paste0(graph_dir,"/","sp_graph",".pdf"),sp_graph,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)





  # plot(pserie_sp$Abundance,pserie_sp$FCR,col=pserie_sp$Treatment,pch=18)
# 
# cor.test(pserie_sp$Abundance,pserie_sp$FCR)

# Scale reads to even depth 
# ps_scale <- ps %>%
#   scale_reads(round = "round") 



# pseq.subset_liquid <- subset_samples(ps, Phase == "liquid")
# s <- rownames(subset(meta(ps), Phase == "liquid"))
# pseq.subset <- prune_samples(s, ps)
# 
# 
# 
# 
# 
# 
# n <- ntaxa(ps)
# topx <- top_taxa(ps, n = 20)
# ranks <- rank_names(ps)  # Taxonomic levels
# taxa  <- taxa(ps)        # Taxa names at the analysed level
# 
# 
# 
# 
require(ade4)
  ###Core microbiome
ps_otu <- as.data.frame(otu_table(ps))
for(i in 1:ncol(ps_otu)){
ps_otu[,i][which(ps_otu[,i]>0)]<-1
};rm(i)
ps_otu2 <- ps_otu
#ps_otu_dist <- dist(ps_otu,method = "binary")#vegdist(ps_otu,method = "jaccard",binary = T,na.rm = T)
ps_otu_dist <- dist.binary(df = ps_otu,method = 1)#vegdist(ps_otu,method = "jaccard",binary = T,na.rm = T
ps_otu$sample <- NA;ps_otu$sample <- row.names(ps_otu)
ss <- merge(ps_otu,prod,by.x="sample",by.y="Reference.number",all=T)
# ps_otu$Treatment <- NA;
# ps_otu$Treatment <- ss$Treatment.y

#ps_otu$Treatment <- 
 # merge(ps_otu,prod,by.x="sample",by.y="Reference.number")

#ps_otu <- as.data.frame(ps_otu[,c(ncol(ps_otu)-1,(ncol(ps_otu)),1:(ncol(ps_otu)-2))])
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
ht <- as.data.frame(t(ss[,c(1:(ncol(ps@otu_table))+1)][which(ss$Treatment=="High feed efficiency"),]))
hrS <- rowSums(ht)
ht$unique <- NA;ht$unique <- hrS

lt <- as.data.frame(t(ss[,c(1:(ncol(ps@otu_table))+1)][which(ss$Treatment=="Low feed efficiency"),]))
lrS <- rowSums(lt)
lt$unique <- NA;lt$unique <- lrS


summarise_group_unique <- as.data.frame(cbind(row.names(ht),ht$unique,lt$unique))
colnames(summarise_group_unique) <- c("seq_id","High_efficiency","Low_efficiency")
summarise_group_unique[,2] <- as.numeric(as.character(summarise_group_unique[,2])) 
summarise_group_unique[,3] <- as.numeric(as.character(summarise_group_unique[,3])) 

summarise_group_unique$delta <- NA; summarise_group_unique$delta <- summarise_group_unique[,2] - summarise_group_unique[,3]
summarise_group_unique$taxon <- NA; summarise_group_unique$taxon <- merge(x = summarise_group_unique,y = asv5,by="seq_id",all=T)$phylum
summarise_group_unique <-summarise_group_unique[!grepl("Unclassified",summarise_group_unique$taxon),]
summarise_group_unique <- summarise_group_unique[order(summarise_group_unique$delta,decreasing = T),]


summarise_group_unique <-   summarise_group_unique[which(
summarise_group_unique$delta < quantile(summarise_group_unique$delta)[2]  | 
summarise_group_unique$delta > quantile(summarise_group_unique$delta)[4]  
),]

summarise_group_unique <- summarise_group_unique[complete.cases(summarise_group_unique),]
summarise_group_unique <- summarise_group_unique[order(summarise_group_unique$delta,decreasing = T),]
#delta_values <- ggplot(data=summarise_group_unique,aes(x=factor(reorder(taxon, -delta)),y=delta,fill="black"))+

summarise_group_unique_rel <- summarise_group_unique
summarise_group_unique_rel$High_efficiency <- (summarise_group_unique_rel$High_efficiency /sum(summarise_group_unique_rel$High_efficiency))*100
summarise_group_unique_rel$Low_efficiency <- (summarise_group_unique_rel$Low_efficiency /sum(summarise_group_unique_rel$Low_efficiency))*100
summarise_group_unique_rel$High_efficiency[which(summarise_group_unique_rel$High_efficiency<0.05)]<- NA
summarise_group_unique_rel$Low_efficiency[which(summarise_group_unique_rel$Low_efficiency<0.05)]<- NA
summarise_group_unique_rel <- summarise_group_unique_rel[complete.cases(summarise_group_unique_rel),]

summarise_group_unique_rel_total <- summarise_group_unique
summarise_group_unique_rel_total$total <- NA; summarise_group_unique_rel_total$total <- summarise_group_unique_rel_total$High_efficiency + summarise_group_unique_rel_total$Low_efficiency
summarise_group_unique_rel_total$total_rel <- NA;summarise_group_unique_rel_total$total_rel <- summarise_group_unique_rel_total$total/sum(summarise_group_unique_rel_total$total)*100
summarise_group_unique_rel_total$total_rel[which(summarise_group_unique_rel_total$total_rel<0.05)]<- NA
summarise_group_unique_rel_total <- 
summarise_group_unique_rel_total[complete.cases(summarise_group_unique_rel_total),]

taxon_list <- unique(summarise_group_unique_rel_total$taxon)
taxon_summary <- do.call(rbind,lapply(1:length(taxon_list),function(i){
  summarise_group_unique_rel_total[which(summarise_group_unique_rel_total$taxon==taxon_list[[i]]),]
 x <- as.data.frame(t(colSums(summarise_group_unique_rel_total[which(summarise_group_unique_rel_total$taxon==taxon_list[[i]]),-c(1,5)])))
 x$taxon <- NA; x$taxon <-as.character(taxon_list[[i]]);x <- x[,c(6,1:5)]
 return(x)
}))
#t(cbind(as.character(taxon_list[[1]]),
 
delta_values <- ggplot(data=taxon_summary,aes(x=factor(reorder(taxon, -delta)),y=delta,fill="black"))+
  
 geom_bar(stat="identity",position="dodge",show.legend = F)+
  xlab("")+
  ylab("Abundance change (High /Low Feed efficiency)")+
  # ylim(0,0.1)+
  geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,angle=90,colour="black"),
        axis.text.y  = element_text(size=60,angle=0,colour="black"))
 # scale_y_continuous(limits = c(0,0.8),breaks=seq(0,0.8,0.05))
#delta_values
ggsave(paste0(graph_dir,"/","delta_values_graphics",".pdf"),delta_values,dpi=200,width =400,height=90,units = "cm",scale=1.2,limitsize = FALSE)

##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
mds <- cmdscale(ps_otu_dist, eig = TRUE, k = 2) 
dpc_dist <- dudi.pco(ps_otu_dist, nf = 3, scannf = FALSE)
scatter(dpc_dist)
for_plot <- data.frame(dpc_dist$tab, group = gsub("\\d$", "", ss$Treatment),phase=ss$Phase)

for_plot_pco <- data.frame(dpc_dist$co[,1:2], group = gsub("\\d$", "", ss$Treatment),phase=ss$Phase)

ggplot(for_plot)+
  geom_point(aes(x = X1,y = X2, color = group,shape=phase))


ps_otu_hclust <- hclust(ps_otu_dist,"ward.D")
ps_otu_hclust$labels


plot(ps_otu_hclust)

pp <- prevalence(ps, detection = 0, sort = TRUE)

# dietswap is a phyloseq object; see above
ps.compositional <- microbiome::transform(ps, "compositional")
# Taxa with over 50% prevance at .2% relative abundance
ps.core <- core(ps.compositional, 
                      detection = .2/100, prevalence = 50/100)


# Define data sets to cross-correlate

x <- sqrt(ps@otu_table)   # OTU Log10 (44 samples x 130 genera)
y <- as.matrix(ps@sam_data[,-c(1:5)]) # Lipids (44 samples x 389 lipids)

# Cross correlate data sets and return a table
correlation.table <- associate(x, y, 
                               method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
colnames(correlation.table) [1:2] <- c("seqid","variable")
seqid<-cbind(as.character(asv5$seq_id),as.character(asv5$final_taxon))
colnames(seqid) <- c("seqid","final_taxon")
correlation.table2 <-merge(correlation.table,seqid,by= "seqid",all=T)
correlation.table2 <-correlation.table2[which(correlation.table2$p.adj<=0.05),]
write.csv(correlation.table2,paste0(out_dir,"/","csv/","taxa_cor.csv"),quote = F,row.names = F)

he_cor <- heat(correlation.table2, "final_taxon", "variable", fill = "Correlation", 
    star = "p.adj", p.adj.threshold = 0.05,colours = c("red","white","blue")) 

ggsave(paste0(graph_dir,"/","abundance_correlation",".pdf"),he_cor,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)

##top sp
N<- 30
top_sp <- as.data.frame(sort(taxa_sums(ps)))
top_sp$seq <- NA; top_sp$seq <- row.names(top_sp);colnames(top_sp) <- c("Freq","seqid")
top_sp$Freq <- top_sp$Freq/sum(top_sp$Freq,na.rm = T)*100
top_sp <- merge(top_sp,seqid,by= "seqid",all=T)
top_sp <- top_sp[1:N,]
top_sp$final_taxon <- factor(top_sp$final_taxon)

abund_values <- ggplot(data=top_sp,aes(x=reorder(final_taxon, -Freq),y=Freq,fill="black",colour="black"))+
  geom_bar(stat="identity",show.legend = F,position = "dodge",colour="black",fill="black")+
  xlab("")+
  ylab("Relative abundance (%)")+
  # ylim(0,0.1)+
 # geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,angle = 90,colour="black"),
        axis.text.y  = element_text(size=60,colour="black"))+
scale_y_continuous(limits = c(0,ceiling(max(top_sp$Freq))),breaks=seq(0,ceiling(max(top_sp$Freq)),1))
#abund_values
ggsave(paste0(graph_dir,"/","top30_abund",".pdf"),abund_values,dpi=300,width =80,height=90,units = "cm",scale=1.2,limitsize = FALSE)




g <- global(ps)
g$Reference.number <- NA; g$Reference.number <- row.names(g)
g_alpha <- merge(g,prod,by="Reference.number",all = T) 

mdf <- melt(dat=g_alpha[, c("Sample.ID","Treatment","diversities_shannon")], #,"dominance_simpson"
             id.vars = c("Sample.ID","Treatment"))

###################




alpha_bp<- ggplot(data=mdf,aes(x=variable,y=,value,fill=Treatment))+
  geom_boxplot(position = "dodge",outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Scaled values")+ xlab("")+
  ggtitle("")+
  scale_fill_manual(labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))+
  
 # scale_fill_manual("legend", values = c("red","blue"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black")#,
        #legend.title=element_blank(),legend.position="none")
)
ggsave(paste0(graph_dir,"/","alpha_bp",".pdf"),alpha_bp,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)



###################




 
 
p_test <- plot_richness(ps, "Treatment", "Phase")

prich <- plot_richness(ps, x="FCR", measures=c("Observed","Shannon", "Simpson"), color="Treatment",shape = "Phase")
prich <- prich + scale_colour_manual(labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))+ 
  geom_point(size = 10) + 
  
theme(panel.background = element_rect(fill = "gray95"),
      text=element_text(size=90),
      axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
      axis.text.y  = element_text(size=60,colour="black"),
      legend.key = element_rect(size =10),
      legend.key.size = unit(4, "cm")
)#,

ggsave(paste0(graph_dir,"/","alpha_diversity",".pdf"),prich,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)


res <- plot_frequencies(sample_data(ps),"Treatment","Phase")


pseq2 <- aggregate_taxa(ps, "Phylum", top = 5) 

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.pcoA <- ordinate(ps.prop, method="PCoA", distance="bray")

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
ord.RDA <- ordinate(ps.prop, method="RDA")
# ord.RDA$Ybar <- as.data.frame(ord.RDA$Ybar)
# ord.RDA$colsum <- as.data.frame(ord.RDA$colsum)
# colnames(ord.RDA$Ybar) <- 
# row.names(ord.RDA$colsum) <- as.character(asv5$seq_id)
# ord.RDA$CA
autoplot(ord.RDA,colour="Treatment")





PCoA_plot <- plot_ordination(
  physeq = ps.prop,
  ordination = ord.nmds.pcoA,
  color = "Treatment",
  shape = "Phase",
  title = "Feed efficiency PCoA"
) + 
  scale_colour_manual(labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))+
  geom_point(aes(color = Treatment), alpha = 0.7, size = 8) +
  geom_point(colour = "grey90", size = 1.5) 

PCoA_plot <- PCoA_plot +  theme(panel.background = element_rect(fill = "gray95"),
                  text=element_text(size=90),
                  axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
                  axis.text.y  = element_text(size=60,colour="black"),
                  legend.key = element_rect(size =10),
                  legend.key.size = unit(4, "cm")
)#,

PCoA_plot <- PCoA_plot +
  geom_text_repel(aes(label=sampleID),show.legend = F,size=10)
                                
ggsave(paste0(graph_dir,"/","PCoA",".pdf"),PCoA_plot,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)


p2 <- plot_ordination(physeq = ps.prop, 
                      ordination = ord.nmds.bray,
                      color="Family",
                      title="Bray NMDS",
                      shape = "Treatment",
                      label = "sampleID",
                      type = "split")
p2

po_CCA <- plot_ordination(physeq = ps.prop,ordination = ord.RDA,
                          type="sites", 
                          color="Treatment",title="Bray NMDS",shape = "Phase")

po_CCA
po <- plot_ordination(physeq = ps.prop,ordination = ord.nmds.bray,
                      type="samples", label="sampleID",
                      color="Treatment",title="Bray NMDS",shape = "Phase")




po <- po + scale_colour_manual(labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))
po <- po + geom_point(size=15,show.legend = F)
po <- po +   scale_shape_manual(values=c(16, 17))
po <- po + scale_size_manual(values=c(90,90))
po <- po +  theme(panel.background = element_rect(fill = "gray95"),
                  text=element_text(size=90),
                  axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
                  axis.text.y  = element_text(size=60,colour="black"),
                  legend.key = element_rect(size =10),
                  legend.key.size = unit(4, "cm")
                  )#,



ggsave(paste0(graph_dir,"/","NMDS_Bray_Curtis",".pdf"),po,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)




#PERMANOVA

set.seed(1)

# Calculate bray curtis distance matrix
ps_bray <- phyloseq::distance(ps.prop, method = "bray")
hc_bray <- hclust(ps_bray,"ward.D")
plot(hc_bray)

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps),colour=prod$Treatment)

# Adonis test
ado_comp <- adonis(ps_bray ~ Treatment, data = sampledf,permutations = 10000)
#R2 variance explained by groups
#http://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html

beta <- betadisper(ps_bray, sampledf$Treatment)
permutest(beta,permutations = 10000)

phoc <- with(sampledf,beta)
TukeyHSD(phoc)


# CAP ordinate
cap_ord <- ordinate(
  physeq = ps.prop, 
  method = "RDA",
  distance = ps_bray,
  formula = ~ pH + ADG + FCR + Average_daily_intake + Total_weight_gain
)


scor = scores(cap_ord, display=c("sp", "cn", "bp"), scaling=2) 
# species


species_centroids <- data.frame(scor$species)
species_centroids$species_names <- rownames(species_centroids) # you can also put label = rownames() etc in the ggplot code. But sometimes you may want it in a column for later



# CAP plot
cap_plot <- plot_ordination(
  physeq = ps.prop, 
  ordination = cap_ord, 
  color = "Treatment", 
  axes = c(1,2)
) + 
  aes(shape = Phase) + 
  geom_point(aes(colour = Treatment), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) + 
  scale_colour_manual(labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))



# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = RDA1, 
                 yend = RDA2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * RDA1, 
                 y = 1.3 * RDA2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot <- cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text_repel(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )
cap_plot2 <- cap_plot
cap_plot <- cap_plot2


teplot <- ggplot(data=species_centroids, aes(label = species_names,x = RDA1,y=RDA2,colour = "orange"))+
  geom_text_repel()







ps_tree2 <- plot_tree(ps, nodelabf=nodeplotboot(80,0,3),
                      color="Treatment", 
                      label.tips="taxa_names", 
                      ladderize="left",
                      shape="Phase")#+
# coord_polar(theta="y")

ps_tree2

plot_tree(ps, color="Treatment", size="abundance", plot.margin=0.4)



# data(atlas1006) 
# atlas1006
# data(dietswap) # Data from http://dx.doi.org/10.1038/ncomms7342
# dietswap
dietswap.compositional <- transform(dietswap, "compositional")
g <- global(atlas1006, index = "gini")
sample_data(atlas1006)$diversity <- global(atlas1006, index = "shannon")[,1]

# Compare age and microbiome diversity
plot_regression(diversity ~ age, meta(atlas1006))
p <- prevalence(dietswap, detection = 0, sort = TRUE)

# Taxa with over 50% prevance at .2% relative abundance
dietswap.core <- core(dietswap.compositional, 
                      detection = .2/100, prevalence = 50/100)
dietswap.core

library(ggplot2, quiet = TRUE)
p <- plot_core(transform(dietswap.core, "compositional"), 
               plot.type = "heatmap", 
               colours = gray(seq(0,1,length=5)),
               prevalences = seq(.05, 1, .05), 
               detections = 10^seq(log10(1e-3), log10(.2), length = 10), 
               horizontal = TRUE) +
  xlab("Detection Threshold (Relative Abundance (%))") 
print(p)  





plot_bar(ps, "Order", fill="Family", facet_grid=Treatment~Phase)
