require("DESeq2")
require("KEGGREST")
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
CowPI_dir <- paste0(dat_dir,"/","CowPI"); if(!file.exists(CowPI_dir)){dir.create(CowPI_dir)} # Copy here the production data

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

CowPI_file2a <- readRDS(paste0(csv_dir,"/","CowPI_file2a",".RDS"))

colData1 = CowPI_file2a[,1:10]#data.frame(sample= CowPI_file2a$Subject,Treatment = as.character(CowPI_file2a$Treatment),row.names=CowPI_file2a$Subject); #colnames(colData1) <-"Treatment"# c("Samples","Treatment")
rownames(colData1) <- CowPI_file2a$Subject
# colData1$Treatment <- sub("High feed efficiency","High_feed_efficiency",  colData1$Treatment)
# colData1$Treatment <- sub("Low feed efficiency","Low_feed_efficiency",  colData1$Treatment)
#colData1$Treatment <- factor(colData1$Treatment,levels =c("Low_feed_efficiency","High_feed_efficiency"))
colData1$Treatment <- factor(colData1$Treatment,levels =c("Low feed efficiency","High feed efficiency"))

#Treatment <- colData1$Treatment
colData1$Phase <- factor(CowPI_file2a$Phase,levels =c("liquid","solid"))
#Treatment <- factor(colData$Treatment,levels =c("High feed efficiency","Low feed efficiency"))
#Treatment <- matrix(Treatment)

countData1 = data.frame(t(CowPI_file2a[,-c(1:10)]))
#row.names(countData1) <- gsub(".","_",as.character(row.names(countData1)),fixed = T)

# 
# for(i in 1:ncol(countData1)){
#   cat(class(countData1[,i]),"\n")
# }
#countData = countData[,order(row.names(CowPI_file2))]
colnames(countData1) <- CowPI_file2a$Subject
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsubread")


# all(rownames(colData1) %in% colnames(countData1))
# all(rownames(colData1) == colnames(countData1))

DES_CowPI <- DESeq2::DESeqDataSetFromMatrix(
  countData = countData1,
  colData = colData1,#DataFrame(Treatment,row.names = colnames(countData)),#,
  design =  ~ Treatment 
)

#Applying geometric means
geoMeans = apply(counts(DES_CowPI), 1, gm_mean)
#Estimating size factors

#abundds2 <- getVarianceStabilizedData(DES_CowPI)
diagdds = estimateSizeFactors(DES_CowPI, geoMeans = geoMeans)
diagdds <- estimateDispersions(diagdds)
diagdds2 <- DESeq2::DESeq(diagdds,fitType="local") 

#dds <- DESeq2::DESeq(DES_CowPI)
res <- results(diagdds2)
rld <- rlog(diagdds2)
# res_rld <- results(diagdds2)
sampleDists <- dist( t( assay(rld) ) )
library("pheatmap")
library("RColorBrewer")

#rld
require(FactoMineR)
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pca  <- prcomp(t(assay(rld)[select,]))
#PCA_DES_RLD <- PCA(t(assay(rld)[select,]),ncp = 5)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )


rld_t_DF <-  CowPI_file2a[,-c(1:10)]#as.data.frame(t(assay(rld)[select,]))#(as.data.frame(t(assay(rld))))
rld_t_DF$Treatment <- colData1$Treatment
rld_t_DF$Phase <- colData1$Phase
rld_t_DF2 <- rld_t_DF

#
rld_t_DF2 <- rld_t_DF2[,-c(253,254)]
rld_t_DF2 <- rld_t_DF2[,-c(252)]
#rld_t_DF2 <- rld_t_DF2/colSums(rld_t_DF2)
for(i in 1:nrow(rld_t_DF2)){
  cat(class(rld_t_DF2[i,]),"\n")
  rld_t_DF2[i,] <- rld_t_DF2[i,]/sum(rld_t_DF2[i,])
};rm(i)

LE_rld <- rld_t_DF2[which(rld_t_DF$Treatment=="Low feed efficiency"),]
LE_rld_L <- rld_t_DF2[which(rld_t_DF$Treatment=="Low feed efficiency" | rld_t_DF$Phase=="liquid"),]
LE_rld_S <- rld_t_DF2[which(rld_t_DF$Treatment=="Low feed efficiency" | rld_t_DF$Phase=="solid"),]
#
HE_rld <- rld_t_DF2[which(rld_t_DF$Treatment=="High feed efficiency"),]
HE_rld_L <- rld_t_DF2[which(rld_t_DF$Treatment=="High feed efficiency" | rld_t_DF$Phase=="liquid"),]
HE_rld_S <- rld_t_DF2[which(rld_t_DF$Treatment=="High feed efficiency" | rld_t_DF$Phase=="solid"),]
z_score_df <- data.frame(
row.names =  names(colMeans(LE_rld,na.rm = T)),
Low_feed = as.numeric(colMeans(LE_rld,na.rm = T)-colMeans(rld_t_DF2,na.rm = T))/as.numeric(apply(rld_t_DF2, 2, sd)),
High_feed = as.numeric(colMeans(HE_rld,na.rm = T)-colMeans(rld_t_DF2,na.rm = T))/as.numeric(apply(rld_t_DF2, 2, sd))


)
z_score_df$Pathways <- row.names(z_score_df)
require(ggplot2);require(reshape2)
nba.m <- melt(z_score_df,"Pathways")
 nba.m$variable <- sub("High_feed","High feed efficiency",nba.m$variable)
 nba.m$variable <- sub("Low_feed","Low feed efficiency",nba.m$variable)
nba.m$variable <- factor(nba.m$variable,levels = c("High feed efficiency","Low feed efficiency"))
s1bp_Z_SCORE <- ggplot(data=nba.m,aes(x=variable,y=value,fill=variable))+
  geom_boxplot(position = "dodge",outlier.colour = NA,
               outlier.fill=NA,outlier.alpha=1,
               outlier.size =NA,na.rm = T,show.legend = F) +
  stat_boxplot(geom ='errorbar') +
  scale_y_continuous(limits = c(-0.4,0.4),breaks=seq(-0.4,0.4,0.05))+
  geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
    #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Z score")+ xlab("")+
  ggtitle("")+
  scale_fill_manual("Groups", values = c("blue","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black")#,
        #legend.title=element_blank(),legend.position="none"
  )+
  coord_flip()
s1bp_Z_SCORE
ggsave(paste0(graph_dir,"/","s1bp_Z_SCORE",".pdf"),s1bp_Z_SCORE,dpi=300,width =100,height=80,units = "cm",scale=1.2,limitsize = FALSE)


# 
# p_hmap <- ggplot(nba.m, aes(variable, Pathways)) + 
#   geom_tile(aes(fill = value),colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")
# ggsave(paste0(graph_dir,"/","p_hmap",".pdf"),p_hmap,dpi=300,width =100,height=200,units = "cm",scale=1.2,limitsize = FALSE)
# 
# p_hmap
################
#PCA APPROAH

#PCA_DES <- plotPCA(rld, intgroup = c("Treatment"),returnData=T)
#3PCA_DES$PC3 <- pca$x[,3]
#PCA_DES <- PCA_DES[,c(1,2,6,3,4,5)]


#####nipals
data_rlog <- as.data.frame(t(assay(rld)))
write.csv(data_rlog,"E:/rlog.csv",row.names = T,quote = F)
#data <- data[,complete.cases(data)]
for(i in 1:ncol(data_rlog)){
 #cat(nrow(data)-length(na.omit(data[,i]))," | ",i,"\n")
cat(sum(data_rlog[,i],na.rm = T)," | ",sum(is.na(data_rlog[,i]))," | ",class(data_rlog[,i])," | ",i,"\n")
  
}
data_rlog <- data_rlog[,which(colSums(data_rlog)>0)]
library(usdm)
library(plsdepot)
z = plsdepot::nipals(Data=data_rlog, comps=5, scaled=T)

require(rgl);require(car)
#scatter3d(x = PCA_DES$PC1, y = PCA_DES$PC2, z = PCA_DES$PC3, 
scatter3d(x = z$scores[,1], y = z$scores[,2], z = z$scores[,3], 
          groups =  factor(paste0(colData1$Treatment,"-",colData1$Phase),levels = 
  c("High feed efficiency-solid","High feed efficiency-liquid","Low feed efficiency-solid", "Low feed efficiency-liquid")),
          # xlab= paste("Axis 1",round(percentVar[1]*100,2),"%"),
          # ylab= paste("Axis 2",round(percentVar[2]*100,2),"%"),
          # zlab= paste("Axis 3",round(percentVar[3]*100,2),"%"),
  xlab= paste("Axis 1",round(z$values$percentage[1],2),"%"),
  ylab= paste("Axis 2",round(z$values$percentage[2],2),"%"),
  zlab= paste("Axis 3",round(z$values$percentage[3],2),"%"),
  
            surface=FALSE, grid = FALSE, ellipsoid = TRUE,
          labels = row.names(PCA_DES),
          point.col =  c("black", "black", "black"),id.n=length(row.names(PCA_DES)),
          surface.col = c("darkblue", "green", "orange","red"),
          axis.col = c("black", "black", "black")
)
require(car);require(rgl);require(ggrepel)
legend3d("topright", legend = c("High feed efficiency - solid","High feed efficiency - liquid",
                                "Low feed efficiency - solid","Low feed efficiency - liquid"),
         pch = 16, col = c("darkblue", "green", "orange","red"), cex=2, inset=c(0.01))
require(ggrepel)

PCA_DES2 <- plotPCA(rld, intgroup = c("Treatment"),returnData=F)
PCA_DES2 <- PCA_DES2 +
  geom_point(size=5,aes(shape=colData1$Phase,color=colData1$Treatment))+
  scale_color_manual(values = c("red","blue"))+
  scale_shape_manual("Phase",values=c(16,17))+
  geom_text_repel(aes(label=colData1$sampleID), size=5,colour="black")
PCA_DES2
###########################################################
hclust_CowPi <- hclust(sampleDists,"ward.D")
plot(hclust_CowPi)
# hclust_CowPi$order
# merge(hclust_CowPi$order,CowPI_file2,by.x="row.names",by.y="Subject")
###CLuster
require(ggdendro)
#We will color the labels according to countries(group_info[,1])
hc_d <- dendro_data(as.dendrogram(hclust_CowPi))
hc_d$labels$Type <- 
  merge(hc_d$labels,colData1,by.y="Subject",by.x="label")[,"Phase"]
hc_d$labels$Treatment <-  merge(hc_d$labels,colData1,by.y="Subject",by.x="label")[,"Treatment"]
hc_d$labels$SubjectID <-  merge(hc_d$labels,colData1,by.y="Subject",by.x="label")[,"sampleID"]
hc_d$labels$SubjectID <-  paste0(hc_d$labels$SubjectID," - ",hc_d$labels$Type)
#  merge(hc_d$labels,labels,by.y="row.names",by.x="label")[,5]

hc_d$labels
hc_d$labels$color <- NA
hc_d$labels$color[which(hc_d$labels$Treatment=="Low feed efficiency")] <- "red"
hc_d$labels$color[which(hc_d$labels$Treatment=="High feed efficiency")] <- "blue"
#hc_d$labels$color=cols[hc_d$labels$Type]

## Plot clusters
p1 <- ggplot(data = segment(hc_d)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  scale_x_discrete(labels=label(hc_d)$Type) +
  geom_point(data=hc_d$label, aes(x = x, y = y, color = color), alpha = 1) +
  #p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
  #scale_color_manual(values = cols)+
  scale_color_manual(name="Groups",labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))+
  ylab("Distance (Metabolic pathways: Euclidean)") + theme_bw()+
  theme(axis.text.y = element_text(color = hc_d$labels$Treatment),
        axis.title.y = element_blank())
p1 <- p1 + geom_text(data=label(hc_d),
                     aes(label=SubjectID, x=x, y=-1, colour=hc_d$labels$color))
p1

library(vegan)
permanova <- adonis(sampleDists ~
                      colData1$Treatment,
                    #otu_table(ps3) ~ sample_data(ps3)$Treatment,
                    permutations=10000, method = "bray")
beta <- betadisper(sampleDists, colData1$Treatment)
perm_beta <- permutest(beta,permutations = 10000)
plot(beta)

################33


sampleDistMatrix <- as.matrix( sampleDists )
rS <- data.frame(rownames(sampleDistMatrix));row.names(rS) <-rownames(sampleDistMatrix)
p_label <- paste0(merge(rS,rld@colData,by.x="row.names",by.y="Subject")$sampleID,"-",
                  merge(rS,rld@colData,by.x="row.names",by.y="Subject")$Phase)           
rownames(sampleDistMatrix) <- p_label
# #colnames(sampleDistMatrix) <- p_label
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# df <- as.data.frame(colData(diagdds2)[,c("Treatment","Phase")])
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          show_colnames=FALSE,
#          annotation_colors = list(Treatment= c(Low_feed_efficiency="red",High_feed_efficiency= "blue")),
#          annotation_col=df,cluster_cols=FALSE)
require(jmvcore)

resOrdered <- res[order(res$pvalue),]
res05 <- resOrdered[which(resOrdered$padj <= 0.05),]
require(ggplot2);require(stringr)
res05$metabolic_pathway <- row.names(res05)
res05 <- res05[order(res05$log2FoldChange,decreasing = T),]
res05$metabolic_pathway <-gsub("."," ",res05$metabolic_pathway,fixed = T)
res05$KEGG_ID <- NA
res05$Class_raw <- NA
res05$Class <- NA
for(i in 1:nrow(res05)){
  cat(i,"\n")
  query <- names(keggFind("pathway",query = str_trim(res05$metabolic_pathway[[i]])))
  if(is.null(query)){
    res05$Class_raw[[i]] <- NA
    res05$KEGG_ID[[i]] <- NA
  } else {
  res05$Class_raw[[i]] <- keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS
if(length(query)>1){
  res05$KEGG_ID[[i]] <- paste(query[[1]],"*")
} else {
    res05$KEGG_ID[[i]] <- query
    }
  }
};rm(i)

res05$ID <- 1:nrow(res05)

#View(cbind(res05$metabolic_pathway,res05$Class_raw))
res05$Class_raw[[5]] <- paste0(res05$metabolic_pathway[[5]],"*")
res05$Class_raw[[9]] <- "Genetic Information Processing; Folding"
res05$Class_raw[[17]] <- "Genetic Information Processing; Folding"
res05$Class_raw[[18]] <- paste0(res05$metabolic_pathway[[18]],"*")
res05$Class_raw[[19]] <- "Genetic Information Processing; Folding"
res05$Class_raw[[21]] <- paste0(res05$metabolic_pathway[[21]],"*")
res05$Class_raw[[22]] <- paste0(res05$metabolic_pathway[[22]],"*")
res05$Class_raw[[32]] <- paste0(res05$metabolic_pathway[[32]],"*")
res05$Class_raw[[33]] <- paste0(res05$metabolic_pathway[[33]],"*")
res05$Class_raw[[38]] <- paste0(res05$metabolic_pathway[[38]],"*")
res05$Class_raw[[39]] <- paste0(res05$metabolic_pathway[[39]],"*")

class_raw <- str_split(res05$Class_raw,";")
class_raw<- lapply(1:nrow(res05),function(i){
  x <- class_raw[[i]][1]
return(x)
  })
res05$Class <- unlist(class_raw)
#res05_a <- res05[which(res05$log2FoldChange > 0.1 | res05$log2FoldChange < -0.1 ),]
res05_a <- res05[which(res05$Class !="Human Diseases"),]
write.csv(res05,paste0(csv_dir,"/","res05",".csv"))

table(res05_a$Class)

      
table(res05_a$Class[which(res05_a$log2FoldChange<0)])[order(table(res05_a$Class[which(res05_a$log2FoldChange<0)]),decreasing = T)]
table(res05_a$Class[which(res05_a$log2FoldChange>0)])[order(table(res05_a$Class[which(res05_a$log2FoldChange>0)]),decreasing = T)]


exp_plot_LF <- ggplot(as.data.frame(res05_a), aes(x=reorder(metabolic_pathway, -log2FoldChange) ,y=log2FoldChange,fill=Class,colour=Class)) + 
  #geom_point(size=20) + 
  geom_bar(stat="identity")+
  geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
  ylab("Log2foldchange")+
  xlab("")+
  #geom_line()+
  theme(text=element_text(size=69),
        legend.text=element_text(face="italic",size=60,colour="black"),
        axis.text.x  = element_text(face="italic",angle = -90,size=60, hjust = 0, vjust=0.5),
        axis.text.y  = element_text(face="italic",size=60,colour="black"))+
  scale_y_continuous(limits = c(-2.2,1.1),breaks=seq(-2.2,1.1,0.2))+
  coord_flip()

exp_plot_LF

ggsave(paste0(graph_dir,"/","exp_plot_LF",".pdf"),exp_plot_LF,dpi=300,width =180,height=140,units = "cm",scale=1.2,limitsize = FALSE)

 test <- as.data.frame(cbind(as.character(CowPI_file2a$sampleID),
                             as.character(CowPI_file2a$Subject),
        as.character(CowPI_file2a$Treatment),
        as.character(CowPI_file2a$Phase),
        CowPI_file2a$Limonene.and.pinene.degradation

        )
)
 test$V5 <- as.numeric(as.character((test$V5)))
sum_cowpi <- data.frame(colSums(CowPI_file2a[,-c(1:10)]))
sum_cowpi$Pathway <- row.names(sum_cowpi)
colnames(sum_cowpi) <- c("Total","Pathway")
row.names(sum_cowpi) <- 1:nrow(sum_cowpi)
#colnames(sum_cowpi)[1] <- "sum"
#sum_cowpi <- sum_cowpi[order(sum_cowpi$sum,decreasing = T),]
#LE <- CowPI_file2a[which(CowPI_file2a$Treatment=="Low feed efficiency"),-c(1:10)]
sum_cowpi$LE <- colSums(CowPI_file2a[,-c(1:10)][which(CowPI_file2a$Treatment=="Low feed efficiency"),])
sum_cowpi$HE <- colSums(CowPI_file2a[,-c(1:10)][which(CowPI_file2a$Treatment=="High feed efficiency"),])
sum_cowpi$MEAN <- colMeans(CowPI_file2a[,-c(1:10)])
sum_cowpi$M_LE <- colMeans(CowPI_file2a[,-c(1:10)][which(CowPI_file2a$Treatment=="Low feed efficiency"),])
sum_cowpi$M_HE <- colMeans(CowPI_file2a[,-c(1:10)][which(CowPI_file2a$Treatment=="High feed efficiency"),])

# sum_cowpi$LE <- merge(sum_cowpi,
#                 data.frame(colSums(CowPI_file2a[,-c(1:10)][which(CowPI_file2a$Treatment=="Low feed efficiency"),])),
#                 by.x="Pathway",by.y="row.names")[,3]
# sum_cowpi$HE <- merge(sum_cowpi,
#                       data.frame(colSums(CowPI_file2a[,-c(1:10)][which(CowPI_file2a$Treatment=="High feed efficiency"),])),
#                       by.x="Pathway",by.y="row.names")[,3]

sum_cowpi <-sum_cowpi[,c(2,1,3,4)]


mean(subset(test,test$V3=="High feed efficiency")$V5)
     mean(subset(test,test$V3=="Low feed efficiency")$V5)
        
  log2(mean(subset(test,test$V3=="High feed efficiency")$V5)/
 mean(subset(test,test$V3=="Low feed efficiency")$V5)
 )
#install.packages("apeglm") 
  # require(apeglm)
  # drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
  # plotMA(res)

  cowPi_assoc <- microbiome::associate(as.data.frame(t(assay(rld)[select,])), CowPI_file2a[,c(5:9)], 
                                       method = "spearman",
                                       mode = "table", p.adj.threshold = 0.05, n.signif = 1,
                                       p.adj.method = "fdr")
  cowPi_assoc_plot <- microbiome::heat(cowPi_assoc, "X1", "X2", fill = "Correlation", colours = c("darkred","red","white","blue","darkblue"),
       star = "p.adj", p.adj.threshold = 0.05,plot.values = F,star.size=40,legend.text = "Spearman"
       
  )+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black")) +
  #scale_fill_continuous(guide = guide_colorbar(direction = "horizontal")) +
  theme(legend.position="right",legend.direction = "vertical",legend.key.size =  unit(1.5, "in"))
  
  write.csv(cowPi_assoc,paste0(out_dir,"/","csv/","correlation_met_cowpi_deseq2.csv"),quote = F,row.names = F)
  ggsave(paste0(graph_dir,"/","heat_met_abund",".pdf"),cowPi_assoc_plot,dpi=300,width =160,height=130,units = "cm",scale=1.2,limitsize = FALSE)
  
  #http://www.bioconductor.org/packages//2.13/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
  ################################
  library(MASS)
  require(caret)
  require(plyr)
  CowPI_file6 <- as.data.frame(t(assay(rld)[select,]))
  CowPI_file6$Subject <- row.names(CowPI_file6)
  CowPI_file2a$Subject <- as.character(CowPI_file2a$Subject)
  
  
  CowPI_file6$FCR <-  join(CowPI_file6,CowPI_file2a,"Subject")$FCR
  
  
  
  grp = data.frame(Subject =join(CowPI_file6,CowPI_file2a,"Subject")$Subject,
    grp =factor(paste0(join(CowPI_file6,CowPI_file2a,"Subject")$Treatment," - ",
                      join(CowPI_file6,CowPI_file2a,"Subject")$Phase)
               ,levels = 
                 c("High feed efficiency - solid","High feed efficiency - liquid",
                   "Low feed efficiency - solid", "Low feed efficiency - liquid"))
  )
  grp$Subject <- as.character(grp$Subject)
  # CowPI_file6$grp <- grp
  #CowPI_file2a$Subject  <- NULL
  

  cowPI <- list(production = CowPI_file2a[,c(1:10)],
  grp = grp,
  CowPI_file6 = CowPI_file6)
  
  saveRDS(cowPI,paste0(out_dir,"/","csv/","cowPI_R_FCR.RDS"))
  
  
  
  if(!file.exists(paste0(out_dir,"/","csv/","cowPI_step_list_des.RDS"))){
    # Set up repeated k-fold cross-validation
    set.seed(300)
    
    train.control <- trainControl(method = "cv", number = 10)
    # Train the model
    step.model <- train(FCR ~., data = CowPI_file6,
                        method = "pls", #="leapseq"
                        preproc= c("zv","center","scale"),
                        trControl=train.control,
                        #metric="Accuracy",
                        # trControl=(trainControl(method="repeatedcv", repeats=3,
                        # classProbs=TRUE, summaryFunction=twoClassSummary)),
                        tuneLength=20
                      # metric="ROC"
                       #tuneGrid = data.frame(nvmax = 1:10),
#                        trControl = train.control
    )
    plot(step.model$finalModel, ncomp=15, asp=1)
    
    
    plot(step.model)

            step.model$results
    step.model$bestTune
    summary(step.model$finalModel)
    
    coef(step.model$finalModel,15)
    plot(step.model$finalModel)
    summary(step.model$finalModel)
    
    pred <- predict(step.model, newdata = CowPI_file6)
    #acc
    pred_df <- data.frame(Subject=row.names(CowPI_file6) ,pred = pred, FCR = CowPI_file6$FCR)
    pred_df$Subject <- as.character(pred_df$Subject)
    #row.names(pred_df) <- row.names(CowPI_file6)
    pred_df$Treatment <- NA; pred_df$Treatment <- join(pred_df,grp,"Subject")$grp
    
    library(ggplot2)
    
    # ggplot(data = pred_df)+
    #   geom_line(aes(x = FCR, y = pred_df))+
    #   geom_point(aes(x = FCR, y = pred_df))
    cor.test(pred_df$pred,pred_df$FCR)
    RMSE(pred = pred_df$pred,obs = pred_df$FCR)
    
    coef(step.model$finalModel,15)
    
  
    require(ggpubr)
  
    importance <- varImp(step.model)$importance#arRF$
    importance$PA <- row.names(importance)
    importance <- importance[,c(2,1)]
    importance <- importance[order(importance$Overall,decreasing = TRUE),]
    importance <- importance[1:25,]
    
   
plot(varImp(step.model))
    ggscatterhist_plot <-ggscatterhist(
      pred_df, x = "FCR", y = "pred",group="Treatment",
      color = "Treatment", size = 3, alpha = 0.6,
      palette = c("darkblue", "green", "orange","red"),
      #margin.params = list(fill = "group", color = "black", size = 0.2,
      margin.plot = "boxplot", centrality.para = "mean",
      marginal=T)
    ggscatterhist_plot
    cowPI_step <- list(step.model= step.model,
                       pred = pred,
                       ggscatterhist_plot= ggscatterhist_plot,
                       varImp = varImp(step.model),
                       importance=importance)
    saveRDS(cowPI_step,paste0(out_dir,"/","csv/","cowPI_step_DESEQ_list.RDS"))
    #+
    #  geom_point(data = CowPI_file6, aes(x=FCR, y = dist))
    
    
  } else {
    
    cowPI_step <- readRDS(paste0(out_dir,"/","csv/","cowPI_step_list.RDS"))
  }
  
 