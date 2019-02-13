##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##COWPI OUTPUT ANALYSIS USING DESeq2 approach    #
##################################################

require("DESeq2")
require("phyloseq")
require("KEGGREST")
require("plyr")
require("BBmisc")
require("ggplot2")
##################################################

#calling colnames

findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == thecolumnname)
}

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
ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))
ps3 <-  readRDS(paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))
ps3_Dist <- phyloseq::distance(ps3, method = "jaccard",binary = TRUE)

LE_Ids <- ps@sam_data[which(ps@sam_data$Treatment== "Low feed efficiency")]$Subject
HE_Ids <- ps@sam_data[which(ps@sam_data$Treatment== "High feed efficiency")]$Subject
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
# pca  <- prcomp(t(assay(rld)[select,]))
# #PCA_DES_RLD <- PCA(t(assay(rld)[select,]),ncp = 5)
# percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
# 
# 
# rld_t_DF <-  CowPI_file2a[,-c(1:10)]#as.data.frame(t(assay(rld)[select,]))#(as.data.frame(t(assay(rld))))
# rld_t_DF$Treatment <- colData1$Treatment
# rld_t_DF$Phase <- colData1$Phase
# rld_t_DF2 <- rld_t_DF
# 
# #
# rld_t_DF2 <- rld_t_DF2[,-c(243,244)] #-c(253,254)
# rld_t_DF2 <- rld_t_DF2[,-c(242)] #-c(252)
# #rld_t_DF2 <- rld_t_DF2/colSums(rld_t_DF2)
# for(i in 1:nrow(rld_t_DF2)){
#   cat(class(rld_t_DF2[i,]),"\n")
#   rld_t_DF2[i,] <- rld_t_DF2[i,]/sum(rld_t_DF2[i,])
# };rm(i)
# 
# LE_rld <- rld_t_DF2[which(rld_t_DF$Treatment=="Low feed efficiency"),]
# LE_rld_L <- rld_t_DF2[which(rld_t_DF$Treatment=="Low feed efficiency" | rld_t_DF$Phase=="liquid"),]
# LE_rld_S <- rld_t_DF2[which(rld_t_DF$Treatment=="Low feed efficiency" | rld_t_DF$Phase=="solid"),]
# #
# HE_rld <- rld_t_DF2[which(rld_t_DF$Treatment=="High feed efficiency"),]
# HE_rld_L <- rld_t_DF2[which(rld_t_DF$Treatment=="High feed efficiency" | rld_t_DF$Phase=="liquid"),]
# HE_rld_S <- rld_t_DF2[which(rld_t_DF$Treatment=="High feed efficiency" | rld_t_DF$Phase=="solid"),]
# z_score_df <- data.frame(
# row.names =  names(colMeans(LE_rld,na.rm = T)),
# Low_feed = as.numeric(colMeans(LE_rld,na.rm = T)-colMeans(rld_t_DF2,na.rm = T))/as.numeric(apply(rld_t_DF2, 2, sd)),
# High_feed = as.numeric(colMeans(HE_rld,na.rm = T)-colMeans(rld_t_DF2,na.rm = T))/as.numeric(apply(rld_t_DF2, 2, sd))
# 
# 
# )
# z_score_df$Pathways <- row.names(z_score_df)
# require(ggplot2);require(reshape2)
# nba.m <- melt(z_score_df,"Pathways")
#  nba.m$variable <- sub("High_feed","High feed efficiency",nba.m$variable)
#  nba.m$variable <- sub("Low_feed","Low feed efficiency",nba.m$variable)
# nba.m$variable <- factor(nba.m$variable,levels = c("High feed efficiency","Low feed efficiency"))
# s1bp_Z_SCORE <- ggplot(data=nba.m,aes(x=variable,y=value,fill=variable))+
#   geom_boxplot(position = "dodge",outlier.colour = NA,
#                outlier.fill=NA,outlier.alpha=1,
#                outlier.size =NA,na.rm = T,show.legend = F) +
#   stat_boxplot(geom ='errorbar') +
#   scale_y_continuous(limits = c(-0.4,0.4),breaks=seq(-0.4,0.4,0.05))+
#   geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
#     #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
#   #stat_summary(fun.y=mean, geom="point")+
#   #guides(fill=FALSE)+
#   #ylim(c(0,100))+
#   ylab("Z score")+ xlab("")+
#   ggtitle("")+
#   scale_fill_manual("Groups", values = c("blue","red"))+
#   theme(panel.background = element_rect(fill = "gray95"),
#         text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
#         axis.text.y  = element_text(size=60,colour="black")#,
#         #legend.title=element_blank(),legend.position="none"
#   )+
#   coord_flip()
# s1bp_Z_SCORE
# ggsave(paste0(graph_dir,"/","s1bp_Z_SCORE",".pdf"),s1bp_Z_SCORE,dpi=300,width =100,height=80,units = "cm",scale=1.2,limitsize = FALSE)
# 

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
write.csv(data_rlog,"E:/rlog2.csv",row.names = T,quote = F)
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
hclust_CowPi <-  hclust(sampleDists,"ward.D")
plot(hclust_CowPi)
# hclust_CowPi$order
# merge(hclust_CowPi$order,CowPI_file2,by.x="row.names",by.y="Subject")
###CLuster
require(ggdendro)
#We will color the labels according to countries(group_info[,1])
hc_d <- dendro_data(as.dendrogram(hclust_CowPi))
colnames(hc_d$labels)[3] <- "Subject"
hc_d$labels$Type <- dplyr::left_join(hc_d$labels,colData1,"Subject")$Phase
hc_d$labels$Treatment <-  dplyr::left_join(hc_d$labels,colData1,"Subject")$Treatment
hc_d$labels$SubjectID <-  dplyr::left_join(hc_d$labels,colData1,"Subject")$sampleID
hc_d$labels$SubjectID <-  paste0(hc_d$labels$SubjectID," - ",hc_d$labels$Type)



# hc_d$labels$Type <- 
#   merge(hc_d$labels,colData1,by.y="Subject",by.x="label")[,"Phase"]
# hc_d$labels$Treatment <-  merge(hc_d$labels,colData1,by.y="Subject",by.x="label")[,"Treatment"]
# hc_d$labels$SubjectID <-  merge(hc_d$labels,colData1,by.y="Subject",by.x="label")[,"sampleID"]
# hc_d$labels$SubjectID <-  paste0(hc_d$labels$SubjectID," - ",hc_d$labels$Type)
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

sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rS <- data.frame(rownames(sampleDistMatrix));row.names(rS) <-rownames(sampleDistMatrix)

colnames(rS) <- "Subject"
p_label <- paste0(dplyr::left_join(rS,as.data.frame(rld@colData),"Subject")$sampleID,"-",dplyr::left_join(rS,as.data.frame(rld@colData),"Subject")$Phase)

#hc_d$labels$Type <-left_join(hc_d$labels,labels,"V1")[,8]
# 
# 
# p_label <- paste0(merge(rS,rld@colData,by.x="row.names",by.y="Subject")$sampleID,"-",
#                   merge(rS,rld@colData,by.x="row.names",by.y="Subject")$Phase)           
ps3_Dist_m <- as.matrix(ps3_Dist)

if(length(match(rownames(as.matrix(ps3_Dist)),rownames(sampleDistMatrix))) < length(rownames(sampleDistMatrix))){
  cat("check row names","\n")
} else {
  rownames(sampleDistMatrix) <- p_label
  rownames(ps3_Dist_m) <- p_label
}


vegan::mantel(ps3_Dist_m,sampleDistMatrix,permutations = 10000,na.rm = T,parallel = 2,method = "pearson")
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
res05$Description <- NA
res05$Brite  <- NA
res05$Brite_description <- NA
###############################################
for(i in 1:nrow(res05)){
  
  cat(i,"\n")
  
  ##PATHWAYS
  query <- names(keggFind("pathway",query = str_trim(res05$metabolic_pathway[[i]])))
  
  
  if(is.null(query)){
    res05$Class_raw[[i]] <- NA
    res05$KEGG_ID[[i]] <- NA
    res05$Description[[i]] <- NA
    
  } else {
    #
    res05$Class_raw[[i]] <- if(is.null(keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS)){res05$Class_raw[[i]] <- "CHECK" } else { res05$Class_raw[[i]] <- keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS}
    res05$Description[[i]] <- if(is.null(keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$DESCRIPTION)){res05$Description[[i]] <- "CHECK" } else { res05$Description[[i]] <- paste(keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$DESCRIPTION,collapse="|")}
    
    #res05$Class_raw[[i]] <-# keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS
    if(length(query)>1){
      res05$KEGG_ID[[i]] <- paste(query[[1]],"*")
    } else {
      res05$KEGG_ID[[i]] <- query
    }
  }
  ##KO
  query_brite <- names(keggFind("brite",query = str_trim(res05$metabolic_pathway[[i]])))
  #cat(query_brite,"\n")
  if(is.null(query_brite)){
    
    
    res05$Brite[[i]] <- NA
    res05$Brite_description[[i]] <- NA  
    # res05$Brite[[i]] <- query_brite[[1]]
    # res05$Brite_description[[i]] <- paste0(keggGet(names(keggFind("brite",query=query_brite[[1]]))),collapse = "|")
    # 
    
  } else {
  
    if(length(grep("br:br",query_brite))>0){
      res05$Brite[[i]] <- NA
      res05$Brite_description[[i]] <- NA
      
    } else {
      an.error.occured <- F
      
      tryCatch( { work <- keggGet(names(keggFind("brite",query=query_brite[[1]])))}
                , error = function(e) {an.error.occured <<- TRUE})
      
      if(an.error.occured==T | is.null(keggGet(names(keggFind("brite",query=query_brite[[1]])))) ){
        res05$Brite[[i]] <- NA
        res05$Brite_description[[i]] <- NA 
      } else {
        res05$Brite[[i]] <- query_brite[[1]]
        res05$Brite_description[[i]] <- paste0(keggGet(names(keggFind("brite",query=query_brite[[1]]))),collapse = "|")
      } 
    }
  }
};rm(i)
###############################################

res05$ID <- 1:nrow(res05)
#View(res05)
#View(cbind(res05$metabolic_pathway,res05$Class_raw))
 res05$Class_raw[[1]] <- "Organismal Systems; Immune system*"
 res05$Class_raw[[16]] <- "Peptidases*"
 res05$Class_raw[[22]] <- "Metabolism; Metabolism of cofactors and vitamins*"
 res05$Class_raw[[25]] <- "Other transporters*"
 res05$Class_raw[[28]] <- "Biosynthesis and biodegradation of secondary metabolites*"
 res05$Class_raw[[29]] <- "Transcription machinery*"
# res05$Class_raw[[5]] <- paste0(res05$metabolic_pathway[[5]],"*")
# res05$Class_raw[[9]] <- "Genetic Information Processing; Folding"
# res05$Class_raw[[17]] <- "Genetic Information Processing; Folding"
# res05$Class_raw[[18]] <- paste0(res05$metabolic_pathway[[18]],"*")
# res05$Class_raw[[19]] <- "Genetic Information Processing; Folding"
# res05$Class_raw[[21]] <- paste0(res05$metabolic_pathway[[21]],"*")
# res05$Class_raw[[22]] <- paste0(res05$metabolic_pathway[[22]],"*")
# res05$Class_raw[[32]] <- paste0(res05$metabolic_pathway[[32]],"*")
# res05$Class_raw[[33]] <- paste0(res05$metabolic_pathway[[33]],"*")
# res05$Class_raw[[38]] <- paste0(res05$metabolic_pathway[[38]],"*")
# res05$Class_raw[[39]] <- paste0(res05$metabolic_pathway[[39]],"*")

class_raw <- str_split(res05$Class_raw,";")
class_raw<- lapply(1:nrow(res05),function(i){
  x <- class_raw[[i]][1]
return(x)
  })
res05$Class <- unlist(class_raw)
# 
# ##########ADDING COUNTS##############
 res05$count_per_sample=NA
 res05$count_per_HE=NA
 res05$count_per_LE=NA
 res05$count_per_SO=NA
 res05$count_per_LI=NA
####
 res05$COUNT=NA
 res05$HE=NA
 res05$HE_SD=NA
 res05$LE=NA
 res05$LE_SD=NA
 res05$SO=NA
 res05$SO_SD=NA
 res05$LI=NA
 res05$LI_SD=NA
 res05$HE_LE_DIR=NA
 ####
 
 #res05$metabolic_pathway <-gsub("."," ",res05$metabolic_pathway,fixed = T)
for(i in 1:nrow(res05)){

 res05$count_per_sample[[i]] <- sum(CowPI_file2a[,names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                  ]>0)
 
 res05$count_per_HE[[i]] <- sum(CowPI_file2a[CowPI_file2a$Treatment=="High feed efficiency",
             names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                      ]>0)

 res05$count_per_LE[[i]] <- sum(CowPI_file2a[CowPI_file2a$Treatment=="Low feed efficiency",
                  names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                  ]>0)
 
 res05$count_per_SO[[i]] <- sum(CowPI_file2a[CowPI_file2a$Phase=="solid",
                                          names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                          ]>0)
 
 res05$count_per_LI[[i]] <- sum(CowPI_file2a[CowPI_file2a$Phase=="liquid",
                                         names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                         ]>0)
 
 

 
 res05$COUNT[[i]] <- mean(CowPI_file2a[,names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                  ])
 res05$HE[[i]]    <-  res05$COUNT[[i]] <- mean(CowPI_file2a[CowPI_file2a$Treatment=="High feed efficiency",
                                                            names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                                            ])
 res05$HE_SD[[i]]    <-  res05$COUNT[[i]] <- sd(CowPI_file2a[CowPI_file2a$Treatment=="High feed efficiency",
                                                            names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                                            ])
 res05$LE[[i]]    <-  res05$COUNT[[i]] <- mean(CowPI_file2a[CowPI_file2a$Treatment=="Low feed efficiency",
                                                            names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                                            ])
 res05$LE_SD[[i]]    <-  res05$COUNT[[i]] <- sd(CowPI_file2a[CowPI_file2a$Treatment=="Low feed efficiency",
                                                            names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                                            ])
  res05$SO[[i]]    <-  res05$COUNT[[i]] <- mean(CowPI_file2a[CowPI_file2a$Phase=="solid",
                                                            names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                                            ])
  res05$SO_SD[[i]]    <-  res05$COUNT[[i]] <- sd(CowPI_file2a[CowPI_file2a$Phase=="solid",
                                                             names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                                             ])
  
 res05$LI[[i]]    <-  res05$COUNT[[i]] <- mean(CowPI_file2a[CowPI_file2a$Phase=="liquid",
                                                            names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                                            ])
 res05$LI_SD[[i]]    <-  res05$COUNT[[i]] <- sd(CowPI_file2a[CowPI_file2a$Phase=="liquid",
                                                            names(CowPI_file2a)[names(CowPI_file2a) == as.character(row.names(res05)[i])]      
                                                            ])
 
};rm(i)
 
 res05$HE_LE_DIR <-  res05$HE > res05$LE
 res05$HE_LE_DIR[which(res05$HE_LE_DIR==T)] <- "HE"
 res05$HE_LE_DIR[which(res05$HE_LE_DIR==F)] <- "LE"
 res05$SO_LI_DIR <-  res05$SO > res05$LI
 res05$SO_LI_DIR[which(res05$SO_LI_DIR==T)] <- "SO"
 res05$SO_LI_DIR[which(res05$SO_LI_DIR==F)] <- "LI"
 res05$SO_LI_DIR[which(res05$SO_LI_DIR==F)] <- "LI"
 res05$status <- paste(res05$HE_LE_DIR,"-",res05$SO_LI_DIR)
 res05 <- res05[,-13]
  #res05_a <- res05[which(res05$log2FoldChange > 0.1 | res05$log2FoldChange < -0.1 ),]
res05_a <- res05[which(res05$Class !="Human Diseases"),]

#res05
write.csv(res05,paste0(csv_dir,"/","meta",".csv"),quote = T,row.names = F)
write.csv(res05_a,paste0(csv_dir,"/","res05a",".csv"),quote = T,row.names = F)

#   table(res05_a$Class)
# 
#       
# table(res05_a$Class[which(res05_a$log2FoldChange<0)])[order(table(res05_a$Class[which(res05_a$log2FoldChange<0)]),decreasing = T)]
# table(res05_a$Class[which(res05_a$log2FoldChange>0)])[order(table(res05_a$Class[which(res05_a$log2FoldChange>0)]),decreasing = T)]

res05_a$Class[[1]] <- "NO KEGG PATH AVAILABLE"
res05_a$Class[[13]] <- "NO KEGG PATH AVAILABLE"
#res05_a$Class[[16]] <- "NO KEGG PATH AVAILABLE"
res05_a$Class[[22]] <- "NO KEGG PATH AVAILABLE"
res05_a$Class[[25]] <- "NO KEGG PATH AVAILABLE"
res05_a$Class[[26]] <- "NO KEGG PATH AVAILABLE"
#res05_a$Class[[28]] <- "NO KEGG PATH AVAILABLE"
res05_a$Class[[29]] <- "NO KEGG PATH AVAILABLE"

exp_plot_LF <- ggplot(as.data.frame(res05_a), aes(x=reorder(metabolic_pathway, -log2FoldChange) ,y=log2FoldChange,fill=Class,colour=Class)) + 
  #geom_point(size=20) + 
  geom_bar(stat="identity")+
  geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
  ylab("Log2foldchange")+
  xlab("")+
  #geom_line()+
  geom_hline(yintercept=-0.08, color='red',linetype="dashed")+
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


# mean(subset(test,test$V3=="High feed efficiency")$V5)
#      mean(subset(test,test$V3=="Low feed efficiency")$V5)
#         
#   log2(mean(subset(test,test$V3=="High feed efficiency")$V5)/
#  mean(subset(test,test$V3=="Low feed efficiency")$V5)
#  )
#install.packages("apeglm") 
  # require(apeglm)
  # drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
  # plotMA(res)
########################################################################################3
####Metabolic pathways
  
  
met_df <-   as.data.frame(t(as.data.frame(t(assay(rld)[select,]))))
an_count <-met_df[,as.character(ps@sam_data$Subject)]; 
for(i in 1:ncol(an_count)){
  for(j in 1:nrow(an_count)){
    cat(i,"|",j,"\n")
    if(is.na(an_count[j,i])){ an_count[j,i] <- 0} else if(an_count[j,i]>0){an_count[j,i] <- 1} else {an_count[j,i] <- 0}
  };rm(j)
};rm(i)

animalsID <- unique(ps@sam_data$sampleID)

sample_status <- lapply(1:length(animalsID),function(j){
  
  x <- rowSums(an_count[,as.character(ps@sam_data$Subject[which(ps@sam_data$sampleID==animalsID[[j]])])])
  x[which(x>0)] <- 1
  x <- as.data.frame(x)
  colnames(x) <- animalsID[[j]]
  return(x)
})
sample_status <- do.call(cbind,sample_status)


met_df$animal_per_sample <- NA; met_df$animal_per_sample <- rowSums(sample_status)
met_df$animal_per_HE <- NA; met_df$animal_per_HE <- rowSums(sample_status[,as.character(unique(ps@sam_data$sampleID[which(ps@sam_data$Treatment=="High feed efficiency")]))])
met_df$animal_per_LE <- NA; met_df$animal_per_LE <- rowSums(sample_status[,as.character(unique(ps@sam_data$sampleID[which(ps@sam_data$Treatment=="Low feed efficiency")]))])

#####
met_df$count_per_sample <- NA; met_df$count_per_sample <- rowSums(an_count)
met_df$count_per_HE <- NA; met_df$count_per_HE <- rowSums(an_count[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="High feed efficiency")])])
met_df$count_per_LE <- NA; met_df$count_per_LE <- rowSums(an_count[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="Low feed efficiency")])])
met_df$mean_HE <- NA; met_df$mean_HE <- rowMeans(met_df[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="High feed efficiency")])])
met_df$sd_HE <- NA; met_df$sd_HE <- matrixStats::rowSds(as.matrix(met_df[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="High feed efficiency")])]))
met_df$mean_LE <- NA; met_df$mean_LE <- rowMeans(met_df[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="Low feed efficiency")])])
met_df$sd_LE <- NA; met_df$sd_LE <- matrixStats::rowSds(as.matrix(met_df[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="Low feed efficiency")])]))
#####
met_df$count_per_SO <- NA; met_df$count_per_SO <- rowSums(sample_status[,as.character(ps@sam_data$sampleID[which(ps@sam_data$Phase=="solid")])])
met_df$count_per_LI <- NA; met_df$count_per_LI <- rowSums(sample_status[,as.character(ps@sam_data$sampleID[which(ps@sam_data$Phase=="liquid")])])
met_df$mean_SO <- NA; met_df$mean_SO <- rowMeans(met_df[,as.character(ps@sam_data$Subject[which(ps@sam_data$Phase=="solid")])])
met_df$sd_SO <- NA; met_df$sd_SO <- matrixStats::rowSds(as.matrix(met_df[,as.character(ps@sam_data$Subject[which(ps@sam_data$Phase=="solid")])]))
met_df$mean_LI <- NA; met_df$mean_LI <- rowMeans(met_df[,as.character(ps@sam_data$Subject[which(ps@sam_data$Phase=="liquid")])])
met_df$sd_LI <- NA; met_df$sd_LI <- matrixStats::rowSds(as.matrix(met_df[,as.character(ps@sam_data$Subject[which(ps@sam_data$Phase=="liquid")])]))
#####
met_df$status_eff <- NA; met_df$status_eff <-met_df$mean_HE > met_df$mean_LE
met_df$status_phase <- NA; met_df$status_phase <-met_df$mean_SO > met_df$mean_LI

met_df$status_eff[which(met_df$status_eff==TRUE)] <- "HE"
met_df$status_eff[which(met_df$status_eff==FALSE)] <- "LE"
met_df$status_phase[which(met_df$status_phase==TRUE)] <- "SO"
met_df$status_phase[which(met_df$status_phase==FALSE)] <- "LI"

met_df$status <- NA; met_df$status <-  paste(met_df$status_eff,"-",met_df$status_phase)


met_df$metabolic_pathway <- row.names(met_df)
met_df$metabolic_pathway <-gsub("."," ",met_df$metabolic_pathway,fixed = T)
met_df$KEGG_ID <- NA
met_df$Class_raw <- NA
met_df$Class <- NA
met_df$Description<- NA
met_df$Brite<- NA
met_df$Brite_description<- NA


for(i in 1:nrow(met_df)){
  
  cat(i,"\n")
  
  ##PATHWAYS
  query <- names(keggFind("pathway",query = str_trim(met_df$metabolic_pathway[[i]])))
  
  
  if(is.null(query)){
    met_df$Class_raw[[i]] <- NA
    met_df$KEGG_ID[[i]] <- NA
    met_df$Description[[i]] <- NA
    
  } else {
    #
    met_df$Class_raw[[i]] <- if(is.null(keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS)){met_df$Class_raw[[i]] <- "CHECK" } else { met_df$Class_raw[[i]] <- keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS}
    met_df$Description[[i]] <- if(is.null(keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$DESCRIPTION)){met_df$Description[[i]] <- "CHECK" } else { met_df$Description[[i]] <- paste(keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$DESCRIPTION,collapse=";;;")}
    
    #met_df$Class_raw[[i]] <-# keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS
    if(length(query)>1){
      met_df$KEGG_ID[[i]] <- paste(query[[1]],"*")
    } else {
      met_df$KEGG_ID[[i]] <- query
    }
  }
  ##KO
  query_brite <- names(keggFind("brite",query = str_trim(met_df$metabolic_pathway[[i]])))
  #cat(query_brite,"\n")
  if(is.null(query_brite)){
    
    
    met_df$Brite[[i]] <- NA
    met_df$Brite_description[[i]] <- NA  
    # met_df$Brite[[i]] <- query_brite[[1]]
    # met_df$Brite_description[[i]] <- paste0(keggGet(names(keggFind("brite",query=query_brite[[1]]))),collapse = "|")
    # 
    
  } else {
    
    if(length(grep("br:br",query_brite))>0){
      met_df$Brite[[i]] <- query_brite[[1]]
      met_df$Brite_description[[i]] <- NA
      
    } else {
      an.error.occured <- F
      
      tryCatch( { work <- keggGet(names(keggFind("brite",query=query_brite[[1]])))}
                , error = function(e) {an.error.occured <<- TRUE})
      
      if(an.error.occured==T | is.null(keggGet(names(keggFind("brite",query=query_brite[[1]])))) ){
        met_df$Brite[[i]] <- NA
        met_df$Brite_description[[i]] <- NA 
      } else {
        met_df$Brite[[i]] <- query_brite[[1]]
        met_df$Brite_description[[i]] <- paste0(keggGet(names(keggFind("brite",query=query_brite[[1]]))),collapse = ";;;")
      } 
    }
  }
};rm(i)

#colnames(met_df[,(ps@sam_data[which(ps@sam_data$Treatment=="High feed efficiency"),]$Subject)])


#met_df$ID <- 1:nrow(met_df)



# met_df <- met_df[,c(findcolnumber(met_df,"metabolic_pathway"),
#                     findcolnumber(met_df,"KEGG_ID"),
#                     findcolnumber(met_df,"Class_raw"),
#                     findcolnumber(met_df,"Class"),
#                     findcolnumber(met_df,row.names(ps@otu_table))
#                     )]


met_df <- met_df[,c(71:76,53:70,1:52,77:78)]
met_df2 <- met_df[,1:77]
met_df2 <- met_df2[,c(2:6,77,7:23,1,25:76)]
met_df3 <- (cbind(met_df$metabolic_pathway,met_df$Class_raw,met_df$Description,met_df$Brite, met_df$Brite_description))
write.table(met_df2,paste0(out_dir,"/","csv/","metabolic_pathways.csv"),quote = T,row.names = F,sep="|")
write.csv(met_df2,paste0(out_dir,"/","csv/","metabolic_pathways.csv"),quote = T,row.names = F)
write.table(met_df3,paste0(out_dir,"/","csv/","metabolic_pathways_brite.csv"),quote = T,row.names = F,sep="|")


#######################################################################################
########################################################################################
met_df_1 <- met_df
met_df_1 <- met_df_1 %>% 
  dplyr::filter(!str_detect(met_df_1$Class_raw, 'Human Diseases'))

met_df_2 <- met_df
met_df_2 <- met_df_2[which(is.na(met_df_2$Class_raw)),]

met_df_filtered <- rbind(met_df_1,met_df_2)

met_df_filtered2 <- met_df_filtered[,1:77]
met_df_filtered2 <- met_df_filtered2[,c(2:6,77,7:23,1,25:76)]

write.csv(met_df_filtered,paste0(out_dir,"/","csv/","metabolic_pathways.csv_filtered.csv"),quote = T,row.names = F)
write.csv(met_df_filtered2,paste0(out_dir,"/","csv/","metabolic_pathways.csv_filtered_NB.csv"),quote = T,row.names = F)


met_df_filtered3 <- as.data.frame(t(met_df_filtered2[,25:76]))
colnames(met_df_filtered3) <- met_df_filtered2$metabolic_pathway
rownames(met_df_filtered3) <- colnames(met_df_filtered2[,25:76])
met_df_filtered3 <- met_df_filtered3[complete.cases(met_df_filtered3),]

library(usdm)
library(plsdepot)
z2 = plsdepot::nipals(Data = met_df_filtered3, comps=5, scaled=F)

require(rgl);require(car)
#scatter3d(x = PCA_DES$PC1, y = PCA_DES$PC2, z = PCA_DES$PC3, 
scatter3d(x = z2$scores[,1], y = z2$scores[,2], z = z2$scores[,3], 
          groups =  factor(paste0(colData1$Treatment,"-",colData1$Phase),levels = 
                             c("High feed efficiency-solid","High feed efficiency-liquid","Low feed efficiency-solid", "Low feed efficiency-liquid")),
          # xlab= paste("Axis 1",round(percentVar[1]*100,2),"%"),
          # ylab= paste("Axis 2",round(percentVar[2]*100,2),"%"),
          # zlab= paste("Axis 3",round(percentVar[3]*100,2),"%"),
          xlab= paste("Axis 1",round(z2$values$percentage[1],2)*100,"%"),
          ylab= paste("Axis 2",round(z2$values$percentage[2],2)*100,"%"),
          zlab= paste("Axis 3",round(z2$values$percentage[3],2)*100,"%"),
          
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
########################################################################################
########################################################################################

met_df2_rel_abund <- met_df2[,c(1,3,6,25:76)]

for(i in 4:ncol(met_df2_rel_abund)){
  met_df2_rel_abund[,i] <- met_df2_rel_abund[,i]/sum(met_df2_rel_abund[,i])
};rm(i)


####
#####
met_df2_rel_abund$count_per_sample <- NA; met_df2_rel_abund$count_per_sample <- rowSums(an_count)
met_df2_rel_abund$count_per_HE <- NA; met_df2_rel_abund$count_per_HE <- rowSums(an_count[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="High feed efficiency")])])
met_df2_rel_abund$count_per_LE <- NA; met_df2_rel_abund$count_per_LE <- rowSums(an_count[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="Low feed efficiency")])])
met_df2_rel_abund$mean_HE <- NA; met_df2_rel_abund$mean_HE <- rowMeans(met_df2_rel_abund[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="High feed efficiency")])])
met_df2_rel_abund$sd_HE <- NA; met_df2_rel_abund$sd_HE <- matrixStats::rowSds(as.matrix(met_df2_rel_abund[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="High feed efficiency")])]))
met_df2_rel_abund$mean_LE <- NA; met_df2_rel_abund$mean_LE <- rowMeans(met_df2_rel_abund[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="Low feed efficiency")])])
met_df2_rel_abund$sd_LE <- NA; met_df2_rel_abund$sd_LE <- matrixStats::rowSds(as.matrix(met_df2_rel_abund[,as.character(ps@sam_data$Subject[which(ps@sam_data$Treatment=="Low feed efficiency")])]))
#####
met_df2_rel_abund$count_per_SO <- NA; met_df2_rel_abund$count_per_SO <- rowSums(sample_status[,as.character(ps@sam_data$sampleID[which(ps@sam_data$Phase=="solid")])])
met_df2_rel_abund$count_per_LI <- NA; met_df2_rel_abund$count_per_LI <- rowSums(sample_status[,as.character(ps@sam_data$sampleID[which(ps@sam_data$Phase=="liquid")])])
met_df2_rel_abund$mean_SO <- NA; met_df2_rel_abund$mean_SO <- rowMeans(met_df2_rel_abund[,as.character(ps@sam_data$Subject[which(ps@sam_data$Phase=="solid")])])
met_df2_rel_abund$sd_SO <- NA; met_df2_rel_abund$sd_SO <- matrixStats::rowSds(as.matrix(met_df2_rel_abund[,as.character(ps@sam_data$Subject[which(ps@sam_data$Phase=="solid")])]))
met_df2_rel_abund$mean_LI <- NA; met_df2_rel_abund$mean_LI <- rowMeans(met_df2_rel_abund[,as.character(ps@sam_data$Subject[which(ps@sam_data$Phase=="liquid")])])
met_df2_rel_abund$sd_LI <- NA; met_df2_rel_abund$sd_LI <- matrixStats::rowSds(as.matrix(met_df2_rel_abund[,as.character(ps@sam_data$Subject[which(ps@sam_data$Phase=="liquid")])]))
#####
met_df2_rel_abund$status_eff <- NA; met_df2_rel_abund$status_eff <-met_df2_rel_abund$mean_HE > met_df2_rel_abund$mean_LE
met_df2_rel_abund$status_phase <- NA; met_df2_rel_abund$status_phase <-met_df2_rel_abund$mean_SO > met_df2_rel_abund$mean_LI

met_df2_rel_abund$status_eff[which(met_df2_rel_abund$status_eff==TRUE)] <- "HE"
met_df2_rel_abund$status_eff[which(met_df2_rel_abund$status_eff==FALSE)] <- "LE"
met_df2_rel_abund$status_phase[which(met_df2_rel_abund$status_phase==TRUE)] <- "SO"
met_df2_rel_abund$status_phase[which(met_df2_rel_abund$status_phase==FALSE)] <- "LI"
met_df2_rel_abund$status <- NA; met_df2_rel_abund$status <- paste(met_df2_rel_abund$status_eff,"-",met_df2_rel_abund$status_phase)

write.table(met_df2_rel_abund,paste0(out_dir,"/","csv/","metabolic_pathways_rel_abund.csv"),quote = T,row.names = F,sep=",")

########################################################################################

########################################################################################
########################################################################################








####Correlation of pathways and production data
  cowPi_assoc <- microbiome::associate(as.data.frame(t(assay(rld)[select,])), CowPI_file2a[,c(5:9)], 
                                       method = "spearman",
                                       mode = "table", p.adj.threshold = 0.05, n.signif = 1,
                                       p.adj.method = "fdr")
  
#####Attaching KEGG and cohort information
  cowPi_assoc_plot <- microbiome::heat(cowPi_assoc, "X1", "X2", fill = "Correlation", colours = c("darkred","red","white","blue","darkblue"),
       star = "p.adj", p.adj.threshold = 0.05,plot.values = F,star.size=40,legend.text = "Spearman"
       
  )+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black")) +
  #scale_fill_continuous(guide = guide_colorbar(direction = "horizontal")) +
  theme(legend.position="right",legend.direction = "vertical",legend.key.size =  unit(1.5, "in"))
  
  
  ggsave(paste0(graph_dir,"/","heat_met_abund",".pdf"),cowPi_assoc_plot,dpi=300,width =160,height=130,units = "cm",scale=1.2,limitsize = FALSE)
  
  
  
  
  #write.csv(cowPi_assoc,paste0(out_dir,"/","csv/","correlation_met_cowpi_deseq2.csv"),quote = F,row.names = F)
  
  #################################
  cowPi_assoc$metabolic_pathway <- NA
  cowPi_assoc$metabolic_pathway <- gsub("."," ",cowPi_assoc$X1,fixed = T)
  ####################
  
  # cowPi_assoc$KEGG_ID <- NA
  # cowPi_assoc$Class_raw <- NA
  # cowPi_assoc$Class <- NA
  # cowPi_assoc$Description <- NA
  # cowPi_assoc$Brite <- NA
  
  ########################################################
  ########################################################
  # colnames(cowPi_assoc)[1] <- "metabolic_pathway"  
  # colnames(cowPi_assoc)[2] <- "variable"  
  cowPi_assoc$KEGG_ID <- dplyr::left_join(cowPi_assoc,met_df2,"metabolic_pathway")$KEGG_ID.y
  cowPi_assoc$Class_raw <- dplyr::left_join(cowPi_assoc,met_df2,"metabolic_pathway")$Class_raw.y
  cowPi_assoc$Class <- dplyr::left_join(cowPi_assoc,met_df2,"metabolic_pathway")$Class.y
  cowPi_assoc$Description <- dplyr::left_join(cowPi_assoc,met_df2,"metabolic_pathway")$Description.y
  cowPi_assoc$Brite <- dplyr::left_join(cowPi_assoc,met_df2,"metabolic_pathway")$Brite.y
  cowPi_assoc$status <- dplyr::left_join(cowPi_assoc,met_df2,"metabolic_pathway")$status.y
  ########################################################  
  ########################################################
  
  #importance <- importance[which(importance$Class_raw !="Human Diseases; Infectious diseases: Bacterial"),]
  cowPi_assoc$p_sig <- NA
  cowPi_assoc$p_sig[which(cowPi_assoc$p.adj<0.05)]<-"*"
  write.csv(cowPi_assoc,paste0(out_dir,"/","csv/","correlation_met_cowpi_deseq2.csv"),quote = T,row.names = F)
  
  cowPi_assoc2 <- cowPi_assoc
  cowPi_assoc2 <- cowPi_assoc2 %>% 
    dplyr::filter(!str_detect(cowPi_assoc2$Class_raw, 'Human Diseases'))
  
  cowPi_assoc3 <- cowPi_assoc
  cowPi_assoc3 <- cowPi_assoc3[which(is.na(cowPi_assoc3$Class_raw)),]
  
  cowPi_assoc_filtered <- rbind(cowPi_assoc2,cowPi_assoc3)
  write.csv(cowPi_assoc_filtered,paste0(out_dir,"/","csv/","correlation_met_cowpi_deseq2_filtered.csv"),quote = T,row.names = F)
  
  #################################
  
  
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
  CowPI_file6 = CowPI_file6,
  met_df = met_df,
  met_df_filtered3= met_df_filtered3,
  cowPi_assoc=cowPi_assoc)
  saveRDS(cowPI,paste0(out_dir,"/","csv/","cowPI_R_FCR.RDS"))
  
  
  
  if(!file.exists(paste0(out_dir,"/","csv/","cowPI_step_list_des.RDS"))){
    # Set up repeated k-fold cross-validation
    set.seed(300)
    
    
    library(stringr)
    
    met_df2 <- met_df
    met_df2 <- met_df2 %>% 
      dplyr::filter(!str_detect(met_df2$Class_raw, 'Human Diseases'))
    
    met_df3 <- met_df
    met_df3 <- met_df3[which(is.na(met_df3$Class_raw)),]
    
    met_df_filtered <- rbind(met_df2,met_df3)
    
    
    met_df_filtered$PA <- met_df_filtered$metabolic_pathway 
    met_df_filtered$PA <- gsub(" ",".",met_df_filtered$PA,fixed = T)
    
    CowPI_file7 <- CowPI_file6
    
    CowPI_file7 <- subset(CowPI_file7, select=met_df_filtered$PA)
    CowPI_file7$FCR <- CowPI_file6$FCR
    
   
    
    train.control <- trainControl(method = "cv", number = 10)
    
    step.model <- train(FCR ~., data = CowPI_file7,
                        method = "pls", #="leapseq"
                        preproc= c("zv","center","scale"),
                        trControl=train.control,
                        #metric="Accuracy",
                        # trControl=(trainControl(method="repeatedcv", repeats=3,
                        # classProbs=TRUE, summaryFunction=twoClassSummary)),
                        tuneLength=30
                      # metric="ROC"
                       #tuneGrid = data.frame(nvmax = 1:10),
#                        trControl = train.control
    )
    plot(step.model$finalModel, ncomp=as.numeric(step.model$bestTune))
    
    confusionMatrix(step.model,pred)
    plot(step.model)

            step.model$results
    step.model$bestTune
    summary(step.model$finalModel)
    
    coef(step.model$finalModel,as.numeric(step.model$bestTune))
    plot(step.model$finalModel)
    summary(step.model$finalModel)
    
    pred <- predict(step.model, newdata = CowPI_file7)
    #acc
    pred_df <- data.frame(Subject=row.names(CowPI_file7) ,pred = pred, FCR = CowPI_file7$FCR)
    pred_df$Subject <- as.character(pred_df$Subject)
    #row.names(pred_df) <- row.names(CowPI_file6)
    pred_df$Treatment <- NA; pred_df$Treatment <- join(pred_df,grp,"Subject")$grp
    
    library(ggplot2)
    
    # ggplot(data = pred_df)+
    #   geom_line(aes(x = FCR, y = pred_df))+
    #   geom_point(aes(x = FCR, y = pred_df))
    caret::R2(pred_df$pred,pred_df$FCR)
    RMSE(pred = pred_df$pred,obs = pred_df$FCR)
    
    coef(step.model$finalModel,5)
    
  
    require(ggpubr)
  
    importance <- varImp(step.model)$importance#arRF$
    importance$PA <- row.names(importance)
    importance <- importance[,c(2,1)]
    importance <- importance[order(importance$Overall,decreasing = TRUE),]
    importance <- importance[which(importance$Overall>=50),]
    importance$metabolic_pathway <- importance$PA 
    importance$metabolic_pathway <- gsub("."," ",importance$metabolic_pathway,fixed = T)
    ####################
    
    importance$KEGG_ID <- NA
    importance$Class_raw <- NA
    importance$Class <- NA
    for(i in 1:nrow(importance)){
      cat(i,"\n")
      query <- names(keggFind("pathway",query = str_trim(importance$metabolic_pathway[[i]])))
      if(is.null(query)){
        importance$Class_raw[[i]] <- NA
        importance$KEGG_ID[[i]] <- NA
      } else {
        importance$Class_raw[[i]] <- if(is.null(keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS)){importance$Class_raw[[i]] <- NA } else { importance$Class_raw[[i]] <- keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS}
        
        #importance$Class_raw[[i]] <- keggGet(names(keggFind("pathway",query=query[[1]])))[[1]]$CLASS
        if(length(query)>1){
          importance$KEGG_ID[[i]] <- paste(query[[1]],"*")
        } else {
          importance$KEGG_ID[[i]] <- query
        }
      }
    };rm(i)
    
    #importance <- importance[which(importance$Class_raw !="Human Diseases; Infectious diseases: Bacterial"),]
    
    ###################
   
    
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
    write.csv(importance,paste0(out_dir,"/","csv/","pls_importance.csv"),quote = T,row.names = F)
    
    
  } else {
    
    cowPI_step <- readRDS(paste0(out_dir,"/","csv/","cowPI_step_list.RDS"))
  }
  
 