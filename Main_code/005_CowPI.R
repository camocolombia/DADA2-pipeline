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
#library("microbiomeViz");
library("stringr")
library("StructFDR")
library("gtools")
require("agricolae")
require("car")
#install_github("JiyuanHu/massMap")

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


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


ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))
ps3 <-  readRDS(paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))

CowPI_file <- read.csv(paste0(CowPI_dir,"/","cowPI2.csv"),header = T,sep = "|")
CowPI_file <- CowPI_file[,which(colSums(CowPI_file)>0)]
CowPI_file2 <- merge(meta(ps),CowPI_file,by.x="Subject",by.y="samples")
CowPI_file2a <-CowPI_file2

#CowPI_file2 <- cbind(CowPI_file2[,c(1:11)],log10(CowPI_file2[,-c(1:11)]))
cowPi_assoc <- microbiome::associate((CowPI_file2[,-c(1:10)]), CowPI_file2[,c(5:9)], 
                      method = "spearman",
                      mode = "table", p.adj.threshold = 0.05, n.signif = 1,
                      p.adj.method = "fdr")


heat(cowPi_assoc, "X1", "X2", fill = "Correlation", colours = c("darkred","red","white","blue","darkblue"),
     star = "p.adj", p.adj.threshold = 0.05) 


cowPi_assoc_FCR <- cowPi_assoc[which(cowPi_assoc$X2=="FCR" & cowPi_assoc$p.adj <= 0.05),]

wilcox_cowPI <- lapply(1:ncol(CowPI_file2[,-c(1:10)]),function(i){ 
cat(i,"\n")
    x <- subset(CowPI_file2,CowPI_file2$Treatment== "Low feed efficiency")
  x <- x[,colnames(CowPI_file2[,-c(1:10)])[i]]
  y <- subset(CowPI_file2,CowPI_file2$Treatment== "High feed efficiency")
  y <- y[,colnames(CowPI_file2[,-c(1:10)])[i]]
  
  x2 <- subset(CowPI_file2,CowPI_file2$Phase== "liquid")
  x2 <- x2[,colnames(CowPI_file2[,-c(1:10)])[i]]
  y2 <- subset(CowPI_file2,CowPI_file2$Phase== "solid")
  y2 <- y2[,colnames(CowPI_file2[,-c(1:10)])[i]]
  
    x <- data.frame(
    Pathway = colnames(CowPI_file2[,-c(1:10)])[i],
    Treatment = wilcox.test(x=x,y=y,paired=F,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value,
    Treatment_padj = NA,
     median_LE = median(x,na.rm=T),
    mean_LE = mean(x,na.rm=T),
    sd_LE = sd(x,na.rm=T),
    sem_LE = sd(x)/sqrt(length(x)),
    medianHE = median(y,na.rm=T),
    meanHE = mean(y,na.rm=T),
    sdHE = sd(y,na.rm=T),
    sem_HE = sd(y)/sqrt(length(y)),
    Phase = wilcox.test(x=x2,y=y2,paired=T,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value,
    Phase_padj = NA,
    median_liquid = median(x2,na.rm=T),
    mean_liquid = mean(x2,na.rm=T),
    sd_liquid = sd(x2,na.rm=T),
    sem_liquid = sd(x2)/sqrt(length(x2)),
    median_solid = median(y2,na.rm=T),
    mean_solid = mean(y2,na.rm=T),
    sd_solid = sd(y2,na.rm=T),
    sem_liquid = sd(y2)/sqrt(length(y2))
)
  return(x)
  
  })

wilcox_cowPI <- do.call(rbind, wilcox_cowPI)
wilcox_cowPI$Treatment_padj <- p.adjust(p = wilcox_cowPI$Treatment, method =  "BH",n = length(wilcox_cowPI$Treatment))
wilcox_cowPI$Phase_padj <- p.adjust(p = wilcox_cowPI$Phase, method =  "BH",n = length(wilcox_cowPI$Phase))


wilcox_cowPI2 <- wilcox_cowPI[which(wilcox_cowPI$Phase_padj <= 0.05),]

require(FactoMineR);require(factoextra)

CowPI_file2$groups <- factor(paste0(CowPI_file2$Treatment," - ",CowPI_file2$Phase),
                             levels =c("High feed efficiency - solid","High feed efficiency - liquid",
                                       "Low feed efficiency - solid","Low feed efficiency - liquid"))

CowPI_file3 <- CowPI_file2[,names(colSums(CowPI_file2[,-c(1:10,263)])
                                  [which(colSums(CowPI_file2[,-c(1:10,263)])
                                         >= quantile(colSums(CowPI_file2[,-c(1:10,263)]),0.1))])]

CowPI_file3$groups <- CowPI_file2$groups
CowPI_file4 <- CowPI_file3
for(i in 1:nrow(CowPI_file4)){
  CowPI_file4[i,-ncol(CowPI_file4)] <- as.numeric(CowPI_file4[i,-ncol(CowPI_file4)])/as.numeric(rowSums(CowPI_file4[i,-ncol(CowPI_file4)]))
};rm(i)


file_to_PCA <- CowPI_file4#CowPI_file2[,-c(1,2,3,4,10)]
file_to_PCA <- CowPI_file4[,-c(ncol(CowPI_file4))]#file_to_PCA[,-c(1:5)]
row.names(file_to_PCA) <- paste0(CowPI_file2$sampleID," - ",CowPI_file2$Phase)


PCA <- FactoMineR::PCA(file_to_PCA,scale.unit = TRUE,ncp=5)#quali.sup = 251  ,ncp = 5)
grp <-  factor(paste0(CowPI_file2$Treatment," - ",CowPI_file2$Phase),
               levels =c("High feed efficiency - solid","High feed efficiency - liquid",
                         "Low feed efficiency - solid","Low feed efficiency - liquid"))

fviz_pca_ind(PCA,
             label = "ind",# "ind", # hide individual labels
             habillage = grp,
             col.var = "black",#CowPI_file2$Treatment, # color by groups
             palette = c("darkblue", "green", "orange","red"),
             addEllipses = TRUE,
             ellipse.type = "convex",
             repel=T# Concentration ellipses
)



###

scatter3d(x = PCA$ind$coord[,1], y = PCA$ind$coord[,2], z = PCA$ind$coord[,3], 
          groups = grp ,
          xlab= paste("Axis 1",round(PCA$eig[1,2],2),"%"),
          ylab= paste("Axis 2",round(PCA$eig[2,2],2),"%"),
          zlab= paste("Axis 3",round(PCA$eig[3,2],2),"%"),
          surface=FALSE, grid = FALSE, ellipsoid = TRUE,
          labels = row.names(PCA$ind$coord),
          point.col =  c("black", "black", "black"),id.n=length(row.names(PCA$ind$coord)),
          surface.col = c("darkblue", "green", "orange","red"),
          axis.col = c("black", "black", "black")
)
require(car);require(rgl)
legend3d("topright", legend = c("High feed efficiency - solid","High feed efficiency - liquid",
                                "Low feed efficiency - solid","Low feed efficiency - liquid"),
         pch = 16, col = c("darkblue", "green", "orange","red"), cex=2, inset=c(0.01))
fviz_pca_var(PCA, col.var = "black")


res.desc <- dimdesc(PCA, axes = c(1:5), proba = 0.05)
# library("corrplot")
# #
# corrplot(cor(file_to_PCA))
# View(cor(file_to_PCA))
# require(caret)
              
#>= quantile(colSums(CowPI_file2[,-c(1:11,261)]),0.1))]

CowPI_file5 <- t(CowPI_file4)
colnames(CowPI_file5) <- row.names(CowPI_file4)


x <-CowPI_file3[,-c(1:10,227)]
x<-x/rowSums(x)
x<-x[,order(colSums(x),decreasing=TRUE)]
#Extract list of top N Taxa
N<-21
taxa_list<-colnames(x)[1:N]
taxa_list <- taxa_list[which(!is.na(taxa_list))]

new_x<-data.frame(x[,colnames(x) %in% taxa_list],
                  
                  Others= if(length(colnames(x))<N){Others=NA} else { rowSums(x[,!colnames(x) %in% taxa_list])}
)
df<-NULL
for (i in 1:dim(new_x)[2]){
  cat(i,"\n")
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),sampleID = CowPI_file2$sampleID,
                  Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],
                  Type=CowPI_file2$Treatment,Phase= CowPI_file2$Phase)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
};rm(i)


colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00");
#N <- 26
library(ggplot2)
p<-ggplot(df,aes(factor(sampleID),Value,fill=Taxa))+
  geom_bar(stat="identity")+
  facet_grid(. ~ Type + Phase, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:(N+1)])
#p<-p+scale_fill_manual(labs="Others",values="white")
p<-p+theme_bw()+ylab("Proportions")
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines"))
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  xlab("Sample")
pdf("TAXAplot.pdf",height=6,width=21)
print(p)
dev.off()
p


# qplot(colSums(CowPI_file3[,-226]))+
#   geom_bar()

sum_CowPI <- data.frame(
Pathways =names(colSums(subset(CowPI_file3[,-227]))),
HE_S_Raw = colSums(subset(CowPI_file3, groups=="High feed efficiency - solid")[,-227]),
HE_L_Raw  = colSums(subset(CowPI_file3, groups=="High feed efficiency - liquid")[,-227]),
LE_S_Raw = colSums(subset(CowPI_file3, groups=="Low feed efficiency - solid")[,-227]),
LE_L_Raw = colSums(subset(CowPI_file3, groups=="Low feed efficiency - liquid")[,-227]),
Raw = colSums(subset(CowPI_file3[,-227])),
HE_S = colSums(subset(CowPI_file3, groups=="High feed efficiency - solid")[,-227])/sum(colSums(CowPI_file3[,-227])),
HE_L = colSums(subset(CowPI_file3, groups=="High feed efficiency - liquid")[,-227])/sum(colSums(CowPI_file3[,-227])),
LE_S = colSums(subset(CowPI_file3, groups=="Low feed efficiency - solid")[,-227])/sum(colSums(CowPI_file3[,-227])),
LE_L = colSums(subset(CowPI_file3, groups=="Low feed efficiency - liquid")[,-227])/sum(colSums(CowPI_file3[,-227])),
Total =colSums(subset(CowPI_file3[,-227]))/sum(colSums(CowPI_file3[,-227]))
)
row.names(sum_CowPI) <- 1:nrow(sum_CowPI)
#colSums(CowPI_file3[,-226])

#####################################
##Wilcoxon using proportional relative abundances
wilcox_cowPI_R <- lapply(1:ncol(CowPI_file4[,-c(227)]),function(i){ 
  cat(i,"\n")
  x <- subset(CowPI_file4,CowPI_file2$Treatment== "Low feed efficiency")
  x <- x[,colnames(CowPI_file4[,-c(227)])[i]]
  y <- subset(CowPI_file4,CowPI_file2$Treatment== "High feed efficiency")
  y <- y[,colnames(CowPI_file4[,-c(227)])[i]]
  #
  x2 <- subset(CowPI_file4,CowPI_file2$Phase== "liquid")
  x2 <- x2[,colnames(CowPI_file4[,-c(227)])[i]]
  y2 <- subset(CowPI_file4,CowPI_file2$Phase== "solid")
  y2 <- y2[,colnames(CowPI_file4[,-c(227)])[i]]
  #
  x <- data.frame(
    Pathway = colnames(CowPI_file4[,-c(227)])[i],
    Treatment = wilcox.test(x=x,y=y,paired=F,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value,
    Treatment_padj = NA,
    median_LE = median(x,na.rm=T),
    mean_LE = mean(x,na.rm=T),
    gmean_LE = gm_mean(x),
    sd_LE = sd(x,na.rm=T),
    sem_LE = sd(x)/sqrt(length(x)),
    median_HE = median(y,na.rm=T),
    mean_HE = mean(y,na.rm=T),
    gmean_HE = gm_mean(y),
    sd_HE = sd(y,na.rm=T),
    sem_HE = sd(y)/sqrt(length(y)),
    fc_HE_LE = NA,
    Phase = wilcox.test(x=x2,y=y2,paired=T,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value,
    Phase_padj = NA,
    median_liquid = median(x2,na.rm=T),
    mean_liquid = mean(x2,na.rm=T),
    gmean_liquid = gm_mean(x2),
    sd_liquid = sd(x2,na.rm=T),
    sem_liquid = sd(x2)/sqrt(length(x2)),
    median_solid = median(y2,na.rm=T),
    mean_solid = mean(y2,na.rm=T),
    gmean_solid = gm_mean(y2),
    sd_solid = sd(y2,na.rm=T),
    sem_solid = sd(y2)/sqrt(length(y2)),
    fc_SO_LI = NA
  )
  return(x)
  
})

wilcox_cowPI_R <- do.call(rbind, wilcox_cowPI_R)

##############
##LOG2FOLDCHANGE


wilcox_cowPI_R$fc_HE_LE <- log2(wilcox_cowPI_R$gmean_HE/wilcox_cowPI_R$gmean_LE);
wilcox_cowPI_R$fc_SO_LI <- log2(wilcox_cowPI_R$gmean_solid/wilcox_cowPI_R$gmean_liquid);

##prop count
for(i in 1:nrow(wilcox_cowPI_R)){
  if(is.infinite(wilcox_cowPI_R$fc_HE_LE[i])){
    if(wilcox_cowPI_R$gmean_HE[i]>0){
      wilcox_cowPI_R$fc_HE_LE[i] <- wilcox_cowPI_R$gmean_HE
    } else { wilcox_cowPI_R$fc_HE_LE[i] <- -(wilcox_cowPI_R$gmean_LE)}
  } else {wilcox_cowPI_R$fc_HE_LE[i] <- wilcox_cowPI_R$fc_HE_LE[i] }
};rm(i)

for(i in 1:nrow(wilcox_cowPI_R)){
  if(is.infinite(wilcox_cowPI_R$fc_SO_LI[i])){
    if(wilcox_cowPI_R$gmean_solid[i]>0){
      wilcox_cowPI_R$fc_SO_LI[i] <- wilcox_cowPI_R$gmean_solid
    } else { wilcox_cowPI_R$fc_SO_LI[i] <- -(wilcox_cowPI_R$gmean_liquid)}
  } else {wilcox_cowPI_R$fc_SO_LI[i] <- wilcox_cowPI_R$fc_SO_LI[i] }
};rm(i)
#############################
#############

wilcox_cowPI_R$Treatment_padj <- p.adjust(p = wilcox_cowPI_R$Treatment, method =  "BH",n = length(wilcox_cowPI_R$Treatment))
wilcox_cowPI_R$Phase_padj <- p.adjust(p = wilcox_cowPI_R$Phase, method =  "BH",n = length(wilcox_cowPI_R$Phase))


wilcox_cowPI2 <- wilcox_cowPI_R[which(wilcox_cowPI_R$Treatment_padj <= 0.05),]
wilcox_cowPI2_Phase <- wilcox_cowPI_R[which(wilcox_cowPI_R$Phase_padj <= 0.05),]


write.csv(wilcox_cowPI2,paste0(out_dir,"/","csv/","wilcoxon_met.csv"),quote = F,row.names = F)

############plot

#Log2 fold change plot
exp_plot2 <- ggplot(wilcox_cowPI2, aes(y=fc_HE_LE, x=reorder(Pathway, -fc_HE_LE),
                                       color=NULL, 
                                       fill=NULL)) + 
 # geom_point(size=25) + 
  geom_bar(stat="identity")+
  geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
  ylab("Log2foldchange")+
  xlab("")+
  #geom_line()+
  theme(text=element_text(size=69),
        legend.text=element_text(face="italic",size=60,colour="black"),
        axis.text.x  = element_text(face="italic",angle = -90,size=60, hjust = 0, vjust=0.5),
        axis.text.y  = element_text(face="italic",size=60,colour="black"))+
  scale_y_continuous(limits = c(-2,2),breaks=seq(-2,2,0.2))+
  coord_flip()
#exp_plot2  
ggsave(paste0(graph_dir,"/","log2foldchange_pathways",".pdf"),exp_plot2,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)

exp_plot2
####
exp_plot2_ph <- ggplot(wilcox_cowPI2_Phase, aes(y=fc_SO_LI, x=reorder(Pathway, -fc_SO_LI),
                                       color=NULL, 
                                       fill=NULL)) + 
  # geom_point(size=25) + 
  geom_bar(stat="identity")+
  geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
  ylab("Log2foldchange")+
  xlab("")+
  #geom_line()+
  theme(text=element_text(size=69),
        legend.text=element_text(face="italic",size=60,colour="black"),
        axis.text.x  = element_text(face="italic",angle = -90,size=60, hjust = 0, vjust=0.5),
        axis.text.y  = element_text(face="italic",size=60,colour="black"))+
  scale_y_continuous(limits = c(-2,2),breaks=seq(-2,2,0.2))+
  coord_flip()
#exp_plot2  
ggsave(paste0(graph_dir,"/","log2foldchange_pathway_phase",".pdf"),exp_plot2_ph,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)


###########################################3
###SOLIDS

wilcox_cowPI_f_ph <- lapply(1:ncol(CowPI_file4[,-c(227)]),function(i){ 
  cat(i,"\n")
  x <- subset(CowPI_file4,CowPI_file2$Treatment== "Low feed efficiency" &
                CowPI_file2$Phase== "liquid" )
  x <- x[,colnames(CowPI_file4[,-c(227)])[i]]
  y <-    subset(CowPI_file4,CowPI_file2$Treatment== "High feed efficiency" &
             CowPI_file2$Phase== "liquid" )
  y <- y[,colnames(CowPI_file4[,-c(227)])[i]]
  #
  x2 <- subset(CowPI_file4,CowPI_file2$Treatment== "Low feed efficiency" &
                 CowPI_file2$Phase== "solid" ) 
  x2 <- x2[,colnames(CowPI_file4[,-c(227)])[i]]
  y2 <- subset(CowPI_file4,CowPI_file2$Treatment== "High feed efficiency" &
                 CowPI_file2$Phase== "solid" ) 
  y2 <- y2[,colnames(CowPI_file4[,-c(227)])[i]]
  #
  x <- data.frame(
    Pathway = colnames(CowPI_file4[,-c(227)])[i],
    Treatment = wilcox.test(x=x,y=y,paired=F,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value,
    Treatment_padj = NA,
    median_LE = median(x,na.rm=T),
    mean_LE = mean(x,na.rm=T),
    gmean_LE = gm_mean(x),
    sd_LE = sd(x,na.rm=T),
    sem_LE = sd(x)/sqrt(length(x)),
    median_HE = median(y,na.rm=T),
    mean_HE = mean(y,na.rm=T),
    gmean_HE = gm_mean(y),
    sd_HE = sd(y,na.rm=T),
    sem_HE = sd(y)/sqrt(length(y)),
    fc_HE_LE = NA,
    Phase = wilcox.test(x=x2,y=y2,paired=T,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value,
    Phase_padj = NA,
    median_liquid = median(x2,na.rm=T),
    mean_liquid = mean(x2,na.rm=T),
    gmean_liquid = gm_mean(x2),
    sd_liquid = sd(x2,na.rm=T),
    sem_liquid = sd(x2)/sqrt(length(x2)),
    median_solid = median(y2,na.rm=T),
    mean_solid = mean(y2,na.rm=T),
    gmean_solid = gm_mean(y2),
    sd_solid = sd(y2,na.rm=T),
    sem_solid = sd(y2)/sqrt(length(y2)),
    fc_SO_LI = NA
  )
  return(x)
  
})

wilcox_cowPI_f_ph <- do.call(rbind, wilcox_cowPI_f_ph)

wilcox_cowPI_f_ph$fc_HE_LE <- log2(wilcox_cowPI_f_ph$gmean_HE/wilcox_cowPI_f_ph$gmean_LE);
wilcox_cowPI_f_ph$fc_SO_LI <- log2(wilcox_cowPI_f_ph$gmean_solid/wilcox_cowPI_f_ph$gmean_liquid);

##prop count
for(i in 1:nrow(wilcox_cowPI_f_ph)){
  if(is.infinite(wilcox_cowPI_f_ph$fc_HE_LE[i])){
    if(wilcox_cowPI_f_ph$gmean_HE[i]>0){
      wilcox_cowPI_f_ph$fc_HE_LE[i] <- wilcox_cowPI_f_ph$gmean_HE
    } else { wilcox_cowPI_f_ph$fc_HE_LE[i] <- -(wilcox_cowPI_f_ph$gmean_LE)}
  } else {wilcox_cowPI_f_ph$fc_HE_LE[i] <- wilcox_cowPI_f_ph$fc_HE_LE[i] }
};rm(i)

for(i in 1:nrow(wilcox_cowPI_f_ph)){
  if(is.infinite(wilcox_cowPI_f_ph$fc_SO_LI[i])){
    if(wilcox_cowPI_f_ph$gmean_solid[i]>0){
      wilcox_cowPI_f_ph$fc_SO_LI[i] <- wilcox_cowPI_f_ph$gmean_solid
    } else { wilcox_cowPI_f_ph$fc_SO_LI[i] <- -(wilcox_cowPI_f_ph$gmean_liquid)}
  } else {wilcox_cowPI_f_ph$fc_SO_LI[i] <- wilcox_cowPI_f_ph$fc_SO_LI[i] }
};rm(i)
#############################
#############

wilcox_cowPI_f_ph$Treatment_padj <- p.adjust(p = wilcox_cowPI_f_ph$Treatment, method =  "BH",n = length(wilcox_cowPI_f_ph$Treatment))
wilcox_cowPI_f_ph$Phase_padj <- p.adjust(p = wilcox_cowPI_f_ph$Phase, method =  "BH",n = length(wilcox_cowPI_f_ph$Phase))


wilcox_cowPI_f_ph2 <- wilcox_cowPI_f_ph[which(wilcox_cowPI_f_ph$Treatment_padj <= 0.05),]
wilcox_cowPI_f_ph_phase <- wilcox_cowPI_f_ph[which(wilcox_cowPI_f_ph$Phase_padj <= 0.05),]



#########################################33



cowPi_assoc <- microbiome::associate((CowPI_file4[,-227]), CowPI_file2[,c(5:9)], 
                                     method = "spearman",
                                     mode = "table", p.adj.threshold = 0.05, n.signif = 1,
                                     p.adj.method = "fdr")
heat(cowPi_assoc, "X1", "X2", fill = "Correlation", colours = c("darkred","red","white","blue","darkblue"),
     star = "p.adj", p.adj.threshold = 0.05) 


write.csv(cowPi_assoc,paste0(out_dir,"/","csv/","correlation_met.csv"),quote = F,row.names = F)

sum_CowPI <- sum_CowPI[order(sum_CowPI$Total,decreasing = T),]
View(sum_CowPI)
CowPI_file4$Treatment <- CowPI_file2$Treatment

ado_comp <- adonis(CowPI_file4[,1:226] ~ Treatment, data = CowPI_file4,permutations = 10000)
anova(betadisper(vegdist(CowPI_file4[,1:226],"mahalanobis"),CowPI_file4$Treatment,bias.adjust = T))
anos <- anosim(CowPI_file4[,1:226],grouping = CowPI_file4$Treatment,distance="gower",permutations = 1000)
anos
CowPI_file2$Subject  <- factor(CowPI_file2$Subject)
row.names(CowPI_file4) <- CowPI_file2$Subject
hclust_CowPi <- hclust(vegdist(CowPI_file4[,1:225],"mahalanobis"),"ward.D")
plot(hclust_CowPi)
# hclust_CowPi$order
# merge(hclust_CowPi$order,CowPI_file2,by.x="row.names",by.y="Subject")
###CLuster
require(ggdendro)
#We will color the labels according to countries(group_info[,1])
hc_d <- dendro_data(as.dendrogram(hclust_CowPi))
hc_d$labels$Type <- 
  merge(hc_d$labels,CowPI_file2,by.y="Subject",by.x="label")[,"Phase"]
hc_d$labels$Treatment <-  merge(hc_d$labels,CowPI_file2,by.y="Subject",by.x="label")[,"Treatment"]
hc_d$labels$SubjectID <-  merge(hc_d$labels,CowPI_file2,by.y="Subject",by.x="label")[,"sampleID"]
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
  ylab("Distance (Metabolic pathways: Mahalanobis)") + theme_bw()+
  theme(axis.text.y = element_text(color = hc_d$labels$Treatment),
        axis.title.y = element_blank())
p1 <- p1 + geom_text(data=label(hc_d),
                     aes(label=SubjectID, x=x, y=-0.8, colour=hc_d$labels$color))
p1


library(MASS)
require(caret)
CowPI_file6 <- CowPI_file4[,1:226]
CowPI_file6$FCR <- CowPI_file2$FCR
cowPI <- list(production = cbind(CowPI_file2[,c(1:10)],CowPI_file4$groups),
              CowPI_file6 = CowPI_file6)
saveRDS(cowPI,paste0(out_dir,"/","csv/","cowPI_R.RDS"))

if(!file.exists(paste0(out_dir,"/","csv/","cowPI_step_list.RDS"))){
# Set up repeated k-fold cross-validation
train.control <- trainControl(method = "cv", number = 10)
# Train the model
step.model <- train(FCR ~., data = CowPI_file6,
                    method = "leapSeq", 
                    tuneGrid = data.frame(nvmax = 1:35),
                    trControl = train.control
)
step.model$results
step.model$bestTune
summary(step.model$finalModel)

coef(step.model$finalModel,6)
plot(step.model$finalModel)
summary(step.model$finalModel)

pred <- predict(step.model, CowPI_file6)
acc
pred <- data.frame(pred = pred, FCR = CowPI_file6$FCR)

library(ggplot2)

ggplot(data = pred)+
  geom_line(aes(x = FCR, y = pred))+
  geom_point(aes(x = FCR, y = pred))
cor.test(pred$pred,pred$FCR)
RMSE(pred = pred$pred,obs = pred$FCR)
RMSE(pred = pred$pred,obs = pred$FCR)
coef(step.model$finalModel,5)


require(ggpubr)

pred$group <- CowPI_file4$groups

ggscatterhist_plot <-ggscatterhist(
  pred, x = "FCR", y = "pred",group="group",
  color = "group", size = 3, alpha = 0.6,
  palette = c("darkblue", "green", "orange","red"),
  #margin.params = list(fill = "group", color = "black", size = 0.2,
  margin.plot = "boxplot", centrality.para = "mean",
  marginal=T)

cowPI_step <- list(step.model= step.model,
                   pred = pred,
                   ggscatterhist_plot= ggscatterhist_plot)
saveRDS(cowPI_step,paste0(out_dir,"/","csv/","cowPI_step_list.RDS"))
  #+
#  geom_point(data = CowPI_file6, aes(x=FCR, y = dist))


} else {
  
  cowPI_step <- readRDS(paste0(out_dir,"/","csv/","cowPI_step_list.RDS"))
}

HE_groups <-CowPI_file2[which(CowPI_file3$groups=="High feed efficiency - solid" | 
                                CowPI_file3$groups=="High feed efficiency - liquid"),]
LE_groups <-CowPI_file2[which(CowPI_file3$groups=="Low feed efficiency - solid" | 
                                CowPI_file3$groups=="Low feed efficiency - liquid"),]
for(i in 1:nrow(HE_groups)){
  HE_groups[i,-227] <- HE_groups[i,-227]/sum(HE_groups[i,-227],na.rm = T)  
  
}
HE_groups <- HE_groups[,-227]
for(i in 1:nrow(LE_groups)){
  LE_groups[i,-227] <- LE_groups[i,-227]/sum(LE_groups[i,-227],na.rm = T)  
  
}
LE_groups <- LE_groups[,-227]


HE_groups <- HE_groups[,order(colSums(HE_groups),decreasing=TRUE)]
LE_groups <- LE_groups[,order(colSums(LE_groups),decreasing=TRUE)]
###########################################################################
  #Extract list of top N Taxa
  N<-10
  taxa_list<-colnames(HE_groups)[1:N]
  taxa_list <- taxa_list[which(!is.na(taxa_list))]
  #Generate a new table with everything added to Others
  #new_x<-data.frame(x[,colnames(x) %in% taxa_list],
  #Others=rowSums(x[,!colnames(x) %in% taxa_list]))
  
  new_HE_groups<-data.frame(HE_groups[,colnames(HE_groups) %in% taxa_list],
                    
                    Others= if(length(colnames(HE_groups))<N){Others=NA} else { rowSums(HE_groups[,!colnames(HE_groups) %in% taxa_list])}
  )
  
  new_HE_groups$Group <- "High feed efficiency"
  

  ##########
  taxa_list<-colnames(LE_groups)[1:N]
  taxa_list <- taxa_list[which(!is.na(taxa_list))]
  #Generate a new table with everything added to Others
  #new_x<-data.frame(x[,colnames(x) %in% taxa_list],
  #Others=rowSums(x[,!colnames(x) %in% taxa_list]))
  
  new_LE_groups<-data.frame(LE_groups[,colnames(LE_groups) %in% taxa_list],
                            
                            Others= if(length(colnames(LE_groups))<N){Others=NA} else { rowSums(LE_groups[,!colnames(LE_groups) %in% taxa_list])}
  )
  
  new_LE_groups$Group <- "Low feed efficiency"
  
  colMeans(new_HE_groups[,-12])
  colMeans(new_LE_groups[,-12])

  
  PCA <- FactoMineR::PCA(CowPI_file3[,-227],scale.unit = TRUE,ncp=5)#quali.sup = 251  ,ncp = 5)
  grp <-  factor(paste0(CowPI_file2$Treatment," - ",CowPI_file2$Phase),
                 levels =c("High feed efficiency - solid","High feed efficiency - liquid",
                           "Low feed efficiency - solid","Low feed efficiency - liquid"))
  
  fviz_pca_ind(PCA,
               label = "ind",# "ind", # hide individual labels
               habillage = grp,
               col.var = "black",#CowPI_file2$Treatment, # color by groups
               palette = c("darkblue", "green", "orange","red"),
               addEllipses = TRUE,
               ellipse.type = "convex",
               repel=T# Concentration ellipses
  )
  

  scatter3d(x = PCA$ind$coord[,1], y = PCA$ind$coord[,2], z = PCA$ind$coord[,3], 
            groups = grp ,
            xlab= paste("Axis 1",round(PCA$eig[1,2],2),"%"),
            ylab= paste("Axis 2",round(PCA$eig[2,2],2),"%"),
            zlab= paste("Axis 3",round(PCA$eig[3,2],2),"%"),
            surface=FALSE, grid = FALSE, ellipsoid = TRUE,
            labels = row.names(PCA$ind$coord),
            point.col =  c("black", "black", "black"),id.n=length(row.names(PCA$ind$coord)),
            surface.col = c("darkblue", "green", "orange","red"),
            axis.col = c("black", "black", "black")
  )
  require(car);require(rgl)
  legend3d("topright", legend = c("High feed efficiency - solid","High feed efficiency - liquid",
                                  "Low feed efficiency - solid","Low feed efficiency - liquid"),
           pch = 16, col = c("darkblue", "green", "orange","red"), cex=2, inset=c(0.01))
  
  saveRDS(CowPI_file2a,paste0(csv_dir,"/","CowPI_file2a",".RDS"))
  