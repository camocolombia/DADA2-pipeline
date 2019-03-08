##################################################
###Chrystian C. Sosa 2018 Chapter 1              #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##################################################
###CALLING LIBRARIES
require(ggplot2)
require(tidyr)
require(reshape2)
require(dplyr)
require(FactoMineR)
require(factoextra)
require(corrplot)
require(FactoMineR)
require(agricolae)
require(missMDA)
 #install.packages("missMDA")
###CALLING FUNCTIONS

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
###DEFINING PATHS
mainDir <- "E:/DADA2"
chapter <- "chapter_2"
dat_dir <- paste0(mainDir,"/",chapter,"/","data"); if(!file.exists(dat_dir)){dir.create(dat_dir)}
production_dir <- paste0(dat_dir,"/","production"); if(!file.exists(production_dir)){dir.create(production_dir)}
production <- read.csv(paste0(production_dir,"/","production_data.csv"),header=T,sep=",")
out_dir <- paste0(mainDir,"/",chapter,"/","outcomes"); if(!file.exists(out_dir)){dir.create(out_dir)}
graph_dir <- paste0(out_dir,"/","graphics"); if(!file.exists(graph_dir)){dir.create(graph_dir)}
csv_dir <- paste0(out_dir,"/","csv"); if(!file.exists(csv_dir)){dir.create(csv_dir)}



samples <- c(
  "640",
  "2121",
  "2126",
  "3112",
  "3802",
  "4791",
  "5940",
  "5941",
  "6060",
  "6069",
  "6071",
  "6072",
  "6088",
  "6091",
  "6092",
  "6100",
  "6105",
  "6113",
  "6466"
)


production <- production[which(production$sampleID %in% samples),]
production <- production[which(production$Phase %in% c("Liquid","Solid")),]
#http://www.sthda.com/english/wiki/normality-test-in-r
# dir <- "E:/Dropbox/Dropbox/Paper_PhD/chapter_1"
# 
# data_set <- read.csv(paste0(dir,"/","dataset.csv"),sep="|",header=T)
###REMOVING COMPLEX DATA
data_set <- production
#data_set <- data_set[complete.cases(data_set),];complete_data <- data_set
data_set <- production[!duplicated(production$sampleID),]
#data_set <- data_set[!duplicated(data_set$Sample_ID),]
#data_set <- data_set[,c(1:30)]

#x <- agricolae::kruskal(data_set$ADG,data_set$Breed)


vars <- c("sampleID","Breed","Phase","ADG_Treatment","ADG","FCR","LW.Gain","Ave.DMI","pH","LW","DW",
          "R_weight_full","R_weight_empty")
my_data <-data_set[,vars]
row.names(my_data) <- my_data[,1]
my_data <- my_data[, -1]
#plot(my_data [,-6],col=my_data[,6])
my_data$ADG_Treatment <- gsub("bottom","Low feed efficiency",my_data$ADG_Treatment); my_data$ADG_Treatment <- gsub("top","High feed efficiency",my_data$ADG_Treatment)
colnames(my_data) <- c("Breed","Rumen fraction","Treatment","ADG","FCR",
                       "Live weight gain","Average dry matter intake","pH","Live weight","Dead weight",
                       "Rumen weight full","Rumen weight empty")
# ++++++++++++++++++++++++++
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
#https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html


###PERFORMING CORRPLOT

pdf(file =paste0(graph_dir,"/","Corrplot.pdf"),width=90, height=80,onefile=T)

corrplot::corrplot(cor(my_data[,-c(1:3)]),method="color",col = col(200),hclust.method = "ward.D",
                   type = "upper", order = "hclust", number.cex = 10,tl.cex=20,cl.cex = 20,
                   addCoef.col = "black", # Add coefficient of correlation
                   tl.col = "black", tl.srt = 90, # Text label color and rotation
                   # Combine with significance
                   #p.mat = p.mat, sig.level = 0.01, insig = "blank", 
                   # hide correlation coefficient on the principal diagonal
                   diag = FALSE)        


### PERFORMING PCA 
#http://www.sthda.com/english/wiki/print.php?id=202


nb <- missMDA::estim_ncpPCA(my_data[,-c(1:3)],ncp.max=10) ## Time consuming, nb = 2
res.comp <- imputePCA(my_data[,-c(1:3)],ncp=as.numeric(nb[1]))

##
#res.pca <- PCA(my_data[,-c(1:3)], graph = FALSE,scale.unit = T)
res.pca <- PCA(res.comp$completeObs, graph = FALSE,scale.unit = T)  
# print(res.pca)
s <-fviz_screeplot(res.pca, ncp=5)
  
  
  a1 <- fviz_pca_biplot(res.pca, repel=T,
                  habillage = my_data$Breed,#factor(my_data[,6]),
                  palette=c("blue","red","green"),
                  addEllipses = TRUE,
                  ellipse.type = "convex",
                  col.var = "black",# alpha.var ="cos2",
                  label = c("ind","var"),
                  labelsize = 25,
                  pointsize=9, title = ""
                  ) +
#    scale_color_brewer(palette="Dark2")+
    
    theme_minimal()
  # a1 <- a1 +
  # theme(text=element_text(size=40),axis.text.x  = element_text(size=40,colour="black"),
  #       axis.text.y  = element_text(size=49,colour="black"))
  a1 <- a1 +theme(text=element_text(size=40),
                  legend.text=element_text(size=40),
                  axis.text.x  = element_text(size=40,colour="black"),
                  axis.text.y  = element_text(size=49,colour="black"))
  
  
  ggsave(paste0(graph_dir,"/","PCA",".pdf"),a1,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)
    # 
  res.desc <- dimdesc(res.pca, axes = c(1,2))

  res.desc
  
  #http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/117-hcpc-hierarchical-clustering-on-principal-components-essentials/
# res_HCPC <-  HCPC(res.pca, nb.clust = -1, max = NULL, graph = TRUE)
# 
# fviz_dend(res_HCPC, 
#           cex = 0.7,                     # Label size
#           palette = "jco",               # Color palette see ?ggpubr::ggpar
#           rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
#           rect_border = "jco",           # Rectangle color    
#           labels_track_height = 0.8      # Augment the room for labels
# )
# 
# 
# a2 <- fviz_pca_biplot(res.pca, 
#                 habillage = res_HCPC$data.clust$clust, addEllipses = TRUE,
#                 col.var = "red", alpha.var ="cos2",
#                 label = "var") +
#   scale_color_brewer(palette="Dark2")+
#   theme_minimal()
# 
# 
# a2
  
###PERFORMING SUMMARIZED TABLES
  
summary_list <- lapply(4:ncol(my_data),function(i){
  
  var_n <- colnames(my_data)[i]
 # var_to <- cbind(my_data[,c(i,ncol(my_data))])
  var_to <- cbind(my_data[,c(i,1)])
colnames(var_to) <- c("var","Breed")
  x <- as.data.frame(group_by(var_to, Breed) %>%
              summarise(
    count = n(),
    min = min(var, na.rm = TRUE),
    max = max(var, na.rm = TRUE),
    mean = mean(var, na.rm = TRUE),
    median = median(var, na.rm = TRUE),    
    sd = sd(var, na.rm = TRUE),
    cv = ((sd/ median) *100)
  ))


colnames(x) <- c("Breed","count","min","max","mean","median","sd","cv")
x$var <- NA; x$var <- rep(as.character(var_n),nrow(x))
x <- x[,c(ncol(x),1:(ncol(x)-1))]
return(x)
})


summary_list <- do.call(rbind,summary_list)

write.csv(summary_list,paste0(csv_dir,"/","summary_breed_production.csv"),quote = F)

###PERFORMING MANNWHITNEY - U - TEST

#summary_wilcox <- lapply(1:(ncol(my_data)-1),function(i){
summary_wilcox <- lapply(4:ncol(my_data),function(i){
  var_n <- colnames(my_data)[i]
  # var_to <- cbind(my_data[,c(i,ncol(my_data))])
  var_to <- cbind(my_data[,c(i,1)])
  colnames(var_to) <- c("var","Breed")

  length(unique(var_to[complete.cases(var_to),]$Breed))
  # w <-wilcox.test(var~Breed,data=var_to,paired=F,alternative="two.sided",exact=T,correct=T,mu=0,conf.int=T,conf.level=0.95)
  # 
  # x <- as.data.frame(matrix(ncol = 8,nrow=1))
  # x[,1] <- var_n
  # x[,2] <- w$statistic
  # x[,3] <- w$parameter
  # x[,4] <- w$p.value
  # x[,5] <- w$null.value
  # x[,6] <- w$alternative
  # x[,7] <- w$method
  # x[,8]<-NA; x[,8][which(x[,4] <= 0.05 )] <- "Significant";x[,8][which(x[,4] > 0.05 )] <- "Not significant"
  # 
  
  x <- as.data.frame(matrix(ncol = 10,nrow=1))
  x[,1] <- var_n
  
  if(length(unique(var_to[complete.cases(var_to),]$Breed))>2){
  w <-agricolae::kruskal(var_to$var,var_to$Breed, alpha = 0.05,group=F)

  

  x[,2] <- w$statistic[2]
  x[,3] <- w$statistic[3]
  x[,4] <- w$statistic[1]
  x[,5] <- w$comparison[which(row.names(w$comparison)=="Cheviot - Connemara"),][2]
  x[,6] <- w$comparison[which(row.names(w$comparison)=="Cheviot - Lanark"),][2]
  x[,7] <- w$comparison[which(row.names(w$comparison)=="Connemara - Lanark"),][2]
  x[,8]<-NA; x[,8][which(x[,3] <= 0.05 )] <- "Significant";x[,8][which(x[,3] > 0.05 )] <- "Not significant"
  x[,9]<-NA; x[,9] <- "Kruskal"
  x[,10]<-NA; x[,10] <- length(unique(var_to[complete.cases(var_to),]$Breed))

  } else {
    
    w <-wilcox.test(var~Breed,data=var_to,paired=F,alternative="two.sided",exact=T,correct=T,mu=0,conf.int=T,conf.level=0.95)
  
    x[,2] <- w$statistic[2]
    x[,3] <- w$p.value
    x[,4] <- w$statistic
    x[,5] <- "NO VALID, 2 BREEDS!"
    x[,6] <- "NO VALID, 2 BREEDS!"
    x[,7] <- "NO VALID, 2 BREEDS!"
    x[,8]<-NA; x[,8][which(x[,3] <= 0.05 )] <- "Significant";x[,8][which(x[,3] > 0.05 )] <- "Not significant"
    x[,9]<-NA; x[,9] <- "Wilcoxon rank sum test"
    x[,10]<-NA; x[,10] <- length(unique(var_to[complete.cases(var_to),]$Breed))
    }
  
  colnames(x) <- c("var","Df","p.value","Statistics",
                   "Cheviot - Connemara","Cheviot - Lanark","Connemara - Lanark",
                  "result","Method","n_Breeds")
  return(x)
  
})

summary_wilcox <- do.call(rbind,summary_wilcox)
#summary_wilcox$BF <- NA;summary_wilcox$BF <- p.adjust(p=summary_wilcox$p.value,method = "bonferroni",n=nrow(summary_wilcox))
#summary_wilcox$BH <- NA;summary_wilcox$BH <- p.adjust(p=summary_wilcox$p.value,method = "fdr",n=nrow(summary_wilcox))


  summary_wilcox <- summary_wilcox[order(summary_wilcox$p.value),]

  write.csv(summary_wilcox,paste0(csv_dir,"/","summary_wilcox_production.csv"),quote = F)
  
p_values <- ggplot(data=summary_wilcox,aes(x=reorder(var, p.value),y=p.value,fill="red"))+
  geom_bar(stat="identity",show.legend = F)+
  xlab("")+
  ylab("P value")+
 # ylim(0,0.1)+
  geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
  theme(text=element_text(size=60),
       legend.text=element_text(size=60),
       axis.text.x  = element_text(size=60,colour="black",angle = 90),
       axis.text.y  = element_text(size=60,colour="black"))+
scale_y_continuous(limits = c(0,1),breaks=seq(0,1,0.05))
p_values
ggsave(paste0(graph_dir,"/","P_value_graphics",".pdf"),p_values,dpi=300,width =80,height=90,units = "cm",scale=1.2,limitsize = FALSE)

#my_data$Treatment <- gsub("bottom","Bottom",my_data$Treatment); my_data$Treatment <- gsub("top","Top",my_data$Treatment)


###SCALED BOXPLOTS

labels <- c( "Breed","Rumen fraction","Treatment","Average daily gain  (kg/d)","Feed Conversion ratio",
"Live weight gain (kg)","Dry Mattter Intake (kg)","Rumen pH","Live weight (Kg)","Dead weight (Kg)",
"Rumen weight full (kg)","Rumen weight empty (Kg)")      


my_data_scaled <- as.data.frame(cbind(scale(my_data[,-c(1:4)],center = T,scale = T),my_data[,1]))
colnames(my_data_scaled)[ncol(my_data_scaled)] <- "Group"
my_data_scaled$Group <- as.character(my_data[,1])
my_data_scaled$ID <- NA;my_data_scaled$ID <- row.names(my_data_scaled)
#my_data_scaled$ID <- factor(my_data_scaled$ID)
row.names(my_data_scaled) <- 1:nrow(my_data_scaled)



data_s <- melt(my_data_scaled,id=c("ID","Group"))
data_s$value <- as.numeric(as.character(data_s$value))

p1bp <- ggplot(data=data_s,aes(x=variable,y=value,fill=Group))+
  geom_boxplot(position = "dodge",outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Scaled values")+ xlab("")+
  ggtitle("")+
  scale_fill_manual("legend", values = c("red","blue","green"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black"),
        legend.title=element_blank(),legend.position="none")

ggsave(paste0(graph_dir,"/","Production_BP_scaled",".pdf"),p1bp,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)
p1bp

###

#as.data.frame(cbind(scale(my_data[,-c(ncol(my_data))],center = T,scale = T),my_data[,ncol(my_data)]))


my_data2 <- my_data[,-3]
my_data2$ID <- NA; my_data2$ID <- row.names(my_data2)

data_unscaled <- melt(my_data2,id=c("ID","Breed","Rumen fraction"))
data_unscaled$value <- as.numeric(as.character(data_unscaled$value))



###Boxplots
# 
# p1bp_all <- ggplot(data=data_unscaled,aes(x=variable,y=value,fill=Treatment))+
#   geom_boxplot(position = "dodge",outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
#   stat_boxplot(geom ='errorbar') +
#   #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
#   #stat_summary(fun.y=mean, geom="point")+
#   #guides(fill=FALSE)+
#   #ylim(c(0,100))+
#   ylab("Scaled values")+ xlab("")+
#   ggtitle("")+
#   scale_fill_manual("legend", values = c("blue","red"))+
#   theme(panel.background = element_rect(fill = "gray95"),
#         text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
#         axis.text.y  = element_text(size=60,colour="black"),
#         legend.title=element_blank(),legend.position="none")
# 
# #ggsave(paste0(graph_dir,"/","Production_BP",".pdf"),p1bp,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)
# 
# p1bp_all <- p1bp_all + facet_wrap(. ~ Treatment, scales='free',ncol=2)
# p1bp_all
# ggsave(paste0(graph_dir,"/","Production_BP_ALL",".pdf"),p1bp_all,dpi=300,width =100,height=100,units = "cm",scale=1.2,limitsize = FALSE)
# 

summary_bp <- lapply(4:ncol(my_data),function(i){

  var_n <- colnames(my_data)[i]
  var_to <- cbind(my_data[,c(i,1,2)])
  colnames(var_to) <- c("var","Breed")
  var_to$Breed <- droplevels(var_to$Breed)

  p1a <- ggplot(data=var_to,aes(x=Breed,y=var))+
    geom_boxplot(data=var_to,aes(fill=Breed)) +
    stat_boxplot(geom ='errorbar') +
    stat_summary(fun.y=mean, geom="line", aes(group=1))  +
    stat_summary(fun.y=mean, geom="point")+
    guides(fill=FALSE)+
    #ylim(c(0,100))+
    ylab(labels[[i]])+ xlab("Group")+
    ggtitle("")+
    scale_fill_manual("legend", values = c("blue","red","green"))+
    theme(panel.background = element_rect(fill = "gray95"),
          text=element_text(size=70),axis.text.x  = element_text(size=70,colour="black",angle = 90, hjust = 1),
          axis.text.y  = element_text(size=70,colour="black"),
          legend.title=element_blank(),legend.position="none")

  ggsave(paste0(graph_dir,"/",colnames(my_data)[[i]],"_BP.pdf"),p1a,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)
  
  return(p1a)


})



summary_wilcox 