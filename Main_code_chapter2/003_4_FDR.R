
##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##FDR using Xiao et al., 2017                    #
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
# library("data.table")
# library("DESeq2")
# library("BSDA")
# library("adegenet")
#install.packages("devtools")
# library("devtools")
library("phytools")
library("massMap")
library("structSSI")
#library("microbiomeViz");
library("stringr")
library("StructFDR")
# library("gtools")
require("agricolae")
#install_github("JiyuanHu/massMap")

##################################################
#Defining paths to be analyzed
# mainDir <- "E:/Dropbox/Dropbox/Paper_PhD"
# chapter <- "chapter_1
# mainDir <- "/home/csosa"
# chapter <- "chapter_1"
mainDir <- "E:/DADA2"
chapter <- "chapter_2"
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
CowPI_dir <- paste0(dat_dir,"/","CowPI"); if(!file.exists(CowPI_dir)){dir.create(CowPI_dir)} # Copy here the production data


##################################################
##################################################
#Defining functions to be use in the FDR approaches
# calculate geometric means prior to estimate size factors

gm_mean = function(x, na.rm=TRUE){
  if(sum(x,na.rm = na.rm)!=0){
  z <-  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    } else {
   z <- 0
     }
  return(z)
  }


#Antilog

antilog<-function(lx,base) 
{ 
  lbx<-lx/log(exp(1),base=base) 
  result<-exp(lbx) 
  result 
} 

#Capitalize first letter

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#permutation function
perm.func <- function (X, Y) {
  return(list(X=X, Y=sample(Y)))
}

#Standar error about mean
sem_function <- function(x){
sd(x)/sqrt(length(x))
}
#T.TEST
#   test.func <- function (X, Y) {
#     obj <- apply(X, 1, function(x) {
#       ttest.obj <- t.test(x ~ Y)#wilcox.test(x ~ Y,paired=F)
#       c(ttest.obj$p.value, sign(ttest.obj$statistic))
#     })
#     return(list(p.value=obj[1, ], e.sign=obj[2, ]))
#   }


#Wilcoxon 

test.func <- function (X, Y) {  
  Y <- as.numeric(factor(Y))
  obj <- apply(X, 1, function(x) {                
    p.value <- suppressWarnings(wilcox.test(x ~ Y)$p.value)
    e.sign <- sign(mean(x[Y == 2]) - mean(x[Y == 1]))
    c(p.value, e.sign)          
  })
  return(list(p.value=obj[1, ], e.sign=obj[2, ])) 
}

##################################################
##################################################
#Loading Phyloseq object

# samples_to_use <- factor(c("3112",
#              "6105",
#              "6088",
#              "6060",
#              "3802",
#              "6071",
#              "6100",
#              "6072",
#              "5941",
#              "4791",
#              "2121",
#              "2126"
#   )
# )
ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))
ps <- subset_samples(ps, Phase %in% c("Solid","Liquid"))

ps@sam_data$ADG_Treatment <- factor(ps@sam_data$ADG_Treatment,levels = c("High feed efficiency","Low feed efficiency"))
ps3 <-  readRDS(paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))
ps3 <- subset_samples(ps3, Phase %in% c("Solid","Liquid"))

gc()
write.table(as.data.frame(t(ps@otu_table)),paste0(CowPI_dir,"/","otu_table.tab"),sep = "\t",row.names = T,quote = F)


# ps <- subset_samples(ps, sampleID %in% samples_to_use)
# ps3 <- subset_samples(ps2, sampleID %in% samples_to_use)
###DAPC

# #Crossvalidate DAPC
# xval <- xvalDapc(ps3@otu_table,grp=ps3@sam_data$Treatment,n.pca.max = 300, training.set = 0.9,
#                  result = "groupMean", center = TRUE, scale = FALSE,
#                  n.pca = NULL, n.rep = 30, xval.plot = TRUE)
# 
# 
# #set.seed(4)
# contrib <- adegenet::loadingplot(xval$DAPC$var.contr,
#                                  thres=.07, lab.jitter=1)
# 
# #Saving DAPC results
# dapc.results <- as.data.frame(xval$DAPC$posterior)
# dapc.results$pop <- ps3@sam_data$Treatment
# dapc.results$indNames <- rownames(dapc.results)
# dapc.results <- melt(dapc.results)
# colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
# 
# cols <- RColorBrewer::brewer.pal(n = length(unique(ps3@sam_data$Treatment)), name = "Dark2")
# 
# ###Composition plot
# p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
# p <- p + geom_bar(stat='identity') 
# p <- p + scale_fill_manual(values = cols) 
# p <- p + facet_grid(~Original_Pop, scales = "free")
# p <- p + theme(text=element_text(size=40),
#                legend.text=element_text(size=60),
#                axis.text.x = element_text(angle = 90, hjust = 1,size=60,colour="black"),
#                axis.text.y = element_text(size=49,colour="black"), 
#                strip.text.x = element_text(size = 80))+
#   guides(fill=guide_legend(title="Subpopulation"))
# 
# #p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = ))
# p
# 
# ggsave("E:/compoplot.pdf",p,dpi=600,width =280,height=100,units = "cm",scale=1,limitsize = FALSE)
# 
################################33
#Loading DESeq2 object
DESeq_list <- readRDS(paste0(csv_dir,"/","DESeq_list",".RDS"))

#Loading Agglomorative taxa objects
ps_object <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list",".RDS"))
ps_object2 <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list_trams",".RDS"))
ps_object3 <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list_melt",".RDS"))


# ps_object_names <- names(ps_object)
# ps_object_names2 <- names(ps_object2)
# ps_object_names3 <- names(ps_object3)
# ###############################################
# ps_object <- lapply(1:length(ps_object),function(i){
#   ps_object[[i]] <- subset_samples(ps_object[[i]], sampleID %in% samples_to_use)
#   
# })
# 
# names(ps_object) <- ps_object_names
# ###############################################
# ps_object2 <- lapply(1:length(ps_object2),function(i){
#   ps_object2[[i]] <- subset_samples(ps_object2[[i]], sampleID %in% samples_to_use)
#   
# })
# 
# names(ps_object2) <- ps_object_names2
###############################################
##################################################
##################################################
###Defining FDR approaches

###Benjamini- Hochberg






if(!file.exists(paste0(out_dir,"/","csv/","FDR_list.RDS"))){
#http://rstudio-pubs-static.s3.amazonaws.com/13988_bb11d85b79b2436280de434988558140.html
#Final count oer group and sample
B1 <- rowSums(as.data.frame(t(otu_table(subset_samples(ps, Breed=="Connemara")))))
B2 <- rowSums(as.data.frame(t(otu_table(subset_samples(ps, Breed=="Cheviot")))))
B3 <- rowSums(as.data.frame(t(otu_table(subset_samples(ps, Breed=="Lanark")))))

B1_ps3 <- (as.data.frame(t(otu_table(subset_samples(ps3, Breed=="Connemara")))))
B2_ps3 <- (as.data.frame(t(otu_table(subset_samples(ps3, Breed=="Cheviot")))))
B3_ps3 <- (as.data.frame(t(otu_table(subset_samples(ps3, Breed=="Lanark")))))

  
HE = rowSums(as.data.frame(t(otu_table(subset_samples(ps, ADG_Treatment=="High feed efficiency")))))
LE = rowSums(as.data.frame(t(otu_table(subset_samples(ps, ADG_Treatment=="Low feed efficiency")))))
#Subset per sample
HE_A = as.data.frame(t(otu_table(subset_samples(ps, ADG_Treatment=="High feed efficiency")))) 
LE_A = as.data.frame(t(otu_table(subset_samples(ps, ADG_Treatment=="Low feed efficiency"))))
#Final count per phase
SO_P = rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Solid")))))
LI_P = rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Liquid")))))
#EP_P = rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Epithelial")))))
#Subset per sample
SO_A_P = as.data.frame(t(otu_table(subset_samples(ps, Phase=="Solid")))) 
LI_A_P = as.data.frame(t(otu_table(subset_samples(ps, Phase=="Liquid"))))
#LI_A_P = as.data.frame(t(otu_table(subset_samples(ps, Phase=="Epithelial"))))
############################################################################################
#Final count oer group and sample transformed
HE_3 = rowSums(as.data.frame(t(otu_table(subset_samples(ps3, ADG_Treatment=="High feed efficiency")))))
LE_3 = rowSums(as.data.frame(t(otu_table(subset_samples(ps3, ADG_Treatment=="Low feed efficiency")))))
#Subset per sample
HE_A_3 = as.data.frame(t(otu_table(subset_samples(ps3, ADG_Treatment=="High feed efficiency")))) 
LE_A_3 = as.data.frame(t(otu_table(subset_samples(ps3, ADG_Treatment=="Low feed efficiency"))))
#Final count per phase
SO_P_3 = rowSums(as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Solid")))))
LI_P_3 = rowSums(as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Liquid")))))
#EP_P_3 = rowSums(as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Epithelial")))))
#Subset per sample
SO_A_P_3 = as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Solid")))) 
LI_A_P_3 = as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Liquid"))))
#EP_A_P_3 = as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Epithelial"))))

###SUBSET PER TREATMENT
SO_OTU <-otu_table(subset_samples(ps, Phase=="Solid"))
LI_OTU <-otu_table(subset_samples(ps, Phase=="Liquid"))
# EP_OTU <-otu_table(subset_samples(ps, Phase=="Epithelial"))
HE_OTU <-otu_table(subset_samples(ps, ADG_Treatment=="High feed efficiency"))
LE_OTU <-otu_table(subset_samples(ps, ADG_Treatment=="Low feed efficiency"))

B1_OTU <-otu_table(subset_samples(ps, Breed=="Cheviot"))
B2_OTU <-otu_table(subset_samples(ps, Breed=="Connemara"))
B3_OTU <-otu_table(subset_samples(ps, Breed=="Lanark"))

df <- data.frame(seqId=names(HE),
                 ###########################################
                 ###total counts
                 animal_per_sample=NA,
                 # animal_per_HE=NA,
                 # animal_per_LE=NA,
                 animal_per_B1=NA,
                 animal_per_B2=NA,
                 animal_per_B3=NA,
                 count_per_sample=NA,
                 # count_per_HE=NA,
                 # count_per_LE=NA,
                 count_per_B1=NA,
                 count_per_B2=NA,
                 count_per_B3=NA,
                 count_per_SO=NA,
                 count_per_LI=NA,
                 # count_per_EP=NA,
                sample_count= rowSums(t(otu_table(ps))>0),
                B1_gmean3 = as.numeric(lapply(1:nrow(B1_ps3),function(i){x2 <- gm_mean(B1_ps3[i,],na.rm = T);return(x2)})),
                B2_gmean3 = as.numeric(lapply(1:nrow(B2_ps3),function(i){x2 <- gm_mean(B2_ps3[i,],na.rm = T);return(x2)})),
                B3_gmean3 = as.numeric(lapply(1:nrow(B3_ps3),function(i){x2 <- gm_mean(B3_ps3[i,],na.rm = T);return(x2)})),
                
                B1_SD_gmean = as.numeric(lapply(1:nrow(B1_ps3),function(i){x2 <- sem_function(as.numeric(B1_ps3[i,]));return(x2)})),
                B2_SD_gmean = as.numeric(lapply(1:nrow(B1_ps3),function(i){x2 <- sem_function(as.numeric(B2_ps3[i,]));return(x2)})),
                B3_SD_gmean = as.numeric(lapply(1:nrow(B1_ps3),function(i){x2 <- sem_function(as.numeric(B3_ps3[i,]));return(x2)})),
                pvalue_Breed = NA,
                pvalue_phase = NA
                
                
                #  sample_count_HE = rowSums(as.data.frame(t(otu_table(subset_samples(ps,ADG_Treatment=="High feed efficiency"))))>0),
                #  sample_count_LE = rowSums(as.data.frame(t(otu_table(subset_samples(ps,ADG_Treatment=="Low feed efficiency"))))>0),
                #  
                # sample_count_LI = rowSums(as.data.frame(t(otu_table(subset_samples(ps,Phase=="Liquid"))))>0),
                #  sample_count_SO =rowSums(as.data.frame(t(otu_table(subset_samples(ps,Phase=="Solid"))))>0), 
                #  #sample_count_EP =rowSums(as.data.frame(t(otu_table(subset_samples(ps,Phase=="Epithelial"))))>0), 
                #  
                #  sample_count_HE_SO = rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Solid" & ADG_Treatment=="High feed efficiency"))))>0),
                #  sample_count_HE_LI = rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Liquid" & ADG_Treatment=="High feed efficiency"))))>0),
                # #sample_count_HE_EP = rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Epithelial" & ADG_ADG_Treatment=="High feed efficiency"))))>0),
                #  
                #  sample_count_LE_SO = rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Solid" & ADG_Treatment=="Low feed efficiency"))))>0),
                #  sample_count_LE_LI = rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Liquid" & ADG_Treatment=="Low feed efficiency"))))>0),  
                # #sample_count_LE_EP = rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Epithelial" & ADG_ADG_Treatment=="Low feed efficiency"))))>0),
                #  
                #  #a.greaterthan.b_T = HE > LE,
                #  SH =rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Solid" & ADG_Treatment=="High feed efficiency"))))),
                #  SL =rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Solid" & ADG_Treatment=="Low feed efficiency"))))),
                #  
                #  LH =rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Liquid" & ADG_Treatment=="High feed efficiency"))))),
                #  LL =rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Liquid" & ADG_Treatment=="Low feed efficiency"))))),
                # 
                #  # EH =rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Epithelial" & ADG_Treatment=="High feed efficiency"))))),
                #  # EL =rowSums(as.data.frame(t(otu_table(subset_samples(ps, Phase=="Epithelial" & ADG_Treatment=="Low feed efficiency"))))),
                #  # 
                # 
                # HE = HE,
                #  LE = LE,
                #  count_T = HE + LE,
                #  rel_HE = HE/(sum(HE,LE)),
                #  rel_LE = LE/(sum(HE,LE)),
                #  rel_total = (((HE+LE)/(sum(HE,LE))*100)),
                #  #HE_gmean = as.numeric(lapply(1:nrow(HE_A),function(i){x2 <- median(as.numeric((HE_A/sum(HE_A,LE_A))[i,]));return(x2)})),
                #  #LE_gmean = as.numeric(lapply(1:nrow(LE_A),function(i){x2 <- median(as.numeric((LE_A/sum(HE_A,LE_A))[i,]));return(x2)})),
                #  HE_gmean = as.numeric(lapply(1:nrow(HE_A),function(i){x2 <- mean(as.numeric((HE_A/sum(HE_A,LE_A))[i,]));return(x2)})),
                #  LE_gmean = as.numeric(lapply(1:nrow(LE_A),function(i){x2 <- mean(as.numeric((LE_A/sum(HE_A,LE_A))[i,]));return(x2)})),
                #  HE_SD_gmean = as.numeric(lapply(1:nrow(HE_A),function(i){x2 <- sd(as.numeric((HE_A/sum(HE_A,LE_A))[i,]));return(x2)})),
                #  LE_SD_gmean = as.numeric(lapply(1:nrow(LE_A),function(i){x2 <- sd(as.numeric((LE_A/sum(HE_A,LE_A))[i,]));return(x2)})),
                #  
                #  
                #  
                #  # HE_gmean = as.numeric(lapply(1:nrow(HE_A),function(i){x2 <- gm_mean((HE_A/sum(HE_A,LE_A))[i,],na.rm = T);return(x2)})),
                #  #LE_gmean = as.numeric(lapply(1:nrow(LE_A),function(i){x2 <- gm_mean((LE_A/sum(HE_A,LE_A))[i,],na.rm = T);return(x2)})),
                #  
                #  # HE_L = NA,
                #  # LE_L  =NA,
                #  change_Treat =NA,
                #  pvalue_T = NA,
                #  p.adj_T = NA,
                #  #a.greaterthan.b_Ph = SO_P > LI_P,
                #  S_Ph = SO_P,
                #  L_Ph = LI_P,
                #  # E_Ph = EP_P,
                # 
                #  count_Ph = SO_P + LI_P,
                #  rel_SO = SO_P/(sum(SO_P,LI_P)),
                #  rel_LI = LI_P/(sum(SO_P,LI_P)),
                #  rel_total_Ph = (((SO_P+LI_P)/(sum(SO_P,LI_P))*100)),
                #  #SO_gmean = as.numeric(lapply(1:nrow(SO_A_P),function(i){x2 <- gm_mean((SO_A_P/sum(SO_A_P,LI_A_P))[i,],na.rm = T);return(x2)})),
                #  #LI_gmean = as.numeric(lapply(1:nrow(LI_A_P),function(i){x2 <- gm_mean((LI_A_P/sum(SO_A_P,LI_A_P))[i,],na.rm = T);return(x2)})),
                #  SO_gmean = as.numeric(lapply(1:nrow(SO_A_P),function(i){x2 <- mean(as.numeric((SO_A_P/sum(SO_A_P,LI_A_P))[i,]));return(x2)})),
                #  LI_gmean = as.numeric(lapply(1:nrow(LI_A_P),function(i){x2 <- mean(as.numeric((LI_A_P/sum(SO_A_P,LI_A_P))[i,]));return(x2)})),
                # 
                #  SO_SD_gmean = as.numeric(lapply(1:nrow(SO_A_P),function(i){x2 <- sd(as.numeric((SO_A_P/sum(SO_A_P,LI_A_P))[i,]));return(x2)})),
                #  LI_SD_gmean = as.numeric(lapply(1:nrow(LI_A_P),function(i){x2 <- sd(as.numeric((LI_A_P/sum(SO_A_P,LI_A_P))[i,]));return(x2)})),
                #  
                #  #SO_gmean = as.numeric(lapply(1:nrow(SO_A_P),function(i){x2 <- median(as.numeric((SO_A_P/sum(SO_A_P,LI_A_P))[i,]));return(x2)})),
                #  #LI_gmean = as.numeric(lapply(1:nrow(LI_A_P),function(i){x2 <- median(as.numeric((LI_A_P/sum(SO_A_P,LI_A_P))[i,]));return(x2)})),
                #  
                #  # SO_P = NA,
                #  # LI_P  =NA,
                #  change_Phase =NA,
                #  pvalue_Phase = NA,
                #  p.adj_Phase = NA,
                #  ###########################################
                #  ###transformed counts
                #  #a.greaterthan.b_T3 = HE_3 > LE_3,
                #  SH3 =rowSums(as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Solid" & ADG_Treatment=="High feed efficiency"))))),
                #  SL3 =rowSums(as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Solid" & ADG_Treatment=="Low feed efficiency"))))),
                #  LH3 =rowSums(as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Liquid" & ADG_Treatment=="High feed efficiency"))))),
                #  LL3 =rowSums(as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Liquid" & ADG_Treatment=="Low feed efficiency"))))),
                # #EH3 =rowSums(as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Epithelial" & ADG_Treatment=="High feed efficiency"))))),
                # #EL3 =rowSums(as.data.frame(t(otu_table(subset_samples(ps3, Phase=="Epithelial" & ADG_Treatment=="Low feed efficiency"))))),
                # 
                #  HE3 = HE_3,
                #  LE3 = LE_3,
                #  count_T3 = HE_3 + LE_3,
                # rel_total3 = (((HE_3+LE_3)/(sum(HE_3,LE_3))*100)),
                # #HE_gmean3 = as.numeric(lapply(1:nrow(HE_A_3),function(i){x2 <- gm_mean(HE_A_3[i,],na.rm = T);return(x2)})),
                # #LE_gmean3 = as.numeric(lapply(1:nrow(LE_A_3),function(i){x2 <- gm_mean(LE_A_3[i,],na.rm = T);return(x2)})),
                # HE_gmean3 = as.numeric(lapply(1:nrow(HE_A_3),function(i){x2 <- mean(as.numeric(HE_A_3[i,]),na.rm = T);return(x2)})),
                # LE_gmean3 = as.numeric(lapply(1:nrow(LE_A_3),function(i){x2 <- mean(as.numeric(LE_A_3[i,]),na.rm = T);return(x2)})),
                # HE_SD_gmean3 = as.numeric(lapply(1:nrow(HE_A_3),function(i){x2 <- sd(as.numeric(HE_A_3[i,]),na.rm = T);return(x2)})),
                # LE_SD_gmean3 = as.numeric(lapply(1:nrow(LE_A_3),function(i){x2 <- sd(as.numeric(LE_A_3[i,]),na.rm = T);return(x2)})),
                # 
                # 
                # change_Treat3 =NA,
                #  pvalue_T3 = NA,
                #  p.adj_T3 = NA,
                #  #a.greaterthan.b_Ph3 = SO_P_3 > LI_P_3,
                #  S_Ph3 = SO_P_3,
                #  L_Ph3 = LI_P_3,
                #  count_Ph3 = SO_P_3 + LI_P_3,
                #  rel_total_Ph3 = (((SO_P_3+LI_P_3)/(sum(SO_P_3,LI_P_3))*100)),
                # # SO_gmean3 = as.numeric(lapply(1:nrow(SO_A_P_3),function(i){x2 <- gm_mean(SO_A_P_3[i,],na.rm = T);return(x2)})),
                # # LI_gmean3 = as.numeric(lapply(1:nrow(LI_A_P_3),function(i){x2 <- gm_mean(LI_A_P_3[i,],na.rm = T);return(x2)})),
                #  SO_gmean3 = as.numeric(lapply(1:nrow(SO_A_P_3),function(i){x2 <- mean(as.numeric(SO_A_P_3[i,]),na.rm = T);return(x2)})),
                #  LI_gmean3 = as.numeric(lapply(1:nrow(LI_A_P_3),function(i){x2 <- mean(as.numeric(LI_A_P_3[i,]),na.rm = T);return(x2)})),
                # #EP_gmean3 = as.numeric(lapply(1:nrow(EP_A_P_3),function(i){x2 <- mean(as.numeric(EP_A_P_3[i,]),na.rm = T);return(x2)})),
                # 
                #  SO_SD_gmean3 = as.numeric(lapply(1:nrow(SO_A_P_3),function(i){x2 <- sd(as.numeric(SO_A_P_3[i,]),na.rm = T);return(x2)})),
                #  LI_SD_gmean3 = as.numeric(lapply(1:nrow(LI_A_P_3),function(i){x2 <- sd(as.numeric(LI_A_P_3[i,]),na.rm = T);return(x2)})),
                # #EP_SD_gmean3 = as.numeric(lapply(1:nrow(EP_A_P_3),function(i){x2 <- sd(as.numeric(EP_A_P_3[i,]),na.rm = T);return(x2)})),
                # 
                # change_Phase3 =NA,
                # pvalue_Phase3 = NA,
                # p.adj_Phase3 = NA,
                # 
                # rel_ch3_FE = NA,
                # rel_ch3_Phase = NA,
                # # rel_ch3_HE_PHASE = NA,
                # # rel_ch3_LE_PHASE = NA,
                # # rel_ch3_LIQUID= NA,
                # # rel_ch3_SOLID = NA,
                # 
                # status3_FE = NA,
                # status3_Phase = NA
                # # status3_HE_PHASE = NA,
                # # status3_LE_PHASE = NA,
                # # status3_LIQUID= NA,
                # # status3_SOLID = NA
                # 
                # 
                # 
                
)


# 
# df$HE_L <- log2(df$HE_gmean);df$LE_L <- log2(df$LE_gmean);
# df$SO_P <- log2(df$SO_gmean);df$LI_P <- log2(df$LI_gmean);
#df$HE_L[which(is.infinite(df$HE_L))] <- 0;df$LE_L[which(is.infinite(df$LE_L))] <- 0;
#http://rstudio-pubs-static.s3.amazonaws.com/13988_bb11d85b79b2436280de434988558140.html
#df$change_Treat <- df$HE_L - df$LE_L; df$change_Phase <- df$SO_P - df$LI_P;
#df$change_Treat <- log2(df$HE_gmean/df$LE_gmean); df$change_Phase <- log2(df$SO_gmean/df$LI_gmean);
#########################################
#Calculating log2 fold change

# df$change_Treat <- log2((df$HE_gmean/df$LE_gmean)); df$change_Phase <- log2((df$SO_gmean/df$LI_gmean));
# df$change_Treat3 <- log2((df$HE_gmean3/df$LE_gmean3)); df$change_Phase3 <- log2((df$SO_gmean3/df$LI_gmean3));


############counts of animals
val_dat <- sample_data(ps)
animalsID <- unique(sample_data(ps)$sampleID)

#for(i in 1:nrow(df)){
  
  animals_count <- list()
  animals_HE_count <- list()  
  animals_LE_count <- list()
  
  animals_B1_count <- list()
  animals_B2_count <- list()
  animals_B3_count <- list()
    for(k in 1:length(animalsID)){
  cat(k,"\n")

    
animals_count[[k]] <- suppressMessages(rowSums(as.data.frame(t(otu_table(subset_samples(ps,sampleID==animalsID[[k]]))))>0))
animals_count[[k]][which(animals_count[[k]]>0)]<- 1


if(as.character(unique(val_dat[which(val_dat$sampleID==animalsID[[k]]),]$Treatment))== "High feed efficiency"){
  animals_HE_count[[k]] <- suppressMessages(rowSums(as.data.frame(t(otu_table(subset_samples(ps,sampleID==animalsID[[k]] &
  Treatment=="High feed efficiency"))))>0))
  animals_HE_count[[k]][which(animals_HE_count[[k]]>0)]<- 1
  
  animals_LE_count[[k]] <- animals_count[[k]]
  animals_LE_count[[k]][which(animals_count[[k]]>0)]<- 0
  
  
} else if(as.character(unique(val_dat[which(val_dat$sampleID==animalsID[[k]]),]$Treatment))== "Low feed efficiency"){
  animals_HE_count[[k]] <- animals_count[[k]]
  animals_HE_count[[k]][which(animals_HE_count[[k]]>0)]<- 0
  animals_LE_count[[k]] <- suppressMessages(rowSums(as.data.frame(t(otu_table(subset_samples(ps,sampleID==animalsID[[k]] &
  Treatment=="Low feed efficiency"))))>0))
  animals_LE_count[[k]][which(animals_LE_count[[k]]>0)]<- 1
}                                                                            
#       animal_per_sample=NA
# animal_per_HE=NA
# animal_per_LE=NA
############################################ ############################################
############################################ ############################################  
############################################ ############################################ 
############################################ ############################################ 

if(as.character(unique(val_dat[which(val_dat$sampleID==animalsID[[k]]),]$Breed))== "Connemara"){
  animals_B1_count[[k]] <- animals_count[[k]]
  #animals_B1_count[[k]][which(animals_B1_count[[k]]>0)]<- 0
  animals_B1_count[[k]] <- suppressMessages(rowSums(as.data.frame(t(otu_table(subset_samples(ps,sampleID==animalsID[[k]] &
                                                                                               Breed=="Connemara"))))>0))
  animals_B1_count[[k]][which(animals_B1_count[[k]]>0)]<- 1
} else if(as.character(unique(val_dat[which(val_dat$sampleID==animalsID[[k]]),]$Breed))== "Cheviot"){
  animals_B2_count[[k]] <- animals_count[[k]]
  #animals_B1_count[[k]][which(animals_B1_count[[k]]>0)]<- 0
  animals_B2_count[[k]] <- suppressMessages(rowSums(as.data.frame(t(otu_table(subset_samples(ps,sampleID==animalsID[[k]] &
                                                                                               Breed=="Cheviot"))))>0))
  animals_B2_count[[k]][which(animals_B2_count[[k]]>0)]<- 1
  
  
}else if(as.character(unique(val_dat[which(val_dat$sampleID==animalsID[[k]]),]$Breed))== "Lanark"){
  animals_B3_count[[k]] <- animals_count[[k]]
  #animals_B1_count[[k]][which(animals_B1_count[[k]]>0)]<- 0
  animals_B3_count[[k]] <- suppressMessages(rowSums(as.data.frame(t(otu_table(subset_samples(ps,sampleID==animalsID[[k]] &
                                                                                               Breed=="Lanark"))))>0))
  animals_B3_count[[k]][which(animals_B3_count[[k]]>0)]<- 1
  
  
}
};rm(k)

#    };rm(i)
  
  
  
  df$animal_per_sample <- rowSums(do.call(cbind,animals_count))
  df$animal_per_HE <- rowSums(do.call(cbind,animals_HE_count))
  df$animal_per_LE <- rowSums(do.call(cbind,animals_LE_count))
  df$animal_per_B1 <- rowSums(do.call(cbind,animals_B1_count))
  df$animal_per_B2 <- rowSums(do.call(cbind,animals_B2_count))
  df$animal_per_B3 <- rowSums(do.call(cbind,animals_B3_count))
  
  
for(i in 1:nrow(df)){
cat(i,"\n")
  df$count_per_sample[i] <-sum(otu_table(ps)[,i]>0)
  df$count_per_HE[i] <-sum(HE_OTU[,i]>0)
  df$count_per_LE[i] <-sum(LE_OTU[,i]>0)
  df$count_per_SO[i] <-sum(SO_OTU[,i]>0)
  df$count_per_LI[i] <-sum(LI_OTU[,i]>0)
  df$count_per_B1[i] <-sum(B1_OTU[,i]>0)
  df$count_per_B2[i] <-sum(B2_OTU[,i]>0)
  df$count_per_B3[i] <-sum(B3_OTU[,i]>0)
  
  #df$count_per_EP[i] <-sum(Ep_OTU[,i]>0)
  
  };rm(i)

# log2HE_gmean <- log2(df$HE_gmean); log2LE_gmean <- log2(df$LE_gmean);
# log2HE_gmean[which(log2HE_gmean==-Inf )]  <-0; log2LE_gmean[which(log2LE_gmean==-Inf )]  <-0;
# 
# log2HE_gmean3 <- log2(df$HE_gmean3); log2LE_gmean3 <- log2(df$LE_gmean3);
# log2HE_gmean3[which(log2HE_gmean3==-Inf )]  <-0; log2LE_gmean3[which(log2LE_gmean3==-Inf )]  <-0;
# 
# log2SO_gmean <- log2(df$SO_gmean); log2LI_gmean <- log2(df$LI_gmean);
# log2SO_gmean[which(log2SO_gmean==-Inf )]  <-0; log2LI_gmean[which(log2LI_gmean==-Inf )]  <-0;
# 
# log2SO_gmean3 <- log2(df$SO_gmean3); log2LI_gmean3 <- log2(df$LI_gmean3);
# log2SO_gmean3[which(log2SO_gmean3==-Inf )]  <-0; log2LI_gmean3[which(log2LI_gmean3==-Inf )]  <-0;
# 
# 
# df$change_Treat <-log2HE_gmean - log2LE_gmean; df$change_Phase <-log2SO_gmean3 - log2LI_gmean3;
# df$change_Treat3 <-log2HE_gmean3 - log2LE_gmean3; df$change_Phase3 <-log2SO_gmean3 - log2LI_gmean3;
# 
# 
# for(i in 1:nrow(df)){
#   if(df$LE_gmean[i]==0){
#     df$change_Treat[i] <-  abs(df$change_Treat[i])
#   } else if(df$HE_gmean[i]==0) {
#     df$change_Treat[i] <-  - abs(df$change_Treat[i])
#     }else {
#     df$change_Treat[i] <- df$change_Treat[i]
#   }
# 
# };rm(i)
# 
# for(i in 1:nrow(df)){
#   if(df$LI_gmean[i]==0){
#     df$change_Phase[i] <-  abs(df$change_Phase[i])
#   } else if(df$SO_gmean[i]==0) {
#     df$change_Phase[i] <- abs(df$change_Phase[i])
#   }else {
#     df$change_Phase[i] <- df$change_Phase[i]
#   }
# };rm(i)
# #####################
# for(i in 1:nrow(df)){
#   if(df$LE_gmean3[i]==0){
#     df$change_Treat3[i] <-  abs(df$change_Treat3[i])
#   } else if(df$HE_gmean3[i]==0) {
#     df$change_Treat3[i] <- - abs(df$change_Treat3[i])
#   }else {
#     df$change_Treat3[i] <- df$change_Treat3[i]
#   }
# };rm(i)
# 
# for(i in 1:nrow(df)){
#   if(df$LI_gmean3[i]==0){
#     df$change_Phase3[i] <-  abs(df$change_Phase3[i])
#   } else if(df$SO_gmean3[i]==0) {
#     df$change_Phase3[i] <- - abs(df$change_Phase3[i])
#   }else {
#     df$change_Phase3[i] <- df$change_Phase3[i]
#   }
# };rm(i)
# 
# ###RELATIVE CHANGE
# df$rel_ch3_FE = (df$HE_gmean3 - df$LE_gmean3)/df$HE_gmean3
# df$rel_ch3_Phase = (df$LI_gmean3 - df$SO_gmean3)/df$LI_gmean3
# 
# df$status3_FE = df$HE_gmean3 > df$LE_gmean3
# df$status3_Phase = df$LI_gmean3 > df$SO_gmean3
# 
# 
# 
# 
# ####
# for(j in 1:nrow(df)){
#   cat(j,"\n")
#   #cat(j,"\n")
#   df$status3_FE[[j]] = if(df$status3_FE[[j]]==T){df$status3_FE[[j]]="HE"}else if(df$HE_gmean3 == df$LE_gmean3){
#     df$status3_FE[[j]]="EQUAL"} else {df$status3_FE[[j]]="LE"}
#   
#   df$status3_Phase[[j]] = if(df$status3_Phase[[j]]==T){df$status3_Phase[[j]]="LI"} else if( df$LI_gmean3 == df$SO_gmean3){
#     df$status3_Phase[[j]]="EQUAL"} else {df$status3_Phase[[j]]="SO"}
#   
#  
# };rm(j)

#####################3
# 
# status3_FE = NA,
# status3_Phase = NA

# for(i in 1:nrow(df)){
#   if(is.infinite(df$change_Treat[i]) | is.na(df$change_Treat[i])){
#     if(df$HE[i]>0){
#       df$change_Treat[i] <- log2(df$HE[i])
#     } else { 
#          df$change_Treat[i] <- -(abs(log2(df$LE[i])))}
#   } else {df$change_Treat[i] <- df$change_Treat[i] }
# };rm(i)
# 
# for(i in 1:nrow(df)){
#   if(is.infinite(df$change_Phase[i]) | is.na(df$change_Treat[i])){
#     if(df$S_Ph[i]>0){
#       df$change_Phase[i] <- log2(df$S_Ph[i])
#     } else { df$change_Phase[i] <- - abs(log2(df$L_Ph[i]))}
#   } else {df$change_Phase[i] <- df$change_Phase[i] }
# };rm(i)
# 
# ##prop count
# for(i in 1:nrow(df)){
#   if(is.infinite(df$change_Treat3[i]) | is.na(df$change_Treat3[i])){
#     if(df$HE3[i]>0){
#       df$change_Treat3[i] <- log2(df$HE3[i])
#     } else { df$change_Treat3[i] <- - abs(log2(df$LE3[i]))}
#   } else {df$change_Treat3[i] <- df$change_Treat3[i] }
# };rm(i)
# 
# for(i in 1:nrow(df)){
#   if(is.infinite(df$change_Phase3[i]) | is.na(df$change_Phase3[i])){
#     if(df$S_Ph3[i]>0){
#       df$change_Phase3[i] <- log2(df$S_Ph3[i])
#     } else { df$change_Phase3[i] <- - abs(log2(df$L_Ph3[i]))}
#   } else {df$change_Phase3[i] <- df$change_Phase3[i]}
# };rm(i)
#############################
#Calculating Wilcoxon per ASV
for(i in 1:nrow(df)){ 
  cat(i,"\n")
  x <- data.frame(var=as.numeric(c(B1_ps3[i,],B2_ps3[i,],B3_ps3[i,])),
                  Breed=c(rep("Connemara",ncol(B1_ps3)),rep("Cheviot",ncol(B2_ps3)),rep("Lanark",ncol(B3_ps3))
    ))
  w <-agricolae::kruskal( x$var, x$Breed, alpha = 0.05,group=F)
  

  df$pvalue_Breed[[i]] <- as.numeric(w$statistics[3])
  df$pvalue_phase[[i]] <-  t.test(x=LI_A_P_3[i,],y=SO_A_P_3[i,],alternative = "two.sided",paired = F)$p.value
  # # df$pvalue[[i]] <- t.test(x=HE_A[i,],y=LE_A[i,],alternative = "two.sided",paired = F)$p.value
  # if((df$count_T[[i]])>0){
  #   #x <- as.numeric(log2(HE_A[i,])); x[which(is.infinite(x))]<-0
  #   #y <- as.numeric(log2(LE_A[i,])); y[which(is.infinite(y))]<-0
  #   #Relative abundances change Treatment
  #   x <- as.numeric(HE_A[i,]);y <- as.numeric(LE_A[i,]);
  #   x <- x/sum(x,y);y <- y/sum(x,y)
  #   #x <- x/sum(x);y <- y/sum(y)
  #   if(length(x[which(is.na(x))])==length(x)){ x <- rep(0,length(x)) }
  #   if(length(y[which(is.na(y))])==length(y)){ y <- rep(0,length(y)) }
  #   #Relative abundances change Phase
  #   x2 <- as.numeric(SO_A_P[i,]);y2 <- as.numeric(LI_A_P[i,]);
  #   x2 <- x2/sum(x2,y2);y2 <- y2/sum(x2,y2)
  #   #x2 <- x2/sum(x2);y2 <- y2/sum(y2)
  #   
  #   if(length(x[which(is.na(x2))])==length(x2)){ x2 <- rep(0,length(x2)) }
  #   if(length(y2[which(is.na(y2))])==length(y2)){ y2 <- rep(0,length(y2)) }
  #   
  #   #     df$pvalue[[i]] <- z.test(x = x, y = y, # Two samples with normal distribution
  #   #            alt = "two.sided",          # Dos colas
  #   #           mu = 0,                     # H_0: mu_1 - mu_2 = 0
  #   #            sigma.x = sd(x),     # desviación estándar m
  #   #            sigma.y = sd(y),     # desviación estandar n
  #   #            conf.level = 0.95)$p.value
  #   # #    df$pvalue[[i]] <- t.test(x=x,y=y,alternative = "two.sided",paired = F)$p.value
  #   df$pvalue_T[[i]] <-  wilcox.test(x=x,y=y,paired=F,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value
  #   df$pvalue_Phase[[i]] <-  wilcox.test(x=x2,y=y2,paired=T,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value
  #   #Prop
  #   df$pvalue_T3[[i]] <-  wilcox.test(x=as.numeric(HE_A_3[i,]),y=as.numeric(LE_A_3[i,]),paired=F,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value
  #   df$pvalue_Phase3[[i]] <-  wilcox.test(x=as.numeric(SO_A_P_3[i,]),y=as.numeric(LI_A_P_3[i,]),paired=T,alternative="two.sided",exact=T,correct=F,conf.int=T,conf.level=0.95)$p.value
  #   
  # } else {
  #   df$pvalue_T[[i]] <- NA 
  #   df$pvalue_Phase[[i]] <- NA 
  #   df$pvalue_T3[[i]] <- NA 
  #   df$pvalue_Phase3[[i]] <- NA 
  # } 
};rm(i)
#ggplot(df, aes(fc, log2fc, colour = a.greaterthan.b), size = 8) + geom_point()
#df <- df[which(!is.na(df$pvalue)),]
# df$p.adj_T <- p.adjust(p = df$pvalue_T, method =  "BH",n = length(df$pvalue_T))
# df$p.adj_Phase <- p.adjust(p = df$pvalue_Phase, method =  "BH",n = length(df$pvalue_Phase))
#prop
df$p.adj_Breed <- p.adjust(p = df$pvalue_Breed, method =  "BH")#,n = length(df$pvalue_Breed))
df$p.adj_Phase3 <- p.adjust(p = df$pvalue_phase, method =  "BH")#,n = length(df$pvalue_phase))

#df <- df[complete.cases(df),]
tax_table_df <- as.data.frame(ps@tax_table)
tax_table_df$seqId <- NA;
tax_table_df$seqId <- as.character(row.names(tax_table_df))
df$seqId <- as.character(df$seqId)
df <- dplyr::full_join(df,tax_table_df,by="seqId")
#df <- merge(df,ps@tax_table,by="row.names",all=F,sort=F)
df_complete <- df

df_filtered <- df_complete[which(df_complete$p.adj_Breed <= 0.05),]

write.csv(df_complete,paste0(csv_dir,"/","df_logchange.csv"),row.names = F,quote = F)
write.csv(df_filtered,paste0(csv_dir,"/","df_logchange_filtered.csv"),row.names = F,quote = F)

# df_log2fold <- df
# df_log2fold <- df_log2fold[which(df_log2fold$change_Treat3 !=0),]
#f
# length(subset(df_log2fold$change_Treat3,df_log2fold$change_Treat3>0))/length(df_log2fold$change_Treat3)
# length(subset(df_log2fold$change_Treat3,df_log2fold$change_Treat3<0))/length(df_log2fold$change_Treat3)

################################################################
#Relative percentage summary


# #count
# comp_list <- list(
#   # SH = (df$SH/(sum(df$SH,df$SL))),
#   # SL = (df$SL/(sum(df$SH,df$SL))),
#   # LH = (df$LH/(sum(df$LH,df$LL))),
#   # LL = (df$LL/(sum(df$LH,df$LL)))
#   SH = (df$SH/(sum(df$SH,df$LH))),
#   LH = (df$LH/(sum(df$SH,df$LH))),
#   SL = (df$SL/(sum(df$SL,df$LL))),
#   LL = (df$LL/(sum(df$SL,df$LL)))
# )
# comp_list  <- melt(comp_list)
# comp_list$L1 <- factor(comp_list$L1,levels = c("LL","LH","SL","SH"))
# comp_list$Treatment <- NA;
# comp_list$Treatment[which(comp_list$L1=="LL" | comp_list$L1=="SL")] <- "Low feed efficiency"
# comp_list$Treatment[which(comp_list$L1=="LH" | comp_list$L1=="SH")] <- "High feed efficiency"
# 
# comp_list2 <-  list(
#   # SH = (df$SH/(sum(df$SH,df$SL))),
#   # SL = (df$SL/(sum(df$SH,df$SL))),
#   # LH = (df$LH/(sum(df$LH,df$LL))),
#   # LL = (df$LL/(sum(df$LH,df$LL)))
#   SH = df$SH3,
#   LH = df$LH3,
#   SL = df$SL3,
#   LL = df$LL3
# # )
# comp_list2  <- melt(comp_list2)
# comp_list2$L1 <- factor(comp_list2$L1,levels = c("LL","LH","SL","SH"))
# comp_list2$Treatment <- NA;
# comp_list2$Treatment[which(comp_list2$L1=="LL" | comp_list2$L1=="SL")] <- "Low feed efficiency"
# comp_list2$Treatment[which(comp_list2$L1=="LH" | comp_list2$L1=="SH")] <- "High feed efficiency"
# 
#   
# comp_list2 <- comp_list2[which(comp_list2$value>0),]
# 
# ##proportional  relative abundance transformation
# 
# #ggplot(comp_list, aes(value, fill = L1,colour = L1)) +
# comp_list_plot <- ggplot(comp_list2, aes(x= L1, y=value, fill = Treatment)) +
# geom_boxplot(na.rm=T)+#+
# ylim(c(0,0.0001))
  #  geom_freqpoly(bins=400)+
  #scale_x_continuous(limits = c(0,0.0001),breaks=seq(0,0.0001,0.00001))

#kruskal.test(L1 ~ value,data=comp_list)
#kruskal(comp_list$value, comp_list$L1, console = TRUE)
# out<-with(comp_list,kruskal(y=value,trt=L1,main="comp_list",group=F, p.adj="BH"))
# 
# #print(out$comparison)
# #prop
# out2<-with(comp_list2,kruskal(y=value,trt=L1,main="comp_list2",group=F, p.adj="BH"))
# 
# 
# rm(HE,LE,HE_A,LE_A,SO_P,LI_P,SO_A_P,LI_A_P,HE_3,LE_3,HE_A_3,LE_A_3,SO_P_3,LI_P_3,SO_A_P_3,LI_A_P_3,x,x2,y,y2);gc()
# #print(out2$comparison)

########################################################################################################################
# #Wilcoxon summary test
# df_pvalues <- data.frame(
#   # High_vs_Low = wilcox.test(x=(df$HE/(sum(df$HE,df$LE))),y=(df$LE/sum(df$HE,df$LE)),paired=F,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
#   # Solid_vs_liquid = wilcox.test(x=(df$S_Ph/(sum(df$S_Ph,df$L_Ph))),y=(df$L_Ph/sum(df$S_Ph,df$L_Ph)),paired=T,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
#   # Solids = wilcox.test(x=(df$SH/(sum(df$SH,df$SL))),y=(df$SL/sum(df$SH,df$SL)),paired=F,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
#   # Liquids = wilcox.test(x=(df$LH/(sum(df$LH,df$LL))),y=(df$LL/sum(df$LH,df$LL)),paired=F,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
#   # High_S_vs_L = wilcox.test(x=(df$SH/(sum(df$SH,df$LH))),y=(df$LH/sum(df$SH,df$LH)),paired=T,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
#   # Low_S_vs_L =wilcox.test(x=(df$SL/(sum(df$SL,df$LL))),y=(df$LL/sum(df$SL,df$LL)),paired=T,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
# #prop
# High_vs_Low3 = wilcox.test(x=df$HE3,y=df$LE3,paired=F,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
# Solid_vs_liquid3 = wilcox.test(x=df$S_Ph3,y=df$L_Ph3,paired=T,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
# Solids3 = wilcox.test(x=df$SH3,y=df$SL3,paired=F,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
# Liquids3 = wilcox.test(x=df$LH3,y=df$LL3,paired=F,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
# High_S_vs_L3 = wilcox.test(x=df$SH3,y=df$LH3,paired=T,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value,
# Low_S_vs_L3 =wilcox.test(x=df$SL3,y=df$LL3,paired=T,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value
# 
#   
#   )

# df_pvalues[2,] <- p.adjust(p = df_pvalues[1,] , method =  "BH",n = ncol(df_pvalues))
# df_pvalues <- as.data.frame(t(df_pvalues))
# df_pvalues$sig <- NA;
# df_pvalues$sig[which(df_pvalues[,2] <= 0.05)] <-"**"; df_pvalues$sig[which(df_pvalues[,2] > 0.05)] <-"-"; 
#colnames(df_pvalues) <- c("p.value","p.adj","sig")
#write.csv(df_pvalues,paste0(csv_dir,"/","wilcoxon_pvalues.csv"),row.names = F,quote = F)
#write.csv(df_log2fold,paste0(csv_dir,"/","df_logchange_no0.csv"),row.names = F,quote = F)


#sum(comp_list2$value)
########################################
# #Log2foldchange fixed per change in HE vs LE plot
# df_p <- df
# df_p$rel_HE <- df_p$rel_HE*100;df_p$rel_LE <- df_p$rel_LE*100
# df_p <- df_p[which(df_p$p.adj_T3 <= 0.05),]
# 
# df_p <- df_p[which(df_p$change_Treat3 > 0.5 | df_p$change_Treat3 < 0.5),]
# df_p <- df_p[which(!is.na(df_p$Genus)),]
# df_p$anti <- antilog(df_p$change_Treat3,2)
# #Log2 fold change plot
# exp_plot2 <- ggplot(df_p, aes(y=change_Treat3, x=reorder(Genus, -change_Treat3), color=Phylum, fill=Phylum)) + 
#   geom_point(size=25) + 
#     #geom_bar(stat="identity")+
#   geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
#   ylab("Log2foldchange")+
#   xlab("")+
#   #geom_line()+
#   theme(text=element_text(size=69),
#         legend.text=element_text(face="italic",size=60,colour="black"),
#         axis.text.x  = element_text(face="italic",angle = -90,size=60, hjust = 0, vjust=0.5),
#         axis.text.y  = element_text(face="italic",size=60,colour="black"))+
#   scale_y_continuous(limits = c(-7,7),breaks=seq(-7,7,0.5))+
#   coord_flip()
# #exp_plot2  
# ggsave(paste0(graph_dir,"/","log2foldchange",".pdf"),exp_plot2,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)

########################################
# #Log2foldchange fixed per change in Solid vs Liquid
# df_p <- df
# df_p$rel_LI <- df_p$rel_LI*100;df_p$rel_SO <- df_p$rel_SO*100
# df_p <- df_p[which(df_p$p.adj_Phase3 <= 0.05),]
# 
# df_p <- df_p[which(df_p$change_Treat3 > 0.5 | df_p$change_Treat3 < 0.5),]
# df_p <- df_p[which(!is.na(df_p$Genus)),]
# df_p$anti <- antilog(df_p$change_Treat3,2)
# #Log2 fold change plot
# exp_plot2a <- ggplot(df_p, aes(y=change_Phase3, x=reorder(Genus, -change_Phase3), color=Phylum, fill=Phylum)) + 
#   geom_point(size=25) + 
#   #geom_bar(stat="identity")+
#   geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
#   ylab("Log2foldchange")+
#   xlab("")+
#   #geom_line()+
#   theme(text=element_text(size=60),
#         legend.text=element_text(face="italic",size=60,colour="black"),
#         axis.text.x  = element_text(face="italic",angle = -90,size=60, hjust = 0, vjust=0.5),
#         axis.text.y  = element_text(face="italic",size=60,colour="black"))+
#   scale_y_continuous(limits = c(-7,7),breaks=seq(-7,7,0.5))+
#   coord_flip()
# #exp_plot2  
# ggsave(paste0(graph_dir,"/","log2foldchange_phase",".pdf"),exp_plot2a,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)
# ###
# 
# #https://gitlab.com/s-schmeier/crc-study-2017/blob/master/scripts/qpcr/calc-fc.R
# 
# ##################################################
##################################################
# 
# #TreeFDR
# 
# tree.fdr.obj <- TreeFDR(X=t(ps3@otu_table), Y=ps3@sam_data$Breed, tree=ps3@phy_tree, test.func, perm.func,B=100)
# BH.p.adj <- p.adjust(tree.fdr.obj$p.unadj, 'fdr')
# 
# # Performance measure: compare TreeFDR and BH
# tree.p.adj <- tree.fdr.obj$p.adj
# BH.p.adj <- p.adjust(tree.fdr.obj$p.unadj, 'fdr')
# # names(tree.p.adj) <- names(tree.fdr.obj$p.unadj)
# # names(BH.p.adj) <- names(tree.fdr.obj$p.unadj)
# X <- t(ps3@otu_table)
# X  <- X[,1]
# 
# #X <- as.data.frame(merge(X,tree.p.adj,"row.names"))
# X <- full_join(data.frame(seqId=row.names(X),X= X[,1]),
# data.frame(seqId=names(tree.p.adj),tree.p.adj=as.numeric(tree.p.adj)),by="seqId"
# )
# X <- X[,-2]
# colnames(X) <-c("ROWNAME","StructFDR")
# row.names(X) <- X$ROWNAME
# 
# 
# 
# #X <- merge(X,BH.p.adj,"row.names")
# X <-  full_join(X,
#                 data.frame(ROWNAME=names(BH.p.adj),BH=as.numeric(BH.p.adj)),by="ROWNAME"
# )
# #row.names(X) <- X$Row.names
# #X <- X[,-c(1,2)]
# #colnames(X) <-c("StructFDR","BH")
# 
# #X <- merge(X,ps3@tax_table,"row.names")
# colnames(tax_table_df)[8] <- "ROWNAME"
# #tax_table_df$ROWNAME
# X$ROWNAME <- as.character(X$ROWNAME)
# X <- full_join(X,tax_table_df,"ROWNAME")
# 
# X2 <- X[which(X$StructFDR<=0.05),]
# row.names(X2) <- X2$ROWNAME
# X2$Genus <- droplevels(X2$Genus)
# df_sp <- df
# row.names(df_sp) <- df_sp$Row.names
# val <- X2[,c(2,3)];val$seqId <- row.names(val)
# #df_sp <- merge(df_sp,X2[,c(2,3)],by.x="Row.names",by.y="row.names")
# df_sp <- full_join(df_sp,val,"seqId")
# #df_sp <- df_sp[,-1]
# df_sp$Row.names <- paste0(df_sp$seqId," - ",df_sp$Family," - ",df_sp$Genus)
# #df_sp <- df_sp[which(df_sp$change_Treat3 > 1 | df_sp$change_Treat3 < - 1),]
# #df_sp <- df_sp[which(!is.na(df_sp$Genus)),]
# df_sp <- df_sp[which(!is.na(df_sp$StructFDR)),]
# 
# write.csv(df_sp,paste0(out_dir,"/","csv/","taxa__dfr_TREEFDR.csv"),quote = F,row.names = F)
#     
# 
# df_sp_THR <- as.data.frame(cbind(df_sp$seqId,df_sp$status3_FE,df_sp$status3_Phase,df_sp$HE_gmean3,df_sp$LE_gmean3,df_sp$change_Treat3))
# colnames(df_sp_THR) <-c("seqId","status_FE","status_Phase","HE_gmean3","LE_gmean3","change_Treat3")
# df_sp_THR$HE_gmean3 <- as.numeric(as.character(df_sp_THR$HE_gmean3))
# df_sp_THR$LE_gmean3 <- as.numeric(as.character(df_sp_THR$LE_gmean3))
# df_sp_THR$change_Treat3 <- as.numeric(as.character(df_sp_THR$change_Treat3))
# 
# 
# df_sp_THR$a.greaterthan.b_T3 <- 
#   df_sp_THR$HE_gmean3 > df_sp_THR$LE_gmean3
# exp_plot3 <- ggplot(df_sp, aes(y=change_Treat3, x=reorder(Row.names, -change_Treat3), color=Phylum, fill=Phylum)) + 
#   #geom_point(size=25) + 
#   geom_bar(stat="identity")+
#   geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
#   ylab("Log2foldchange")+
#   xlab("")+
# #  geom_hline()+
#   theme(text=element_text(size=60),
#         legend.text=element_text(face="italic",size=60,colour="black"),
#         axis.text.x  = element_text(face="italic",angle = -90,size=60, hjust = 0, vjust=0.5),
#         axis.text.y  = element_text(face="italic",size=60,colour="black"))+
#   scale_y_continuous(limits = c(-15,4),breaks=seq(-15,4,0.5))+
#   coord_flip()
# #exp_plot2  
# ggsave(paste0(graph_dir,"/","log2foldchange_structFDR",".pdf"),exp_plot3,dpi=300,width =110,height=80,units = "cm",scale=1.2,limitsize = FALSE)
# 
# 
# sum(table(df_sp$Phylum)[order(table(df_sp$Phylum),decreasing=T)])
# table(df_sp$Family[which(df_sp$change_Treat3<0)])[order(table(df_sp$Family[which(df_sp$change_Treat3<0)]),decreasing = T)]
# table(df_sp$Family[which(df_sp$change_Treat3>0)])[order(table(df_sp$Family[which(df_sp$change_Treat3>0)]),decreasing = T)]
#   #table(X2$Genus[order(X2$Genus,decreasing = T)])
# 
# 
# ###DEseq2
# tree.fdr.obj2 <- TreeFDR(X=DESeq_list$abundds2, Y=ps3@sam_data$Treatment, tree=ps3@phy_tree, test.func, perm.func,B=100)
# BH.p.adj2 <- p.adjust(tree.fdr.obj$p.unadj, 'fdr')
# 
# # Performance measure: compare TreeFDR and BH
# tree.p.adj2 <- tree.fdr.obj2$p.adj
# BH.p.adj2 <- p.adjust(tree.fdr.obj2$p.unadj, 'fdr')
# names(tree.p.adj2) <- names(tree.fdr.obj2$p.unadj)
# names(BH.p.adj2) <- names(tree.fdr.obj2$p.unadj)
# X2_0 <- t(ps3@otu_table)
# X2_0  <- X2_0[,1]
# 
# X2_0 <- full_join(data.frame(seqId=row.names(X2_0),X=X2_0[,1]),
#                data.frame(seqId=names(tree.p.adj2),tree.p.adj=as.numeric(tree.p.adj2)),by="seqId"
# )
# X2_0 <- X2_0[,-2]
# colnames(X2_0) <-c("ROWNAME","StructFDR")
# row.names(X2_0) <- X2_0$ROWNAME
# X2_0 <-  full_join(X2_0,
#                 data.frame(ROWNAME=names(BH.p.adj2),BH=as.numeric(BH.p.adj2)),by="ROWNAME"
# )
# 
# 
# X2_0$ROWNAME <- as.character(X2_0$ROWNAME)
# X2_0 <- full_join(X2_0,tax_table_df,"ROWNAME")
# 
# X2 <- X[which(X$StructFDR<=0.05),]
# row.names(X2) <- X2$ROWNAME
# 
# X2_0_FINAL <- X2_0[which(X2_0$StructFDR<0.05),]
# 
# ##################33
# par(mfrow=c(1, 2))
# te <- df$p.adj_T3;names(te) <- row.names(te)
# 
# # plot(ps3@phy_tree, type = 'fan', edge.color = "gray", cex=0.2, tip.color = "black", 
# #      show.tip.label = F, label.offset=0.06) 
# # #tiplabels(text="", tip=which(BH.p.adj <= 0.05), frame="n", pch=4, col="blue")
# # tiplabels(text="", tip=which(te <= 0.05), frame="n", pch=4, col="blue")
# # #
# # #
# # plot(ps3@phy_tree, type = 'fan', edge.color = "gray", cex=0.2, tip.color = "black", 
# #      show.tip.label = F, label.offset=0.06) 
# # tiplabels(text="", tip=which(tree.p.adj <= 0.05), frame="n", pch=4, col="red")
# 
# 
# ###Structured FDR
# el <- phy_tree(ps3)$edge
# el0 <- el
# el0 <- el0[nrow(el):1, ]
# abundds <- t(ps3@otu_table)
# 
# short_names <- row.names(abundds)#
# el_names <- c(short_names, seq_len(phy_tree(ps)$Nnode))
# el[, 1] <- el_names[el0[, 1]]
# el[, 2] <- el_names[as.numeric(el0[, 2])]#
# #unadj_p <- treePValues(el, abundds, sample_data(ps)$Treatment)
# unadj_p <- treePValues(el, as.data.frame(t(ps3@otu_table)), sample_data(ps)$Treatment)
# hfdr_res <- hFDR.adjust(unadj_p, el,alpha = .75)
# 
# tax <- tax_table(ps)[, c("Family", "Genus")] %>%
#   data.frame()
# tax$seq <- short_names
# 
# hfdr_res@p.vals$seq <- rownames(hfdr_res@p.vals)
# tax2 <- tax %>%
#   left_join(hfdr_res@p.vals) %>%
#   arrange(adjp) #%>%
# 
# 
# short_names <- row.names(DESeq_list$abundds2)#
# el_names <- c(short_names, seq_len(phy_tree(ps)$Nnode))
# el[, 1] <- el_names[el0[, 1]]
# el[, 2] <- el_names[as.numeric(el0[, 2])]#
# #unadj_p <- treePValues(el, abundds, sample_data(ps)$Treatment)
# unadj_p2 <- treePValues(el, DESeq_list$abundds2, sample_data(ps)$Treatment)
# 
# tax11 <- tax_table(ps)[, c("Family", "Genus")] %>%
#   data.frame()
# tax11$seq <- short_names
# hfdr_res2 <- hFDR.adjust(unadj_p2, el,alpha = .75)
# 
# hfdr_res2@p.vals$seq <- rownames(hfdr_res2@p.vals)
# tax3 <- tax11 %>%
#   left_join(hfdr_res2@p.vals) %>%
#   arrange(adjp) #%>%
# 
# #####################################################
# #Comparing adjusted p values
# 
# p.adj_T <- as.numeric(df$p.adj_T);names(p.adj_T) <- df$Row.names
# p.adj_T3 <- as.numeric(df$p.adj_T3);names(p.adj_T3) <- df$Row.names
# dat <- sapply(seq(0, 0.5, len=100), function (x) {
#   c(BH_rel=sum(p.adj_T <= x),
#     BH_prop=sum(p.adj_T3<= x),
#     TreeFDR=sum(tree.p.adj <= x),
#     TreeFDR_DS=sum(tree.p.adj2 <= x),
#     StructFDR = sum(tax2$adjp <=x,na.rm = T),
#     StructFDR_DS = sum(tax3$adjp <=x,na.rm = T)) #StructFDR = sum(tax2$padj <=x,na.rm = T
# })
# colnames(dat) <- seq(0, 0.5, len=100)
# dat <- melt(dat)
# colnames(dat) <- c('Method', 'FDR_cutoff', 'OTU_number')
# dat$Method <- factor(dat$Method, levels=c('BH_rel',"BH_prop",'TreeFDR',"TreeFDR_DS","StructFDR","StructFDR_DS"))
# FDR_comparison_plot <- ggplot(dat, aes_string(x='FDR_cutoff', y='OTU_number', group = 'Method', 
#                        shape='Method', linetype='Method')) +
#   geom_line(size=0.9) +
#   geom_point(size=20) +
#   geom_vline(xintercept=0.05,hfdr_res,alpha = 0.5,size=1,linetype = 2,colour="red")+
#   ylab('Number of differential OTUs') +
#   xlab('FDR level') +
#   scale_shape_manual(values=c(15,12, 1,9,6,8)) +
#   scale_x_continuous(limits = c(0,0.5),breaks=seq(0,0.5,0.05))+
#     scale_y_continuous(limits = c(0,300),breaks=seq(0,300,9))+
#   theme(text=element_text(size=60),
#         legend.text=element_text(size=60),
#         axis.text.x  = element_text(size=60,colour="black"),
#         axis.text.y  = element_text(size=60,colour="black"))
# ggsave(paste0(graph_dir,"/","FDR_comparison_plot",".pdf"),FDR_comparison_plot,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)
# 
# #####################################################
#####################################################
#####################################################
#####################################################
#Coverage per ASVs
# rel_abund_total <- melt(df[,c("seqId","HE","LE","count_T")],"seqId")
# rel_abund_total <-left_join(rel_abund_total,df[,c("seqId",firstup(levels))],"seqId"
#                             )

# 
# rel_abund_total <- merge(rel_abund_total,df[,c("seqId",firstup(levels))],"Row.names")


# coverage_taxP2<- ggplot(data=rel_abund_total,aes(x=variable,y=value,fill=Phylum))+
#   geom_bar(stat="identity",show.legend = T,position = "dodge")+
#   xlab("")+
#   ylab("Coverage (%)")+
#   # ylim(0,0.1)+
#   #geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
#   theme(text=element_text(size=60),
#         legend.text=element_text(size=60),
#         axis.text.x  = element_text(size=60,colour="black"),
#         axis.text.y  = element_text(size=60,colour="black"))+
#   coord_flip()+
#   facet_wrap(~Kingdom)
# #  scale_y_continuous(limits = c(0,45),breaks=seq(0,45,5))
# coverage_taxP2
#ggsave(paste0(graph_dir,"/","Taxonomy_coverage_current_perc",".pdf"),coverage_tax_curr_plot2,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)

###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################

#FDR per taxonomic level

tax_levels_df <- 
rbind(ps_object3$pserie_kingdom[,c(1:18,20)],
      ps_object3$pserie_phylum[,c(1:18,21)],
      ps_object3$pserie_class[,c(1:18,22)],
      ps_object3$pserie_order[,c(1:18,23)],
      ps_object3$pserie_family[,c(1:18,24)],
      ps_object3$pserie_genus[,c(1:18,25)])
tax_levels_df$taxon <- NA
tax_levels_df$taxon <- c(as.character(ps_object3$pserie_kingdom[,19]),
                         as.character(ps_object3$pserie_phylum[,20]),
                         as.character(ps_object3$pserie_class[,21]),
                         as.character(ps_object3$pserie_order[,22]),
                         as.character(ps_object3$pserie_family[,23]),
                         as.character(ps_object3$ pserie_genus[,24]))
tax_levels_df$tax_level <-  as.factor(tax_levels_df$tax_level)
tax_levels_df$tax_level <- factor(tax_levels_df$tax_level,
                           levels = c("Kingdom","Phylum","Class","Order","Family","Genus"))


coverage_taxP2<- ggplot(data=tax_levels_df,aes(x=tax_level,y=Abundance,fill=Breed))+
 # geom_bar(stat="identity",show.legend = T,position = "dodge")+
  geom_boxplot(position = "dodge",outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
  stat_boxplot(geom ='errorbar') +
  xlab("")+
  ylab("Proportional relative abundance (%)")+
  scale_fill_manual("legend", values = c("red","blue","green"))+
    # ylim(0,0.1)+
  #geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
  theme(panel.background = element_rect(fill = "gray95"),text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black",angle=90),
        axis.text.y  = element_text(size=60,colour="black"),
        legend.title=element_blank(),legend.position="none"
        )#+

  #coord_flip()
#  scale_y_continuous(limits = c(0,45),breaks=seq(0,45,5))
t_levels <- firstup(levels)
t_levels <- t_levels[1:6]

#################################################################################################################################
#################################################################################################################################

#WILCOXON PER TAXONOMICAL RANK


#################################################################################################################################
#################################################################################################################################

dfs_tax_level <- lapply(1:length(t_levels), function(i){
  cat("  ","\n");cat(i," |",as.character(t_levels[[i]]),"\n");  cat("  ","\n");
X_Lev   <- subset(tax_levels_df,tax_levels_df$tax_level==t_levels[[i]])
taxa <- unique(X_Lev$taxon) 

dfs <- lapply(1:length(taxa),function(j){
cat(i," | ",j," |",taxa[[j]],"\n")
taxon_obj <- taxa[[j]]
tax_subset <- subset(X_Lev,X_Lev$taxon==taxon_obj 
                      )#$Abundance
#tax_subset <- tax_subset[duplicated(tax_subset$sampleID),]


#tax_subset <- tax_subset[duplicated(tax_subset$sampleID),]
HE_S <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                 X_Lev$ADG_Treatment=="High feed efficiency" &
                 X_Lev$Phase=="Solid" )$Abundance
HE_L <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                 X_Lev$ADG_Treatment=="High feed efficiency" &
                 X_Lev$Phase=="Liquid" )$Abundance
LE_S <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                 X_Lev$ADG_Treatment=="Low feed efficiency" &
                 X_Lev$Phase=="Solid" )$Abundance
LE_L <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                 X_Lev$ADG_Treatment=="Low feed efficiency" &
                 X_Lev$Phase=="Liquid" )$Abundance
HE_Sub <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                 X_Lev$ADG_Treatment=="High feed efficiency")$Abundance

LE_Sub <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                   X_Lev$ADG_Treatment=="Low feed efficiency")$Abundance

S_Sub <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                   X_Lev$Phase=="Solid")$Abundance

L_Sub <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                   X_Lev$Phase=="Liquid")$Abundance

B1_Sub <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                  X_Lev$Breed=="Cheviot")$Abundance
B2_Sub <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                  X_Lev$Breed=="Connemara")$Abundance
B3_Sub <- subset(X_Lev,X_Lev$taxon==taxon_obj & 
                  X_Lev$Breed=="Lanark")$Abundance

x <- data.frame(var=as.numeric(c(B1_Sub,B2_Sub,B3_Sub)),
                Breed=c(rep("Connemara",length(B1_Sub)),rep("Cheviot",length(B2_Sub)),rep("Lanark",length(B3_Sub))
                ))
w <-agricolae::kruskal( x$var, x$Breed, alpha = 0.05,group=F)
require(plotrix)
df_pvalues <- data.frame(
  # High_mean = paste0(round(mean(HE_Sub),3),"±",round(std.error(HE_Sub),3)),
  # Low_mean = paste0(round(mean(LE_Sub),3),"±",round(std.error(LE_Sub),3)),
  
High_vs_Low3 = if(sum(HE_Sub)>0 & sum(LE_Sub)>0){High_vs_Low3 = wilcox.test(x=HE_Sub,y=LE_Sub,paired=F,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value} else {High_vs_Low3 = 1},
#Solid_vs_liquid3 =if(sum(S_Sub)>0 & sum(L_Sub)>0){Solid_vs_liquid3 = wilcox.test(x=S_Sub,y=L_Sub,paired=T,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value} else {Solid_vs_liquid3 = 1},
Solids3 = if(sum(HE_S)>0 & sum(LE_S)>0){Solids3 = wilcox.test(x=HE_S,y=LE_S,paired=F,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value} else {Solids3 = 1},
Liquids3 =if(sum(HE_L)>0 & sum(LE_L)>0){Liquids3 = wilcox.test(x=HE_L,y=LE_L,paired=F,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value} else {Liquids3 = 1},
#High_S_vs_L3 = if(sum(HE_S)>0 & sum(HE_L)>0){High_S_vs_L3 = wilcox.test(x=HE_S,y=HE_L,paired=T,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value} else {High_S_vs_L3 = 1},
#Low_S_vs_L3 = if(sum(LE_S)>0 & sum(LE_L)>0){Low_S_vs_L3 = wilcox.test(x=LE_S,y=LE_L,paired=T,alternative="two.sided",exact=F,correct=F,conf.int=T,conf.level=0.95)$p.value} else {Low_S_vs_L3 = 1},
Breed_pvalue = as.numeric(w$statistics[3]),

count_per_sample = length(unique(tax_subset[which(tax_subset$Abundance>0),]$Sample)),
  #if(sum(tax_subset!=0)>nrow(ps@otu_table)) {count_per_sample =nrow(ps@otu_table)} else {count_per_sample = sum(tax_subset!=0)},#sum(HE_S>0) + sum(HE_L>0) + sum(LE_S>0) + sum(LE_L>0),
count_per_sample_HE = length(unique(as.character(subset(X_Lev,X_Lev$taxon==taxon_obj & X_Lev$Treatment=="High feed efficiency"  & Abundance >0)$Sample))),#sum(HE_Sub>0),
count_per_sample_LE=  length(unique(as.character(subset(X_Lev,X_Lev$taxon==taxon_obj & X_Lev$Treatment=="Low feed efficiency" & Abundance >0)$Sample))),#sum(LE_Sub>0),
count_per_sample_SO = length(unique(as.character(subset(X_Lev,X_Lev$taxon==taxon_obj & X_Lev$Phase=="Solid" & Abundance >0)$Sample))),#sum(S_Sub>0),
count_per_sample_LI = length(unique(as.character(subset(X_Lev,X_Lev$taxon==taxon_obj & X_Lev$Phase=="Liquid" & Abundance >0)$Sample))),#sum(L_Sub>0)

count_per_sample_B1 = length(unique(as.character(subset(X_Lev,X_Lev$taxon==taxon_obj & X_Lev$Breed=="Connemara" & Abundance >0)$Sample))),#sum(S_Sub>0),
count_per_sample_B2 = length(unique(as.character(subset(X_Lev,X_Lev$taxon==taxon_obj & X_Lev$Breed=="Cheviot" & Abundance >0)$Sample))),#sum(S_Sub>0),
count_per_sample_B3 = length(unique(as.character(subset(X_Lev,X_Lev$taxon==taxon_obj & X_Lev$Breed=="Lanark" & Abundance >0)$Sample)))#sum(S_Sub>0),


)

df_pvalues[2,] <- NA
#df_pvalues[2,1:2] <-NA 
#df_pvalues[2,3:8] <- p.adjust(p = df_pvalues[1,3:8] , method =  "BH",n = 6)
df_pvalues[2,1:4] <- p.adjust(p = df_pvalues[1,1:4] , method =  "BH",n = 4)
df_pvalues <- as.data.frame(t(df_pvalues))
#df_pvalues[,2] <- NA;
df_pvalues[,2]  <- as.numeric(as.character(df_pvalues[,2]))
df_pvalues$sig <- NA
#df_pvalues$sig[which(df_pvalues[,2] <= 0.05)] <-"**"; df_pvalues$sig[which(df_pvalues[,2] > 0.05)] <-"-"; 
colnames(df_pvalues) <- c("p.value","p.adj","sig")#,"HE_counts","LE_counts","Solid_counts","Liquid_counts")
df_pvalues$taxon <- taxon_obj
df_pvalues$comparison <- row.names(df_pvalues)
df_pvalues$taxlevel <- t_levels[[i]]
df_pvalues <- df_pvalues[,c(6,4,5,1:3)]

return(df_pvalues)

})

dfs <- do.call(rbind,dfs)

#dfs$count_per_sample <- NA

return(dfs)
cat("  ","\n")

})


dfs_tax_level2 <- do.call(rbind,dfs_tax_level)
#dfs_tax_level3 <- dfs_tax_level2
#dfs_tax_level2 <- dfs_tax_level3
dfs_tax_level <- dfs_tax_level2[which(dfs_tax_level2$comparison !="count_per_sample" &
                                        dfs_tax_level2$comparison !="count_per_sample_HE" &
                                        dfs_tax_level2$comparison !="count_per_sample_LE" &
                                        dfs_tax_level2$comparison !="count_per_sample_SO" &
                                        dfs_tax_level2$comparison !="count_per_sample_LI" &
                                        dfs_tax_level2$comparison !="count_per_sample_B1" &
                                        dfs_tax_level2$comparison !="count_per_sample_B2" &
                                        dfs_tax_level2$comparison !="count_per_sample_B3" 
                                        
                                         ),]
dfs_tax_level_count <- dfs_tax_level2[which(dfs_tax_level2$comparison !="High_vs_Low3" &
                                        dfs_tax_level2$comparison !="Solid_vs_liquid3" &
                                        dfs_tax_level2$comparison !="Solids3" &
                                        dfs_tax_level2$comparison !="Liquids3" &
                                        dfs_tax_level2$comparison !="High_S_vs_L3" &
                                        dfs_tax_level2$comparison !="Low_S_vs_L3" &
                                        dfs_tax_level2$comparison !="Breed_pvalue" 
                                          
),]

dfs_tax_level_count <- dfs_tax_level_count[,1:4]
colnames(dfs_tax_level_count)[4] <- "count"

row.names(dfs_tax_level) <- 1: nrow(dfs_tax_level)
dfs_tax_level$p.adj_total <- NA
dfs_tax_level$p.adj_total <- p.adjust(p = na.omit(dfs_tax_level$p.value), method =  "BH",n = length(na.omit(dfs_tax_level$p.value)))
dfs_tax_level$sig_total <- NA;
dfs_tax_level$sig_total[which(dfs_tax_level$p.adj_total <= 0.05)] <-"**"; dfs_tax_level$sig_total[which(dfs_tax_level$p.adj_total > 0.05)]  <-"-";

dfs_tax_level <- dfs_tax_level[which(dfs_tax_level$sig_total=="**"),]

#dfs_tax_level[which(dfs_tax_level$taxlevel=="Phylum"),]
coverage_taxP3<- ggplot(data=dfs_tax_level[which(dfs_tax_level$taxlevel=="Phylum"),],
                        aes(x=comparison,y=p.adj_total,fill=taxon,colour=taxon))+
   geom_bar(stat="identity",show.legend = T,position = "dodge")+
  #geom_boxplot(position = "dodge",outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
  stat_boxplot(geom ='errorbar') +
  xlab("")+
  ylab("Proportional relative abundance (%)")+
  #scale_fill_manual("legend", values = c("red","blue"))+
  # ylim(0,0.1)+
  #geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
  theme(panel.background = element_rect(fill = "gray95"),text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black",angle=90),
        axis.text.y  = element_text(size=60,colour="black"))+
  coord_flip()#,
#        legend.title=element_blank(),legend.position="none"
#+
coverage_taxP3
write.csv(dfs_tax_level,paste0(out_dir,"/","csv/","taxa_fdr2.csv"),quote = F,row.names = F)
write.csv(dfs_tax_level_count,paste0(out_dir,"/","csv/","taxa_fdr2_count.csv"),quote = F,row.names = F)
####
############################################################################################################
############################################################################################################
############################################################################################################
#http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html
test_df <- data.frame(Taxa=as.character(row.names(ps@tax_table)), ps@tax_table)
test_df$Taxa <- as.character(test_df$Taxa)

composition_list_tax <- lapply(1:6,function(j){
  cat("  ","\n");cat(firstup(levels)[j],"\n")
x <- ps_object[[j]]@otu_table
x<-x/rowSums(x)
x<-x[,order(colSums(x),decreasing=TRUE)]
#Extract list of top N Taxa
N<-16
taxa_list<-colnames(x)[1:N]
taxa_list <- taxa_list[which(!is.na(taxa_list))]
#Generate a new table with everything added to Others
#new_x<-data.frame(x[,colnames(x) %in% taxa_list],
#Others=rowSums(x[,!colnames(x) %in% taxa_list]))

new_x<-data.frame(x[,colnames(x) %in% taxa_list],
                  
                  Others= if(length(colnames(x))<N){Others=NA} else { rowSums(x[,!colnames(x) %in% taxa_list])}
)
                   
#You can change the Type=grouping_info[,1] should you desire any other grouping of panels
df<-NULL
for (i in 1:dim(new_x)[2]){
  cat(i,"\n")
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),sampleID = ps_object[[j]]@sam_data$sampleID,
                  Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],
                  Type=ps_object[[j]]@sam_data$ADG_Treatment,Phase= ps_object[[j]]@sam_data$Phase,
                  Breed=ps_object[[j]]@sam_data$Breed)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
};rm(i)
df$Taxa <- as.character(df$Taxa)


df <- left_join(df,test_df,by="Taxa")
#df <- merge(df,ps@tax_table,by.x="Taxa",by.y="row.names",all=T)
df <- df[which(!is.na(df$Value)),]
df <- df[,c(1:7,match(firstup(levels[[j]]),colnames(df)))]
colnames(df)[ncol(df)] <- "Taxon"
#for(i in 5:11){
  df[,8] <- as.character(df[,8])
  df[,8][is.na(df[,8])] <- "Others"
  df[,8] <- factor(df[,8],unique(df[,8]))
#};rm(i)

return(df)
})

names(composition_list_tax) <- firstup(levels)[1:6]
#############################################
###BARPLOT
#############################################
lapply(1:length(composition_list_tax),function(i){
  cat(i,"\n")
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00");
N <- 17
  library(ggplot2)
p<-ggplot(composition_list_tax[[i]],aes(factor(sampleID),Value,fill=Taxon))+
  geom_bar(stat="identity")+
  facet_grid(. ~ Breed + Phase, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values=colours[1:(N+1)])
p<-p+theme_bw()+ylab("Relative abundances")
p<-p+ scale_y_continuous(expand = c(0,0))+
  theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines"))
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=60),
           text=element_text(size=60),
           axis.text.y  = element_text(size=60,colour="black"))+
  xlab("Sample")
  
  ggsave(paste0(graph_dir,"/","compoplot_",names(composition_list_tax)[[i]],".pdf"),p,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)
  
# pdf("TAXAplot.pdf",height=6,width=21)
# print(p)
# dev.off()
# p
})
#############################################
###BOXPLOT
#############################################

lapply(1:length(composition_list_tax),function(i){
  cat(i,"\n")
  colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00");
  N <- 17
  library(ggplot2)
  p<-ggplot(composition_list_tax[[i]],aes(Taxon,Value*100,fill=Breed))+
    # geom_bar(stat="identity")+
    geom_boxplot(aes(fill=Breed)) +
    stat_boxplot(geom ='errorbar') +
    ylab("Relative abundance (%)")+
    # stat_summary(fun.y=mean, geom="line", aes(group=1))  +
    # stat_summary(fun.y=mean, geom="point")+
    scale_fill_manual("Breed", values = c("blue","red","green"))
  
     
  #p<-p+scale_fill_manual(values=colours[1:(N+1)])
  
    
  # p<-p+ scale_y_continuous(expand = c(0,0))+
  #   theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines"))
  p <- p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=60),
             text=element_text(size=60),
             axis.text.y  = element_text(size=60,colour="black"))
  
  #p <-p + facet_grid(. ~ Phase, drop=TRUE,scale="free",space="free_x")
  
  ggsave(paste0(graph_dir,"/","compoplot_",names(composition_list_tax)[[i]],"BP_abundance.pdf"),p,dpi=300,width =150,height=90,units = "cm",scale=1.2,limitsize = FALSE)
  
  # pdf("TAXAplot.pdf",height=6,width=21)
  # print(p)
  # dev.off()
  # p
})
#############################################
#############################################
# p<-ggplot(composition_list_tax$Kingdom,aes(factor(sampleID),Value,fill=Taxon))+
#   geom_bar(stat="identity")+
#   facet_grid(. ~ Type + Phase, drop=TRUE,scale="free",space="free_x")
# p<-p+scale_fill_manual(values=c("gray93","gray85"))
# p<-p+theme_bw()+ylab("Relative abundances")
# p<-p+ scale_y_continuous(expand = c(0,0))+
#   theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines"))
# p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=60),
#            text=element_text(size=60),
#            axis.text.y  = element_text(size=60,colour="black"))+
#   xlab("Sample")
# 
# ggsave(paste0(graph_dir,"/","compoplot_",names(composition_list_tax)[[1]],".pdf"),p,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)

#############################################
#############################################


composition_list_tax2 <- composition_list_tax
for(i in 1:length(composition_list_tax2)){
  composition_list_tax2[[i]]$Value[which(composition_list_tax2[[i]]$Value>0)] <- 1
  };rm(i)



# 
# FDR_list <- list(
#   BH_prop = df,
#   #Tree_FDR_DESeq = X2_0_FINAL,
#    #Tree_FDR_PRA = df_sp,
#   composition_list_tax = composition_list_tax,
#   composition_list_tax_bin = composition_list_tax2
# )
# 
# saveRDS(FDR_list,paste0(out_dir,"/","csv/","FDR_list.RDS"))
# 
# }else {
#   FDR_list <- readRDS(paste0(out_dir,"/","csv/","FDR_list.RDS"))
# }
#   HE_ps <- subset_taxa(ps3,ps3@sam_data$Treatment=="High feed efficiency")
# LE_ps <- subset_taxa(ps3,ps3@sam_data$Treatment=="Low feed efficiency")
# 
# 
# #Average counts per sample
# mean(sample_sums(ps))
# #Standard deviation of count per sample
# sd(sample_sums(ps))
# 
# FDR_list
# 
# #Kingdom total reads
# king_sp_prod<- table(tax_table(ps3)[,1])#table(tax_table(ps3)[,1])
# king_sp_prod <- king_sp_prod[order(-king_sp_prod)]
# king_sp_prod_a <- king_sp_prod
# king_sp_prod_a <- (king_sp_prod_a/sum(table(tax_table(ps3)[,1])))*100
# king_sp_prod_a
# #Kingdom HE
# king_sp_prod_1 <- table(tax_table(HE_ps)[,1])#table(tax_table(ps3)[,1])
# king_sp_prod_1 <- king_sp_prod_1[order(-king_sp_prod_1)]
# king_sp_prod_1_1 <- king_sp_prod_1
# king_sp_prod_1_1 <- (king_sp_prod_1_1/sum(table(tax_table(HE_ps)[,1])))*100
# king_sp_prod_1_1
# #Kingdom LE
# king_sp_prod_2 <- table(tax_table(LE_ps)[,1])#table(tax_table(ps3)[,1])
# king_sp_prod_2 <- king_sp_prod_2[order(-king_sp_prod_2)]
# king_sp_prod_2_2 <- king_sp_prod_2
# king_sp_prod_2_2 <- (king_sp_prod_2_2/sum(table(tax_table(LE_ps)[,1])))*100
# king_sp_prod_2_2
# #########################3
# #Family total reads
# fam_sp_prod<- table(tax_table(ps3)[,5])#table(tax_table(ps3)[,1])
# fam_sp_prod <- fam_sp_prod[order(-fam_sp_prod)]
# fam_sp_prod_a <- fam_sp_prod
# fam_sp_prod_a <- (fam_sp_prod_a/sum(table(tax_table(ps3)[,5])))*100
# fam_sp_prod_a
# #Family HE
# fam_sp_prod_1 <- table(tax_table(HE_ps)[,5])#table(tax_table(ps3)[,1])
# fam_sp_prod_1 <- fam_sp_prod_1[order(-fam_sp_prod_1)]
# fam_sp_prod_1_1 <- fam_sp_prod_1
# fam_sp_prod_1_1 <- (fam_sp_prod_1_1/sum(table(tax_table(HE_ps)[,5])))*100
# fam_sp_prod_1_1
# #Family LE
# fam_sp_prod_2 <- table(tax_table(LE_ps)[,5])#table(tax_table(ps3)[,1])
# fam_sp_prod_2 <- fam_sp_prod_2[order(-fam_sp_prod_2)]
# fam_sp_prod_2_2 <- fam_sp_prod_2
# fam_sp_prod_2_2 <- (fam_sp_prod_2_2/sum(table(tax_table(LE_ps)[,5])))*100
# fam_sp_prod_2_2
# 

# pseq.rel <- microbiome::transform(ps, "compositional")
# core.taxa.standard <- core_members(pseq.rel, detection = 0)
# pseq.core <- core(pseq.rel, detection = 0, prevalence = .5)
# core.taxa <- taxa(pseq.core)
# core.abundance <- sample_sums(core(pseq.rel, detection = .01, prevalence = .05))
# # With compositional (relative) abundances
# det <- c(0, 0.1, 0.5, 2, 5, 20)/100
# prevalences <- seq(.05, 1, .05)
# plot_core(pseq.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")
# 
# ########################################################
########################################################
# 
# # Core with compositionals:
# prevalences <- seq(.05, 1, .05)
# detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
# 
# # Also define gray color palette
# gray <- gray(seq(0,1,length=5))
# p <- plot_core(pseq.rel, plot.type = "heatmap", colours = gray,
#                prevalences = prevalences, detections = detections) +
#   xlab("Detection Threshold (Relative Abundance (%))")
# print(p)    


# Core with absolute counts and horizontal view:
# and minimum population prevalence (given as percentage)
# detections <- 10^seq(log10(1), log10(max(abundances(ps))/10), length = 10)
# 
# library(RColorBrewer)
# p <- plot_core(ps, plot.type = "heatmap", 
#                prevalences = prevalences,
#                detections = detections,
#                colours = rev(brewer.pal(5, "Spectral")),
#                min.prevalence = .2, horizontal = TRUE)
# print(p)
# 
# 
# tab <- global(ps, index = "all")
# sample_data(ps)$diversity <- tab$diversities_shannon
# p <- plot_regression(diversity ~ FCR,data=meta(ps),show.points=T,B=1000,color="red") + xlab("FCR") + ylab("Diversity")
# p
# 
# 
# tab <- microbiome::coverage(ps, threshold = 0.5)
# tab <- inequality(ps)
# tab <- evenness(ps, "all")
# b.pla <- microbiome::divergence(subset_samples(ps, Treatment == "High feed efficiency"))
# b.pla2 <- microbiome::divergence(subset_samples(ps, Treatment == "Low feed efficiency"))
# 
# boxplot(list(HFE = b.pla, LFE = b.pla2))
# 
# 
# # Pick relative abundances (compositional) and sample metadata 
# pseq.rel <- microbiome::transform(ps, "compositional")
# otu <- abundances(pseq.rel)
# meta <- meta(pseq.rel)
# p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "Treatment", size = 3)
# print(p)
# 
# 
# library(vegan)
# permanova <- adonis(t(otu) ~ Treatment,
#                     data = meta, permutations=1000, method = "bray")
# 
# # P-value
# print(as.data.frame(permanova$aov.tab)["Treatment", "Pr(>F)"])
# 
# dist <- vegdist(t(otu))
# anova(betadisper(dist, meta$Treatment))
# 
# 
# coef <- coefficients(permanova)["Treatment1",]
# top.coef <- coef[rev(order(abs(coef)))[1:20]]
# par(mar = c(3, 8, 2, 1))
# barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
# 
# 
# library("randomForest")
# library("plyr") # for the "arrange" function
# library("rfUtilities") # to test model significance
# library("caret")
# #https://github.com/LangilleLab/microbiome_helper/wiki/Random-Forest-Tutorial
