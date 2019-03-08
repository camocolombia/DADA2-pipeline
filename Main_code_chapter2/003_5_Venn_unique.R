
##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##VENN DIAGRAMS                                  #
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
# library("vegan")
# library("ggvegan")
library("dplyr")
library("ggrepel")
library("ggtree")
# library("data.table")
# library("DESeq2")
# library("BSDA")
# library("adegenet")
#install.packages("devtools")
# library("devtools")
# library("phytools")
# library("massMap")
# library("structSSI")
# library("microbiomeViz");
# library("stringr")
# library("StructFDR")
# library("gtools")
library("VennDiagram")
library("ggpubr")
library("tidyverse")
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


##################################################
##################################################
#FDR_list <- readRDS(paste0(out_dir,"/","csv/","FDR_list.RDS"))

ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))
ps@sam_data$Treatment <- factor(ps@sam_data$Treatment,levels = c("High feed efficiency","Low feed efficiency"))
ps3 <-  readRDS(paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))


##################################################
#df <- FDR_list$BH_prop

df <- read.csv(paste0(csv_dir,"/","df_logchange.csv"),header = T)


taxa_df <- data.frame(seqId=row.names(ps3@tax_table),ps3@tax_table)
taxa_df$seqId <- as.character(taxa_df$seqId)

##################################################
##################################################
#Defining dataframe for venn diagrams

B1 <- otu_table(subset_samples(ps, Breed=="Connemara"))
B2 <- otu_table(subset_samples(ps, Breed=="Cheviot"))
B3 <- otu_table(subset_samples(ps, Breed=="Lanark"))

B1_LI <- otu_table(subset_samples(ps, Breed=="Connemara" & Phase=="Liquid"))
B2_LI <- otu_table(subset_samples(ps, Breed=="Cheviot" & Phase=="Liquid"))
B3_LI <- otu_table(subset_samples(ps, Breed=="Lanark" & Phase=="Liquid"))

B1_SO <- otu_table(subset_samples(ps, Breed=="Connemara" & Phase=="Solid"))
B2_SO <- otu_table(subset_samples(ps, Breed=="Cheviot" & Phase=="Solid"))
B3_SO <- otu_table(subset_samples(ps, Breed=="Lanark" & Phase=="Solid"))


LI <- otu_table(subset_samples(ps, Phase=="Liquid"))
SO <- otu_table(subset_samples(ps, Phase=="Solid"))

df_venn2 <-data.frame(
  seqId= colnames(B1),
  B1 = colSums(B1),
  B2 = colSums(B2),
  B3 = colSums(B3),
  SO = colSums(SO),
  LI = colSums(LI),
  B1_LI = colSums(B1_LI),
  B2_LI = colSums(B2_LI),
  B3_LI = colSums(B3_LI),
  B1_SO = colSums(B1_SO),
  B2_SO = colSums(B2_SO),
  B3_SO = colSums(B3_SO)
)
df_venn <- df_venn2
for(i in 2:12){df_venn[,c("seqId","B1","B2","B3","SO","LI",
                         "B1_LI","B2_LI","B3_LI",
                         "B1_SO","B2_SO","B3_SO")]
  df_venn[,i][which(df_venn[,i]>0)] <- 1
}


grid.newpage()
veen_breed <- draw.triple.venn(
  area1 = nrow(subset(df_venn, B1 == 1)), 
  area2 = nrow(subset(df_venn, B2 ==1)),
  area3 = nrow(subset(df_venn, B3 == 1)),

  n12 = nrow(subset(df_venn, B1 == 1 & B2 ==1)), 
  n13 = nrow(subset(df_venn, B1 == 1 & B3 ==1)), 
  n23 = nrow(subset(df_venn, B2 == 1 & B3 ==1)), 
  n123 = nrow(subset(df_venn, B1 == 1 & B2 ==1 & B3 ==1)), 
   category = c("Cheviot","Connemara",
               "Lanark"),
  lty = "blank", 
  cex=6,
  cat.cex = 6,
  fill =  c("blue","red","green")
  )
ggsave(paste0(graph_dir,"/","veen_breed",".pdf"),veen_breed,dpi=300,width =160,height=90,units = "cm",scale=1.2,limitsize = FALSE)

#####################################################
#LIQUID

grid.newpage()
veen_breed_LI <- draw.triple.venn(
  area1 = nrow(subset(df_venn, B1_LI == 1)), 
  area2 = nrow(subset(df_venn, B2_LI ==1)),
  area3 = nrow(subset(df_venn, B3_LI == 1)),
  
  n12 = nrow(subset(df_venn, B1_LI == 1 & B2_LI ==1)), 
  n13 = nrow(subset(df_venn, B1_LI == 1 & B3_LI ==1)), 
  n23 = nrow(subset(df_venn, B2_LI == 1 & B3_LI ==1)), 
  n123 = nrow(subset(df_venn, B1_LI == 1 & B2_LI ==1 & B3_LI ==1)), 
  category = c("Cheviot","Connemara",
               "Lanark"),
  lty = "blank", 
  cex=6,
  cat.cex = 6,
  fill =  c("blue","red","green")
)
ggsave(paste0(graph_dir,"/","veen_breed_liquid",".pdf"),veen_breed_LI,dpi=300,width =160,height=90,units = "cm",scale=1.2,limitsize = FALSE)

#####################################################
#SOLID

grid.newpage()
veen_breed_SO <- draw.triple.venn(
  area1 = nrow(subset(df_venn, B1_SO== 1)), 
  area2 = nrow(subset(df_venn, B2_SO ==1)),
  area3 = nrow(subset(df_venn, B3_SO == 1)),
  
  n12 = nrow(subset(df_venn, B1_SO == 1 & B2_SO ==1)), 
  n13 = nrow(subset(df_venn, B1_SO == 1 & B3_SO ==1)), 
  n23 = nrow(subset(df_venn, B2_SO == 1 & B3_SO==1)), 
  n123 = nrow(subset(df_venn, B1_SO == 1 & B2_SO ==1 & B3_SO ==1)), 
  category = c("Cheviot","Connemara",
               "Lanark"),
  lty = "blank", 
  cex=6,
  cat.cex = 6,
  fill =  c("blue","red","green")
)

ggsave(paste0(graph_dir,"/","veen_breed_solid",".pdf"),veen_breed_SO,dpi=300,width =160,height=90,units = "cm",scale=1.2,limitsize = FALSE)

# df_venn$UNIQ_SOLID <- rowSums(df_venn[,c("SH","SL")]) ; df_venn$UNIQ_SOLID[which(df_venn$UNIQ_SOLID <2)] <- 0;df_venn$UNIQ_SOLID[which(df_venn$UNIQ_SOLID==2)] <- 1
# df_venn$UNIQ_LIQUID <- rowSums(df_venn[,c("LH","LL")]) ; df_venn$UNIQ_LIQUID[which(df_venn$UNIQ_LIQUID <2)] <- 0;df_venn$UNIQ_LIQUID[which(df_venn$UNIQ_LIQUID==2)] <- 1
# 
# df_venn$SOLID <- NA;df_venn$SOLID <- rowSums(df_venn[,c("SH","SL")]); df_venn$SOLID[which(df_venn$SOLID <2)] <- 0;df_venn$SOLID[which(df_venn$SOLID==2)] <- 1
# df_venn$LIQUID <- NA;df_venn$LIQUID <- rowSums(df_venn[,c("LH","LL")]);df_venn$LIQUID[which(df_venn$LIQUID <2)] <- 0;df_venn$LIQUID[which(df_venn$LIQUID==2)] <- 1
# df_venn$SOLID_Vs_LIQUID <- NA;df_venn$SOLID_Vs_LIQUID <- rowSums(df_venn[,c("UNIQ_SOLID","UNIQ_LIQUID")]);df_venn$SOLID_Vs_LIQUID[which(df_venn$SOLID_Vs_LIQUID <2)] <- 0;df_venn$SOLID_Vs_LIQUID[which(df_venn$SOLID_Vs_LIQUID==2)] <- 1
# df_venn$HIGH <- NA;df_venn$HIGH <- rowSums(df_venn[,c("SH","LH")]);df_venn$HIGH[which(df_venn$HIGH <2)] <- 0;df_venn$HIGH[which(df_venn$HIGH==2)] <- 1
# df_venn$LOW <- NA;df_venn$LOW <- rowSums(df_venn[,c("SL","LL")]);df_venn$LOW[which(df_venn$LOW <2)] <- 0;df_venn$LOW[which(df_venn$LOW==2)] <- 1
# df_venn$EFF <- NA;df_venn$EFF <- rowSums(df_venn[,c("HE","LE")]);df_venn$EFF[which(df_venn$EFF <2)] <- 0;df_venn$EFF[which(df_venn$EFF==2)] <- 1
# 
# #colSums(df_venn[,c(2:15)])
df_venn_graph <- left_join(df_venn,taxa_df,by="seqId")

write.csv(df_venn_graph,paste0(out_dir,"/","csv/","df_venn.csv"),quote = F,row.names = F)

 ######################
 prop <- data.frame(
   area1 = nrow(subset(df_venn, B1 == 1)), 
   area2 = nrow(subset(df_venn, B2 ==1)),
   area3 = nrow(subset(df_venn, B3 == 1)),
   
   n12 = nrow(subset(df_venn, B1 == 1 & B2 ==1)), 
   n13 = nrow(subset(df_venn, B1 == 1 & B3 ==1)), 
   n23 = nrow(subset(df_venn, B2 == 1 & B3 ==1)), 
   n123 = nrow(subset(df_venn, B1 == 1 & B2 ==1 & B3  ==1)))
# 
 
# 
# 
# 
category = c("Cheviot","Connemara",
             "Lanark")
prop2 <- prop
# 
 colnames(prop2) <- c(

 "Cheviot",#
   "Connemara",#
   "Lanark",
 "Cheviot-Connemara",
 "Cheviot-Lanark",
 "Connemara-Lanark",
 "All breedds"
 )

 
# 
 prop2$Cheviot_unique <- nrow(subset(df_venn[,1:4], B1 == 1 & B2 !=1 & B3 !=1))# prop2$LE_UNIQUE <- nrow(subset(df_venn, SH != 1 & LH !=1 & SL ==1 & LL==1))
 prop2$Connemara_unique <- nrow(subset(df_venn[,1:4], B1 != 1 & B2 ==1 & B3 !=1))
 prop2$Lanark_unique <- nrow(subset(df_venn[,1:4], B1 != 1 & B2 !=1 & B3 ==1))


 prop2[2,] <- (prop2[1,]/nrow(df_venn))*100


 prop2$ITEM <- NA
 prop2$ITEM[[1]] <- "COUNTS"
 prop2$ITEM[[2]] <- "PERC"

 prop2 <- prop2[,c(ncol(prop2),1:(ncol(prop2)-1))]
 prop2 <-t(prop2)
 colnames(prop2) <- prop2[1,]
 prop2 <- prop2[-1,]
 write.csv(prop2,paste0(out_dir,"/","csv/","Venn_areas_group.csv"),quote = F,row.names = T)


 df_venn_unique <- df_venn_graph[,-c(5:12)]
 B1_sub <-subset(df_venn_unique, B1 == 1 & B2 !=1 & B3 !=1)
 B1_sub$Breed <- "Cheviot"
 
 B2_sub <-subset(df_venn_unique, B1 != 1 & B2 ==1 & B3 !=1)
 B2_sub$Breed <- "Connemara"
 B3_sub <-subset(df_venn_unique, B1 != 1 & B2 !=1 & B3 ==1)
 B3_sub$Breed <- "Lanark"
 
 B_UNIQUE <- rbind(B1_sub,B2_sub,B3_sub)
 B_UNIQUE <- B_UNIQUE[,-c(2,3,4)]
 
 write.csv(B_UNIQUE,paste0(out_dir,"/","csv/","unique_per_breed.csv"),quote = F,row.names = F)

##########High Vs Low

# 
# df_venn_graph <- df_venn
# df_venn_graph$seqId <- as.character(df_venn_graph$seqId)
# #df_venn_graph <- merge(df_venn_graph,ps3@tax_table,by.x="seqId",by.y="row.names",sort=F)
# df_venn_graph <- left_join(df_venn_graph,taxa_df,by="seqId")
#   
# #    # merge(df_venn_graph,ps3@tax_table,by.x="seqId",by.y="row.names",sort=F)
# # df_venn_graph_HE <- df_venn_graph[which(df_venn_graph$HE==1 & df_venn_graph$LE==0),];
# # df_venn_graph_HE$Group <- "High feed efficiency";df_venn_graph_HE <- df_venn_graph_HE[,c("seqId",colnames(ps3@tax_table),"Group")]
# # df_venn_graph_LE <- df_venn_graph[which(df_venn_graph$HE==0 & df_venn_graph$LE==1),]
# # df_venn_graph_LE$Group <- "Low feed efficiency";df_venn_graph_LE <- df_venn_graph_LE[,c("seqId",colnames(ps3@tax_table),"Group")]
# # 
# # df_venn_graph <- rbind(df_venn_graph_HE,df_venn_graph_LE)
# # 
# # write.csv(df_venn_graph,paste0(out_dir,"/","csv/","unique_per_group.csv"),quote = F,row.names = F)
# # 
# # df_venn_graph$Family <-droplevels(df_venn_graph$Family)
# # 
# # df_venn_graph_melt <- melt(table(df_venn_graph$Family,df_venn_graph$Group))
# # 
# # 
# # df_venn_graph_values_HE <- ggplot(data=df_venn_graph_melt,aes(x=reorder(Var1, value),y=value,fill=Var2))+
# #   geom_bar(stat="identity",show.legend = F)+
# #   xlab("Family")+
# #   ylab("Count")+
# #   # ylim(0,0.1)+
# #   scale_fill_manual("legend", values = c("red","blue"))+
# #   theme(text=element_text(size=60),
# #         legend.text=element_text(size=60),
# #         axis.text.x  = element_text(size=60,colour="black",angle=90),
# #         axis.text.y  = element_text(size=60,colour="black"))+
# #   coord_flip()
# # #df_venn_graph_values_HE
# # 
# # 
# # ggsave(paste0(graph_dir,"/","df_venn_graph_values",".pdf"),df_venn_graph_values_HE,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)
# 
# 
# ##########High solid Vs Low solid
# 
# df_venn_graph_S <- df_venn
# df_venn_graph_S$seqId <- as.character(df_venn_graph_S$seqId)
# #df_venn_graph_S <- merge(df_venn_graph_S,ps3@tax_table,by.x="seqId",by.y="row.names")
# df_venn_graph_S <- left_join(df_venn_graph_S,taxa_df,by="seqId")
# 
# df_venn_graph_SH <- df_venn_graph_S[which(df_venn_graph_S$SH==1 & df_venn_graph_S$SL==0),];
# df_venn_graph_SH$Group <- "High feed efficiency";df_venn_graph_SH <- df_venn_graph_SH[,c("seqId",colnames(ps3@tax_table),"Group")]
# df_venn_graph_LH <- df_venn_graph_S[which(df_venn_graph_S$SH==0 & df_venn_graph_S$LH==1),]
# df_venn_graph_LH$Group <- "Low feed efficiency";df_venn_graph_LH <- df_venn_graph_LH[,c("seqId",colnames(ps3@tax_table),"Group")]
# 
# df_venn_graph_S <- rbind(df_venn_graph_SH,df_venn_graph_LH)
# 
# write.csv(df_venn_graph_S,paste0(out_dir,"/","csv/","unique_per_group_solid.csv"),quote = F,row.names = F)
# 
# df_venn_graph_S$Family <-droplevels(df_venn_graph_S$Family)
# 
# df_venn_graph_melt_S <- melt(table(df_venn_graph_S$Family,df_venn_graph_S$Group))
# 
# 
# df_venn_graph_values_S <- ggplot(data=df_venn_graph_melt_S,aes(x=reorder(Var1, value),y=value,fill=Var2))+
#   geom_bar(stat="identity",show.legend = F)+
#   xlab("Family")+
#   ylab("Count")+
#   # ylim(0,0.1)+
#   scale_fill_manual("legend", values = c("red","blue"))+
#   scale_y_discrete(breaks =0:35)+
#   theme(text=element_text(size=60),
#         legend.text=element_text(size=60),
#         axis.text.x  = element_text(size=60,colour="black",angle=90),
#         axis.text.y  = element_text(size=60,colour="black"))+
#   coord_flip()
# #df_venn_graph_values_HE
# 
# ##df_venn_graph_values_S
# ggsave(paste0(graph_dir,"/","df_venn_graph_values_solid",".pdf"),df_venn_graph_values_S,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)
# 
# #########High liquid Vs Low liquid
# 
# 
# df_venn_graph_L <- df_venn
# df_venn_graph_L$seqId <- as.character(df_venn_graph_L$seqId)
# #df_venn_graph_L <- merge(df_venn_graph_L,ps3@tax_table,by.x="seqId",by.y="row.names")
# df_venn_graph_L <- left_join(df_venn_graph_L,taxa_df,by="seqId")
# 
# df_venn_graph_LH <- df_venn_graph_L[which(df_venn_graph_L$LH==1 & df_venn_graph_L$LL==0),];
# df_venn_graph_LH$Group <- "High feed efficiency";df_venn_graph_LH <- df_venn_graph_LH[,c("seqId",colnames(ps3@tax_table),"Group")]
# df_venn_graph_LL <- df_venn_graph_L[which(df_venn_graph_L$LH==0 & df_venn_graph_L$LL==1),]
# df_venn_graph_LL$Group <- "Low feed efficiency";df_venn_graph_LL <- df_venn_graph_LL[,c("seqId",colnames(ps3@tax_table),"Group")]
# 
# df_venn_graph_L <- rbind(df_venn_graph_LH,df_venn_graph_LL)
# 
# write.csv(df_venn_graph_L,paste0(out_dir,"/","csv/","unique_per_group_liquid.csv"),quote = F,row.names = F)
# 
# df_venn_graph_L$Family <-droplevels(df_venn_graph_L$Family)
# 
# df_venn_graph_melt_L <- melt(table(df_venn_graph_L$Family,df_venn_graph_L$Group))
# 
# 
# df_venn_graph_values_L <- ggplot(data=df_venn_graph_melt_L,aes(x=reorder(Var1, value),y=value,fill=Var2))+
#   geom_bar(stat="identity",show.legend = F)+
#   xlab("Family")+
#   ylab("Species count")+
#   # ylim(0,0.1)+
#   scale_fill_manual("legend", values = c("red","blue"))+
#   theme(text=element_text(size=60),
#         legend.text=element_text(size=60),
#         axis.text.x  = element_text(size=60,colour="black",angle=90),
#         axis.text.y  = element_text(size=60,colour="black"))+
#   coord_flip()
# #df_venn_graph_values_HE
# 
# ##df_venn_graph_values_S
# ggsave(paste0(graph_dir,"/","df_venn_graph_values_liquid",".pdf"),df_venn_graph_values_L,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)
# 
# 
# 
# 
# comb <- ggarrange(df_venn_graph_values_HE, df_venn_graph_values_S, df_venn_graph_values_L,# + rremove("x.text"), 
#           labels = c("A", "B", "C"),
#           ncol = 3, nrow = 1)
# 
# 
# 
# ggsave(paste0(graph_dir,"/","df_venn_graph_combined",".pdf"),comb,dpi=400,width =900,height=400,units = "cm",scale=1.2,limitsize = FALSE)
# 
# 
# 
# 
# 
# 
# # 
# # p<-p+scale_fill_manual(values=colours[1:(N+1)])
# # p<-p+theme_bw()+ylab("Relative abundances")
# # p<-p+ scale_y_continuous(expand = c(0,0))+
# #   theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines"))
# # p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=60),
# #            text=element_text(size=60),
# #            axis.text.y  = element_text(size=60,colour="black"))+
# #   xlab("Sample")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# grid.newpage()
# draw.pairwise.venn(502, 508, 465, category = c("High feed efficiency", "Low feed efficiency"), 
#                    lty = rep("blank",2), fill = c("blue", "red"), alpha = rep(0.5, 2), 
#                    cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = FALSE)
# 
# 
# ###########
# grid.newpage()
# veen_four <- draw.quad.venn(
#   area1 = nrow(subset(df_venn, SH == 1)), 
#   area2 = nrow(subset(df_venn, LH ==1)),
#   area3 = nrow(subset(df_venn, SL == 1)),
#   area4 =nrow(subset(df_venn, LL ==1)),
#   n12 = nrow(subset(df_venn, SH == 1 & LH ==1)), 
#   n13 = nrow(subset(df_venn, SH == 1 & SL ==1)), 
#   n14 = nrow(subset(df_venn, SH == 1 & LL ==1)), 
#   n23 = nrow(subset(df_venn, LH == 1 & SL ==1)), 
#   n24 = nrow(subset(df_venn, LH == 1 & LL ==1)),
#   n34 = nrow(subset(df_venn, SL == 1 & LL ==1)), 
#   n123 = nrow(subset(df_venn, SH == 1 & LH ==1 & SL ==1)), 
#   n124 = nrow(subset(df_venn, SH == 1 & LH ==1 & LL ==1)), 
#   n134 = nrow(subset(df_venn, SH == 1 & SL ==1 & LL ==1)), 
#   n234 = nrow(subset(df_venn, LH == 1 & SL ==1 & LL ==1)), 
#   n1234 = nrow(subset(df_venn, SH == 1 & LH ==1 & SL ==1 & LL==1)),
#   category = c("High feed efficiency (Solid)", "High feed efficiency (Liquid)",
#                "Low feed efficiency (Solid)","Low feed efficiency (Liquid)"),
#   lty = "blank", 
#   cex=6,
#   cat.cex = 6,
#   fill = c("darkblue", "green", "orange","red"))
# 
# ggsave(paste0(graph_dir,"/","veen_four",".pdf"),veen_four,dpi=300,width =160,height=90,units = "cm",scale=1.2,limitsize = FALSE)
# 
# 
# df_venn
# ######################
# prop <- data.frame(
#   area1 = nrow(subset(df_venn, SH == 1)), ##
#   area2 = nrow(subset(df_venn, LH ==1)),##
#   area3 = nrow(subset(df_venn, SL == 1)),## 
#   area4 =nrow(subset(df_venn, LL ==1)),##
#   n12 = nrow(subset(df_venn, SH == 1 & LH ==1)),##
#   n13 = nrow(subset(df_venn, SH == 1 & SL ==1)), ##
#   n14 = nrow(subset(df_venn, SH == 1 & LL ==1)), ##
#   n23 = nrow(subset(df_venn, LH == 1 & SL ==1)), ##
#   n24 = nrow(subset(df_venn, LH == 1 & LL ==1)), ##
#   n34 = nrow(subset(df_venn, SL == 1 & LL ==1)), #
#   n123 = nrow(subset(df_venn, SH == 1 & LH ==1 & SL ==1)),# 
#   n124 = nrow(subset(df_venn, SH == 1 & LH ==1 & LL ==1)), #
#   n134 = nrow(subset(df_venn, SH == 1 & SL ==1 & LL ==1)), #
#   n234 = nrow(subset(df_venn, LH == 1 & SL ==1 & LL ==1)), 
#   n1234 = nrow(subset(df_venn, SH == 1 & LH ==1 & SL ==1 & LL==1))
# )
# 
# 
# 
# 
# 
# prop2 <- prop
# 
# colnames(prop2) <- c(
#   "Solid_HE",#
#   "Liquid_HE",#
#   "Solid_LE",#
#   "Liquid_LE",#
#   "SH_LH",#
#   "SH_SL",#
#   "SH_LL",#
#   "LH_SL",#
#   "LH_LL",#
#   "SL_LL",#
#   "SH_LH_SL",
#   "SH_LH_LL",
#   "SH_SL_LL",
#   "LH_SL_LL",
#   "SH_LH_SL_LL"
# )
# 
# 
# prop2$HE_UNIQUE <- nrow(subset(df_venn, SH == 1 & LH ==1 & SL !=1 & LL!=1))# prop2$LE_UNIQUE <- nrow(subset(df_venn, SH != 1 & LH !=1 & SL ==1 & LL==1))
# prop2$SOL_UNIQUE <- nrow(subset(df_venn, SH == 1 & LH !=1 & SL ==1 & LL!=1))
# prop2$LIQ_UNIQUE <- nrow(subset(df_venn, SH != 1 & LH ==1 & SL !=1 & LL==1))

# 
# prop2$HE_SOL_UNIQUE <- nrow(subset(df_venn, SH == 1 & LH !=1 & SL !=1 & LL!=1))
# prop2$HE_LIQ_UNIQUE <- nrow(subset(df_venn, SH != 1 & LH ==1 & SL !=1 & LL!=1))
# 
# prop2$LE_SOL_UNIQUE <- nrow(subset(df_venn, SH == 1 & LH !=1 & SL ==1 & LL!=1))
# prop2$LE_LIQ_UNIQUE <- nrow(subset(df_venn, SH != 1 & LH !=1 & SL !=1 & LL==1))
# 
# 
# prop2[2,] <- (prop2[1,]/nrow(df_venn))*100
# 
# 
# prop2$ITEM <- NA
# prop2$ITEM[[1]] <- "COUNTS"
# prop2$ITEM[[2]] <- "PERC"
# 
# prop2 <- prop2[,c(ncol(prop2),1:(ncol(prop2)-1))]
# prop2 <-t(prop2)
# colnames(prop2) <- prop2[1,]
# prop2 <- prop2[-1,]
# write.csv(prop2,paste0(out_dir,"/","csv/","Venn_areas_group.csv"),quote = F,row.names = T)
# 
# # 
# # 
# # 
# # df_venn2 <- df_venn
# # 
# # for(i in 1:nrow(df_venn2)){
# #   ###SOLID
# #   if(df_venn2$SH[[i]]==1 & df_venn2$SL[[i]]==1){
# #     df_venn2$UNIQ_SOLID[[i]]<- "Shared"
# #   } else if(df_venn2$SH[[i]]==1 & df_venn2$SL[[i]]==0){
# #     df_venn2$UNIQ_SOLID[[i]] <- "High feed efficiency"
# #   } else if(df_venn2$SH[[i]]==0 & df_venn2$SL[[i]]==1){
# #     df_venn2$UNIQ_SOLID[[i]] <- "Low feed efficiency"
# #   } else if(df_venn2$SH[[i]]==0 & df_venn2$SL[[i]]==0){
# #     df_venn2$UNIQ_SOLID[[i]] <- "Not sp"
# #   }
# #   ###LIQUID 
# #   if(df_venn2$LH[[i]]==1 & df_venn2$LL[[i]]==1){
# #     df_venn2$UNIQ_SOLID[[i]]<- "Shared"
# #   } else if(df_venn2$LH[[i]]==1 & df_venn2$LL[[i]]==0){
# #     df_venn2$UNIQ_SOLID[[i]] <- "High feed efficiency"
# #   } else if(df_venn2$LH[[i]]==0 & df_venn2$LL[[i]]==1){
# #     df_venn2$UNIQ_SOLID[[i]] <- "Low feed efficiency"
# #   } else if(df_venn2$LH[[i]]==0 & df_venn2$LL[[i]]==0){
# #     df_venn2$UNIQ_SOLID[[i]] <- "Not sp"
# #   }
# #   ###"SOLID_Vs_LIQUID
# #   if(df_venn2$LH[[i]]==1 & df_venn2$LL[[i]]==1){
# #     df_venn2$UNIQ_SOLID[[i]]<- "Shared"
# #   } else if(df_venn2$LH[[i]]==1 & df_venn2$LL[[i]]==0){
# #     df_venn2$UNIQ_SOLID[[i]] <- "High feed efficiency"
# #   } else if(df_venn2$LH[[i]]==0 & df_venn2$LL[[i]]==1){
# #     df_venn2$UNIQ_SOLID[[i]] <- "Low feed efficiency"
# #   } else if(df_venn2$LH[[i]]==0 & df_venn2$LL[[i]]==0){
# #     df_venn2$UNIQ_SOLID[[i]] <- "Not sp"
# #   }
# #   
# # };rm(i)
# # 
# # df_venn2
# 
