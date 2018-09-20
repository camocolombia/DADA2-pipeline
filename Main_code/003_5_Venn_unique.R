
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
# library("DESeq2")
# library("BSDA")
library("adegenet")
#install.packages("devtools")
library("devtools")
library("phytools")
library("massMap")
library("structSSI")
library("microbiomeViz");
library("stringr")
library("StructFDR")
library("gtools")
library("VennDiagram")
library("ggpubr")
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
FDR_list <- readRDS(paste0(out_dir,"/","csv/","FDR_list.RDS"))

ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))
ps@sam_data$Treatment <- factor(ps@sam_data$Treatment,levels = c("High feed efficiency","Low feed efficiency"))
ps3 <-  readRDS(paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))

##################################################
df <- FDR_list$BH_prop
##################################################
##################################################
#Defining dataframe for venn diagrams

df_venn <- df[,c("seqId","SH","SL","LH","LL","HE","LE")]
for(i in 2:7){
  df_venn[,i][which(df_venn[,i]>0)] <- 1
}

df_venn$UNIQ_SOLID <- rowSums(df_venn[,c("SH","SL")]) ; df_venn$UNIQ_SOLID[which(df_venn$UNIQ_SOLID <2)] <- 0;df_venn$UNIQ_SOLID[which(df_venn$UNIQ_SOLID==2)] <- 1
df_venn$UNIQ_LIQUID <- rowSums(df_venn[,c("LH","LL")]) ; df_venn$UNIQ_LIQUID[which(df_venn$UNIQ_LIQUID <2)] <- 0;df_venn$UNIQ_LIQUID[which(df_venn$UNIQ_LIQUID==2)] <- 1

df_venn$SOLID <- NA;df_venn$SOLID <- rowSums(df_venn[,c("SH","SL")]); df_venn$SOLID[which(df_venn$SOLID <2)] <- 0;df_venn$SOLID[which(df_venn$SOLID==2)] <- 1
df_venn$LIQUID <- NA;df_venn$LIQUID <- rowSums(df_venn[,c("LH","LL")]);df_venn$LIQUID[which(df_venn$LIQUID <2)] <- 0;df_venn$LIQUID[which(df_venn$LIQUID==2)] <- 1
df_venn$SOLID_Vs_LIQUID <- NA;df_venn$SOLID_Vs_LIQUID <- rowSums(df_venn[,c("UNIQ_SOLID","UNIQ_LIQUID")]);df_venn$SOLID_Vs_LIQUID[which(df_venn$SOLID_Vs_LIQUID <2)] <- 0;df_venn$SOLID_Vs_LIQUID[which(df_venn$SOLID_Vs_LIQUID==2)] <- 1
df_venn$HIGH <- NA;df_venn$HIGH <- rowSums(df_venn[,c("SH","LH")]);df_venn$HIGH[which(df_venn$HIGH <2)] <- 0;df_venn$HIGH[which(df_venn$HIGH==2)] <- 1
df_venn$LOW <- NA;df_venn$LOW <- rowSums(df_venn[,c("SL","LL")]);df_venn$LOW[which(df_venn$LOW <2)] <- 0;df_venn$LOW[which(df_venn$LOW==2)] <- 1
df_venn$EFF <- NA;df_venn$EFF <- rowSums(df_venn[,c("HE","LE")]);df_venn$EFF[which(df_venn$EFF <2)] <- 0;df_venn$EFF[which(df_venn$EFF==2)] <- 1

#colSums(df_venn[,c(2:15)])

write.csv(df_venn,paste0(out_dir,"/","csv/","df_venn.csv"),quote = F,row.names = F)
write.csv(df,paste0(out_dir,"/","csv/","final_table.csv"),quote = F,row.names = F)

##########High Vs Low


df_venn_graph <- df_venn
df_venn_graph <- merge(df_venn_graph,ps3@tax_table,by.x="seqId",by.y="row.names")
df_venn_graph_HE <- df_venn_graph[which(df_venn_graph$HE==1 & df_venn_graph$LE==0),];
df_venn_graph_HE$Group <- "High feed efficiency";df_venn_graph_HE <- df_venn_graph_HE[,c("seqId",colnames(ps3@tax_table),"Group")]
df_venn_graph_LE <- df_venn_graph[which(df_venn_graph$HE==0 & df_venn_graph$LE==1),]
df_venn_graph_LE$Group <- "Low feed efficiency";df_venn_graph_LE <- df_venn_graph_LE[,c("seqId",colnames(ps3@tax_table),"Group")]

df_venn_graph <- rbind(df_venn_graph_HE,df_venn_graph_LE)

write.csv(df_venn_graph,paste0(out_dir,"/","csv/","unique_per_group.csv"),quote = F,row.names = F)

df_venn_graph$Family <-droplevels(df_venn_graph$Family)

df_venn_graph_melt <- melt(table(df_venn_graph$Family,df_venn_graph$Group))


df_venn_graph_values_HE <- ggplot(data=df_venn_graph_melt,aes(x=reorder(Var1, value),y=value,fill=Var2))+
  geom_bar(stat="identity",show.legend = F)+
  xlab("Family")+
  ylab("Count")+
  # ylim(0,0.1)+
  scale_fill_manual("legend", values = c("red","blue"))+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black",angle=90),
        axis.text.y  = element_text(size=60,colour="black"))+
  coord_flip()
#df_venn_graph_values_HE


ggsave(paste0(graph_dir,"/","df_venn_graph_values",".pdf"),df_venn_graph_values_HE,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)


##########High solid Vs Low solid


df_venn_graph_S <- df_venn
df_venn_graph_S <- merge(df_venn_graph_S,ps3@tax_table,by.x="seqId",by.y="row.names")
df_venn_graph_SH <- df_venn_graph_S[which(df_venn_graph_S$SH==1 & df_venn_graph_S$SL==0),];
df_venn_graph_SH$Group <- "High feed efficiency";df_venn_graph_SH <- df_venn_graph_SH[,c("seqId",colnames(ps3@tax_table),"Group")]
df_venn_graph_LH <- df_venn_graph_S[which(df_venn_graph_S$SH==0 & df_venn_graph_S$LH==1),]
df_venn_graph_LH$Group <- "Low feed efficiency";df_venn_graph_LH <- df_venn_graph_LH[,c("seqId",colnames(ps3@tax_table),"Group")]

df_venn_graph_S <- rbind(df_venn_graph_SH,df_venn_graph_LH)

write.csv(df_venn_graph_S,paste0(out_dir,"/","csv/","unique_per_group_solid.csv"),quote = F,row.names = F)

df_venn_graph_S$Family <-droplevels(df_venn_graph_S$Family)

df_venn_graph_melt_S <- melt(table(df_venn_graph_S$Family,df_venn_graph_S$Group))


df_venn_graph_values_S <- ggplot(data=df_venn_graph_melt_S,aes(x=reorder(Var1, value),y=value,fill=Var2))+
  geom_bar(stat="identity",show.legend = F)+
  xlab("Family")+
  ylab("Count")+
  # ylim(0,0.1)+
  scale_fill_manual("legend", values = c("red","blue"))+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black",angle=90),
        axis.text.y  = element_text(size=60,colour="black"))+
  coord_flip()
#df_venn_graph_values_HE

##df_venn_graph_values_S
ggsave(paste0(graph_dir,"/","df_venn_graph_values_solid",".pdf"),df_venn_graph_values_S,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)

#########High liquid Vs Low liquid


df_venn_graph_L <- df_venn
df_venn_graph_L <- merge(df_venn_graph_L,ps3@tax_table,by.x="seqId",by.y="row.names")
df_venn_graph_LH <- df_venn_graph_L[which(df_venn_graph_L$LH==1 & df_venn_graph_L$LL==0),];
df_venn_graph_LH$Group <- "High feed efficiency";df_venn_graph_LH <- df_venn_graph_LH[,c("seqId",colnames(ps3@tax_table),"Group")]
df_venn_graph_LL <- df_venn_graph_L[which(df_venn_graph_L$LH==0 & df_venn_graph_L$LL==1),]
df_venn_graph_LL$Group <- "Low feed efficiency";df_venn_graph_LL <- df_venn_graph_LL[,c("seqId",colnames(ps3@tax_table),"Group")]

df_venn_graph_L <- rbind(df_venn_graph_LH,df_venn_graph_LL)

write.csv(df_venn_graph_L,paste0(out_dir,"/","csv/","unique_per_group_liquid.csv"),quote = F,row.names = F)

df_venn_graph_L$Family <-droplevels(df_venn_graph_L$Family)

df_venn_graph_melt_L <- melt(table(df_venn_graph_L$Family,df_venn_graph_L$Group))


df_venn_graph_values_L <- ggplot(data=df_venn_graph_melt_L,aes(x=reorder(Var1, value),y=value,fill=Var2))+
  geom_bar(stat="identity",show.legend = F)+
  xlab("Family")+
  ylab("Species count")+
  # ylim(0,0.1)+
  scale_fill_manual("legend", values = c("red","blue"))+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black",angle=90),
        axis.text.y  = element_text(size=60,colour="black"))+
  coord_flip()
#df_venn_graph_values_HE

##df_venn_graph_values_S
ggsave(paste0(graph_dir,"/","df_venn_graph_values_liquid",".pdf"),df_venn_graph_values_L,dpi=300,width =100,height=90,units = "cm",scale=1.2,limitsize = FALSE)




comb <- ggarrange(df_venn_graph_values_HE, df_venn_graph_values_S, df_venn_graph_values_L,# + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)



ggsave(paste0(graph_dir,"/","df_venn_graph_combined",".pdf"),comb,dpi=400,width =900,height=400,units = "cm",scale=1.2,limitsize = FALSE)






# 
# p<-p+scale_fill_manual(values=colours[1:(N+1)])
# p<-p+theme_bw()+ylab("Relative abundances")
# p<-p+ scale_y_continuous(expand = c(0,0))+
#   theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines"))
# p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=60),
#            text=element_text(size=60),
#            axis.text.y  = element_text(size=60,colour="black"))+
#   xlab("Sample")














grid.newpage()
draw.pairwise.venn(502, 508, 465, category = c("High feed efficiency", "Low feed efficiency"), 
                   lty = rep("blank",2), fill = c("blue", "red"), alpha = rep(0.5, 2), 
                   cat.pos = c(0,0), cat.dist = rep(0.025, 2), scaled = FALSE)


###########
grid.newpage()
veen_four <- draw.quad.venn(
  area1 = nrow(subset(df_venn, SH == 1)), 
  area2 = nrow(subset(df_venn, LH ==1)),
  area3 = nrow(subset(df_venn, SL == 1)),
  area4 =nrow(subset(df_venn, LL ==1)),
  n12 = nrow(subset(df_venn, SH == 1 & LH ==1)), 
  n13 = nrow(subset(df_venn, SH == 1 & SL ==1)), 
  n14 = nrow(subset(df_venn, SH == 1 & LL ==1)), 
  n23 = nrow(subset(df_venn, LH == 1 & SL ==1)), 
  n24 = nrow(subset(df_venn, LH == 1 & LL ==1)),
  n34 = nrow(subset(df_venn, SL == 1 & LL ==1)), 
  n123 = nrow(subset(df_venn, SH == 1 & LH ==1 & SL ==1)), 
  n124 = nrow(subset(df_venn, SH == 1 & LH ==1 & LL ==1)), 
  n134 = nrow(subset(df_venn, SH == 1 & SL ==1 & LL ==1)), 
  n234 = nrow(subset(df_venn, LH == 1 & SL ==1 & LL ==1)), 
  n1234 = nrow(subset(df_venn, SH == 1 & LH ==1 & SL ==1 & LL==1)),
  category = c("High feed efficiency (Solid)", "High feed efficiency (Liquid)",
               "Low feed efficiency (Solid)","Low feed efficiency (Liquid)"),
  lty = "blank", 
  cex=6,
  cat.cex = 6,
  fill = c("darkblue", "green", "orange","red"))

ggsave(paste0(graph_dir,"/","veen_four",".pdf"),veen_four,dpi=300,width =160,height=90,units = "cm",scale=1.2,limitsize = FALSE)


df_venn
######################
prop <- data.frame(
  area1 = nrow(subset(df_venn, SH == 1)), 
  area2 = nrow(subset(df_venn, LH ==1)),
  area3 = nrow(subset(df_venn, SL == 1)),
  area4 =nrow(subset(df_venn, LL ==1)),
  n12 = nrow(subset(df_venn, SH == 1 & LH ==1)), 
  n13 = nrow(subset(df_venn, SH == 1 & SL ==1)), 
  n14 = nrow(subset(df_venn, SH == 1 & LL ==1)), 
  n23 = nrow(subset(df_venn, LH == 1 & SL ==1)), 
  n24 = nrow(subset(df_venn, LH == 1 & LL ==1)),
  n34 = nrow(subset(df_venn, SL == 1 & LL ==1)), 
  n123 = nrow(subset(df_venn, SH == 1 & LH ==1 & SL ==1)), 
  n124 = nrow(subset(df_venn, SH == 1 & LH ==1 & LL ==1)), 
  n134 = nrow(subset(df_venn, SH == 1 & SL ==1 & LL ==1)), 
  n234 = nrow(subset(df_venn, LH == 1 & SL ==1 & LL ==1)), 
  n1234 = nrow(subset(df_venn, SH == 1 & LH ==1 & SL ==1 & LL==1))
)

prop[2,] <- (prop[1,]/545)*100
# 
# 
# 
# df_venn2 <- df_venn
# 
# for(i in 1:nrow(df_venn2)){
#   ###SOLID
#   if(df_venn2$SH[[i]]==1 & df_venn2$SL[[i]]==1){
#     df_venn2$UNIQ_SOLID[[i]]<- "Shared"
#   } else if(df_venn2$SH[[i]]==1 & df_venn2$SL[[i]]==0){
#     df_venn2$UNIQ_SOLID[[i]] <- "High feed efficiency"
#   } else if(df_venn2$SH[[i]]==0 & df_venn2$SL[[i]]==1){
#     df_venn2$UNIQ_SOLID[[i]] <- "Low feed efficiency"
#   } else if(df_venn2$SH[[i]]==0 & df_venn2$SL[[i]]==0){
#     df_venn2$UNIQ_SOLID[[i]] <- "Not sp"
#   }
#   ###LIQUID 
#   if(df_venn2$LH[[i]]==1 & df_venn2$LL[[i]]==1){
#     df_venn2$UNIQ_SOLID[[i]]<- "Shared"
#   } else if(df_venn2$LH[[i]]==1 & df_venn2$LL[[i]]==0){
#     df_venn2$UNIQ_SOLID[[i]] <- "High feed efficiency"
#   } else if(df_venn2$LH[[i]]==0 & df_venn2$LL[[i]]==1){
#     df_venn2$UNIQ_SOLID[[i]] <- "Low feed efficiency"
#   } else if(df_venn2$LH[[i]]==0 & df_venn2$LL[[i]]==0){
#     df_venn2$UNIQ_SOLID[[i]] <- "Not sp"
#   }
#   ###"SOLID_Vs_LIQUID
#   if(df_venn2$LH[[i]]==1 & df_venn2$LL[[i]]==1){
#     df_venn2$UNIQ_SOLID[[i]]<- "Shared"
#   } else if(df_venn2$LH[[i]]==1 & df_venn2$LL[[i]]==0){
#     df_venn2$UNIQ_SOLID[[i]] <- "High feed efficiency"
#   } else if(df_venn2$LH[[i]]==0 & df_venn2$LL[[i]]==1){
#     df_venn2$UNIQ_SOLID[[i]] <- "Low feed efficiency"
#   } else if(df_venn2$LH[[i]]==0 & df_venn2$LL[[i]]==0){
#     df_venn2$UNIQ_SOLID[[i]] <- "Not sp"
#   }
#   
# };rm(i)
# 
# df_venn2