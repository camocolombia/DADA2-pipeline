
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
#install.packages("devtools")
library("devtools")
library("phytools")
library("VennDiagram")
#install_github("JiyuanHu/massMap")

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


#install.packages('VennDiagram')

ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object",".RDS"))

abund <- as.data.frame(t(abundances(ps)))


metadata <- as(sample_data(ps), "data.frame")

g <- global(ps)
g
sample_data(ps)$Shannon <- global(ps, index = "shannon")[,1]
sample_data(ps)$evenness_pielou <- global(ps, index = "pielou")[,1]
sample_data(ps)$Simpson <- global(ps, index = "simpson")[,1]
s1bp <- ggplot(data=sample_data(ps),aes(x=Phase,y=Shannon,fill=Treatment))+
  geom_boxplot(position = "dodge",outlier.colour = NA,
               outlier.fill=NA,outlier.alpha=1,
               outlier.size =NA,na.rm = T,show.legend = F) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Shannon index")+ xlab("")+
  ggtitle("")+
  scale_fill_manual("Groups", values = c("blue","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=90,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=80,colour="black")#,
        #legend.title=element_blank(),legend.position="none"
  )+
  coord_flip()
s1bp
ggsave(paste0(graph_dir,"/","Shannon_group_phase",".pdf"),s1bp,dpi=300,width =100,height=80,units = "cm",scale=1.2,limitsize = FALSE)

###Treatment
s1bp_T <- ggplot(data=sample_data(ps),aes(x=Treatment,y=Shannon,fill=Treatment))+
  geom_boxplot(position = "dodge",outlier.colour = NA,
               outlier.fill=NA,outlier.alpha=1,
               outlier.size =NA,na.rm = T,show.legend = F) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Shannon index")+ xlab("")+
  ggtitle("")+
  scale_fill_manual("Groups", values = c("blue","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=90,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=80,colour="black")#,
        #legend.title=element_blank(),legend.position="none"
  )+
  coord_flip()
s1bp_T
ggsave(paste0(graph_dir,"/","Shannon_group",".pdf"),s1bp_T,dpi=300,width =100,height=80,units = "cm",scale=1.2,limitsize = FALSE)



###
############ DOMINANCCE

s1bp_SI <- ggplot(data=sample_data(ps),aes(x=Phase,y=Simpson,fill=Treatment))+
  geom_boxplot(position = "dodge",outlier.colour = NA,
               outlier.fill=NA,outlier.alpha=1,
               outlier.size =NA,na.rm = T,show.legend = F) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Simpson index")+ xlab("")+
  ggtitle("")+
  scale_fill_manual("Groups", values = c("blue","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=90,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=80,colour="black")#,
        #legend.title=element_blank(),legend.position="none"
  )+
  coord_flip()
s1bp_SI
ggsave(paste0(graph_dir,"/","Simpson_group_phase",".pdf"),s1bp_SI,dpi=300,width =100,height=80,units = "cm",scale=1.2,limitsize = FALSE)

####
s1bp_SI_T <- ggplot(data=sample_data(ps),aes(x=Treatment,y=Simpson,fill=Treatment))+
  geom_boxplot(position = "dodge",outlier.colour = NA,
               outlier.fill=NA,outlier.alpha=1,
               outlier.size =NA,na.rm = T,show.legend = F) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Simpson index")+ xlab("")+
  ggtitle("")+
  scale_fill_manual("Groups", values = c("blue","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),axis.text.x  = element_text(size=90,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=80,colour="black")#,
        #legend.title=element_blank(),legend.position="none"
  )+
  coord_flip()
s1bp_SI_T
ggsave(paste0(graph_dir,"/","Simpsons_group",".pdf"),s1bp_SI_T,dpi=300,width =100,height=80,units = "cm",scale=1.2,limitsize = FALSE)




# cor.test(sample_data(ps)$evenness_pielou, sample_data(ps)$FCR)
require(ggpubr)
sample_data(ps)$groups <- paste0(sample_data(ps)$Treatment," (",sample_data(ps)$Phase,")")
sample_data(ps)$groups <- factor(sample_data(ps)$groups,levels=
                                   c("High feed efficiency (solid)", "High feed efficiency (liquid)",
                                     "Low feed efficiency (solid)","Low feed efficiency (liquid)")                 
)
sp <- ggpubr::ggscatterhist(
  data= sample_data(ps), x = "Shannon", y = "FCR",
  color = "groups",size = 3, alpha = 0.6,
  palette = c("darkblue", "green", "orange","red"),#c("red","blue"),
  margin.params = list(fill = "groups", color = "black", size = 0.2),
  ggtheme = theme_bw()
)
sp

#plot_richness(ps,measures = c("Shannon","Simpson"),color = "Treatment",shape = "Phase")



library(psych)
# describe(subset(sample_data(ps),Treatment=="High feed efficiency")$Simpson,
#          type=2)$mean   
###########################################################################
###########################################################################
###########################################################################


wilcox_list_DIV <- data.frame(
  ###########################################################################
  #FEED EFFICIENCY SIMPSON
SIMP_HE = describe(subset(sample_data(ps),Treatment=="High feed efficiency")$Simpson,
                   type=2)$mean,
SIMP_HE_SE = describe(subset(sample_data(ps),Treatment=="High feed efficiency")$Simpson,
                   type=2)$sd,
SIMP_LE = describe(subset(sample_data(ps),Treatment=="Low feed efficiency")$Simpson,
                   type=2)$mean,
SIMP_LE_SE = describe(subset(sample_data(ps),Treatment=="Low feed efficiency")$Simpson,
                   type=2)$sd,
SIMP_EFF_pvalue = wilcox.test(subset(sample_data(ps),Treatment=="High feed efficiency")$Simpson,
                       subset(sample_data(ps),Treatment=="Low feed efficiency")$Simpson)$p.value,
###########################################################################
#FEED EFFICIENCY SHANNON
SHANN_HE = describe(subset(sample_data(ps),Treatment=="High feed efficiency")$Shannon,
                   type=2)$mean,
SHANN_HE_SD = describe(subset(sample_data(ps),Treatment=="High feed efficiency")$Shannon,
                      type=2)$sd,
SHANN_LE = describe(subset(sample_data(ps),Treatment=="Low feed efficiency")$Shannon,
                   type=2)$mean,
SHANN_LE_SD = describe(subset(sample_data(ps),Treatment=="Low feed efficiency")$Shannon,
                      type=2)$sd,
SHANN_EFF_pvalue = wilcox.test(subset(sample_data(ps),Treatment=="High feed efficiency")$Shannon,
                       subset(sample_data(ps),Treatment=="Low feed efficiency")$Shannon)$p.value,
###########################################################################
#PHASE SIMPSON
SIMP_SO = describe(subset(sample_data(ps),Phase=="solid")$Simpson,
                    type=2)$mean,
SIMP_SO_SD = describe(subset(sample_data(ps),Phase=="solid")$Simpson,
                       type=2)$sd,
SIMP_LI = describe(subset(sample_data(ps),Phase=="liquid")$Simpson,
                    type=2)$mean,
SIMP_LI_SD = describe(subset(sample_data(ps),Phase=="liquid")$Simpson,
                       type=2)$sd,
SIMP_PHASE_pvalue = wilcox.test(subset(sample_data(ps),Phase=="liquid")$Simpson,
                                 subset(sample_data(ps),Phase=="solid")$Simpson)$p.value,
###########################################################################
#PHASE SHANNON
SHANN_SO = describe(subset(sample_data(ps),Phase=="solid")$Shannon,
                    type=2)$mean,
SHANN_SO_SD = describe(subset(sample_data(ps),Phase=="solid")$Shannon,
                       type=2)$sd,
SHANN_LI = describe(subset(sample_data(ps),Phase=="liquid")$Shannon,
                    type=2)$mean,
SHANN_LI_SD = describe(subset(sample_data(ps),Phase=="liquid")$Shannon,
                       type=2)$sd,
SHANN_PHASE_pvalue = wilcox.test(subset(sample_data(ps),Phase=="liquid")$Shannon,
                               subset(sample_data(ps),Phase=="solid")$Shannon)$p.value,
###########################################################################
#PHASE SIMPSON HFE
SIMP_SO_HE = describe(subset(sample_data(ps),Phase=="solid" & Treatment=="High feed efficiency")$Simpson,
                   type=2)$mean,
SIMP_SO_HE_SD = describe(subset(sample_data(ps),Phase=="solid" & Treatment=="High feed efficiency")$Simpson,
                      type=2)$sd,
SIMP_LI_HE = describe(subset(sample_data(ps),Phase=="liquid"& Treatment=="High feed efficiency")$Simpson,
                   type=2)$mean,
SIMP_LI_HE_SD = describe(subset(sample_data(ps),Phase=="liquid" & Treatment=="High feed efficiency")$Simpson,
                      type=2)$sd,
SIMP_PHASE_HE_pvalue = wilcox.test(subset(sample_data(ps),Phase=="liquid"  & Treatment=="High feed efficiency")$Simpson,
                                subset(sample_data(ps),Phase=="solid" & Treatment=="High feed efficiency")$Simpson)$p.value,
###########################################################################
#PHASE SHANNON HFE
SHANN_O_HE = describe(subset(sample_data(ps),Phase=="solid" & Treatment=="High feed efficiency")$Shannon,
                      type=2)$mean,
SHANN_SO_HE_SD = describe(subset(sample_data(ps),Phase=="solid" & Treatment=="High feed efficiency")$Shannon,
                         type=2)$sd,
SHANN_LI_HE = describe(subset(sample_data(ps),Phase=="liquid"& Treatment=="High feed efficiency")$Shannon,
                      type=2)$mean,
SHANN_LI_HE_SD = describe(subset(sample_data(ps),Phase=="liquid" & Treatment=="High feed efficiency")$Shannon,
                         type=2)$sd,
SHANN_PHASE_HE_pvalue = wilcox.test(subset(sample_data(ps),Phase=="liquid"  & Treatment=="High feed efficiency")$Shannon,
                                   subset(sample_data(ps),Phase=="solid" & Treatment=="High feed efficiency")$Shannon)$p.value,
###########################################################################
#PHASE SIMPSON LFE
SIMP_SO_LE = describe(subset(sample_data(ps),Phase=="solid" & Treatment=="Low feed efficiency")$Simpson,
                      type=2)$mean,
SIMP_SO_LE_SD = describe(subset(sample_data(ps),Phase=="solid" & Treatment=="Low feed efficiency")$Simpson,
                         type=2)$sd,
SIMP_LI_LE = describe(subset(sample_data(ps),Phase=="liquid"& Treatment=="Low feed efficiency")$Simpson,
                      type=2)$mean,
SIMP_LI_LE_SD = describe(subset(sample_data(ps),Phase=="liquid" & Treatment=="Low feed efficiency")$Simpson,
                         type=2)$sd,
SIMP_PHASE_LE_pvalue = wilcox.test(subset(sample_data(ps),Phase=="liquid"  & Treatment=="Low feed efficiency")$Simpson,
                                   subset(sample_data(ps),Phase=="solid" & Treatment=="Low feed efficiency")$Simpson)$p.value,
##########################################################################
#PHASE SHANNON LFE
SHANN_SO_LE = describe(subset(sample_data(ps),Phase=="solid" & Treatment=="Low feed efficiency")$Shannon,
                      type=2)$mean,
SHANN_SO_LE_SD = describe(subset(sample_data(ps),Phase=="solid" & Treatment=="Low feed efficiency")$Shannon,
                         type=2)$sd,
SHANN_LI_LE = describe(subset(sample_data(ps),Phase=="liquid"& Treatment=="Low feed efficiency")$Shannon,
                      type=2)$mean,
SHANN_LI_LE_SD = describe(subset(sample_data(ps),Phase=="liquid" & Treatment=="Low feed efficiency")$Shannon,
                         type=2)$sd,
SHANN_PHASE_LE_pvalue = wilcox.test(subset(sample_data(ps),Phase=="liquid"  & Treatment=="Low feed efficiency")$Shannon,
                                   subset(sample_data(ps),Phase=="solid" & Treatment=="Low feed efficiency")$Shannon)$p.value,
###########################################################################
#COMPARING LIQUID FRACTIONS
SIMP_LIQs_pvalue = wilcox.test(subset(sample_data(ps),Phase=="liquid"  & Treatment=="Low feed efficiency")$Simpson,
                                   subset(sample_data(ps),Phase=="liquid" & Treatment=="High feed efficiency")$Simpson)$p.value,
#####
SHANN_LIQs_pvalue = wilcox.test(subset(sample_data(ps),Phase=="liquid"  & Treatment=="Low feed efficiency")$Shannon,
                               subset(sample_data(ps),Phase=="liquid" & Treatment=="High feed efficiency")$Shannon)$p.value,

###########################################################################
#COMPARING SOLID FRACTIONS
SIMP_SOLs_pvalue = wilcox.test(subset(sample_data(ps),Phase=="solid"  & Treatment=="Low feed efficiency")$Simpson,
                               subset(sample_data(ps),Phase=="solid" & Treatment=="High feed efficiency")$Simpson)$p.value,
#####
SHANN_SOLs_pvalue = wilcox.test(subset(sample_data(ps),Phase=="solid"  & Treatment=="Low feed efficiency")$Shannon,
                                subset(sample_data(ps),Phase=="solid" & Treatment=="High feed efficiency")$Shannon)$p.value,

###########################################################################
###########################################################################




SIMP_EFF_FDR = NA,
SHANN_EFF_FDR = NA,
SIMP_PHASE_FDR = NA,
SHANN_PHASE_FDR = NA,
SIMP_PHASE_HE_FDR = NA,
SHANN_PHASE_HE_FDR = NA,
SIMP_PHASE_LE_FDR = NA,
SHANN_PHASE_LE_FDR = NA,
SIMP_LIQs_FDR = NA,
SHANN_LIQs_FDR = NA,
SIMP_SOLs_FDR = NA,
SHANN_SOLs_FDR = NA
# 
# SHAN_EFF = wilcox.test(subset(sample_data(ps),Treatment=="High feed efficiency")$Shannon,
#             subset(sample_data(ps),Treatment=="Low feed efficiency")$Shannon)$p.value,
# 
# SIMP_PHASE = wilcox.test(subset(sample_data(ps),Phase=="liquid")$Simpson,
#             subset(sample_data(ps),Phase=="solid")$Simpson)$p.value,
# 
# SHAN_PHASE = wilcox.test(subset(sample_data(ps),Phase=="liquid")$Shannon,
#             subset(sample_data(ps),Phase=="solid")$Shannon)$p.value
)


wilcox_list_DIV$SIMP_EFF_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                               wilcox_list_DIV$SHANN_EFF_pvalue,
                                               wilcox_list_DIV$SIMP_PHASE_pvalue,
                                               wilcox_list_DIV$SHANN_PHASE_pvalue,
                                               wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                               wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                               wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                               wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                               wilcox_list_DIV$SIMP_LIQs_pvalue,
                                               wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                               wilcox_list_DIV$SIMP_SOLs_pvalue,
                                               wilcox_list_DIV$SHANN_SOLs_pvalue       
), method =  "BH",n = 12)[1]

wilcox_list_DIV$SHANN_EFF_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                wilcox_list_DIV$SHANN_EFF_pvalue,
                                                wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                wilcox_list_DIV$SHANN_SOLs_pvalue    
                                                ), method =  "BH",n = 12)[2]
                                            
                                            
wilcox_list_DIV$SIMP_PHASE_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                 wilcox_list_DIV$SHANN_EFF_pvalue,
                                                 wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                 wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                 wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                 wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                 wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                 wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                 wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                 wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                 wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                 wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[3]

wilcox_list_DIV$SHANN_PHASE_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                  wilcox_list_DIV$SHANN_EFF_pvalue,
                                                  wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                  wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                  wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                  wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                  wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                  wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                  wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                  wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                  wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                  wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[4]                                            
                                            
wilcox_list_DIV$SIMP_PHASE_HE_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                  wilcox_list_DIV$SHANN_EFF_pvalue,
                                                  wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                  wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                  wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                  wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                  wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                  wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                  wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                  wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                  wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                  wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[5]      
                                  
wilcox_list_DIV$SHANN_PHASE_HE_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                    wilcox_list_DIV$SHANN_EFF_pvalue,
                                                    wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                    wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                    wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                    wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                    wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                    wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                    wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                    wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                    wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                    wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[6]      

wilcox_list_DIV$SIMP_PHASE_LE_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                     wilcox_list_DIV$SHANN_EFF_pvalue,
                                                     wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                     wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                     wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                     wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                     wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                     wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                     wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                     wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                     wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                     wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[7]      

wilcox_list_DIV$SHANN_PHASE_LE_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                    wilcox_list_DIV$SHANN_EFF_pvalue,
                                                    wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                    wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                    wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                    wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                    wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                    wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                    wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                    wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                    wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                    wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[8]      

wilcox_list_DIV$SIMP_LIQs_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                    wilcox_list_DIV$SHANN_EFF_pvalue,
                                                    wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                    wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                    wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                    wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                    wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                    wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                    wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                    wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                    wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                    wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[9]      

wilcox_list_DIV$SHANN_LIQs_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                     wilcox_list_DIV$SHANN_EFF_pvalue,
                                                     wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                     wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                     wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                     wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                     wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                     wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                     wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                     wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                     wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                     wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[10]      

wilcox_list_DIV$SIMP_SOLs_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                wilcox_list_DIV$SHANN_EFF_pvalue,
                                                wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[11]      

wilcox_list_DIV$SHANN_SOLs_FDR <- p.adjust(p = c(wilcox_list_DIV$SIMP_EFF_pvalue,
                                                 wilcox_list_DIV$SHANN_EFF_pvalue,
                                                 wilcox_list_DIV$SIMP_PHASE_pvalue,
                                                 wilcox_list_DIV$SHANN_PHASE_pvalue,
                                                 wilcox_list_DIV$SIMP_PHASE_HE_pvalue,
                                                 wilcox_list_DIV$SHANN_PHASE_HE_pvalue,
                                                 wilcox_list_DIV$SIMP_PHASE_LE_pvalue,
                                                 wilcox_list_DIV$SHANN_PHASE_LE_pvalue,                                               
                                                 wilcox_list_DIV$SIMP_LIQs_pvalue,
                                                 wilcox_list_DIV$SHANN_LIQs_pvalue,                                               
                                                 wilcox_list_DIV$SIMP_SOLs_pvalue,
                                                 wilcox_list_DIV$SHANN_SOLs_pvalue    
), method =  "BH",n = 12)[12]      


write.csv(wilcox_list_DIV,paste0(out_dir,"/","csv/","Wilcoxon_alpha.csv"),quote = F,row.names = F)

# wilcox_list_DIV[2,] <- NA
# 
# wilcox_list_DIV[2,] <-  p.adjust(p = wilcox_list_DIV[1,],method =  "BH",n = length(wilcox_list_DIV[1,]))
# 

#install.packages("devtools")
library(devtools)
#install_github("easyGgplot2", "kassambara")
library(easyGgplot2)
ab <- ggplot2.multiplot(s1bp,s1bp_SI, plotlist=NULL, cols=2)

ggsave(paste0(graph_dir,"/","Alpha",".pdf"),ggplot2.multiplot(s1bp,s1bp_SI, plotlist=NULL, cols=2),dpi=300,width =200,height=100,units = "cm",scale=1.2,limitsize = FALSE)
