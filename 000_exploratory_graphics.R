require(ggplot2);require(reshape)
mainDir <- "E:/Dropbox/Dropbox/Paper_PhD"
chapter <- "chapter_1"
dat_dir <- paste0(mainDir,"/",chapter,"/","data"); if(!file.exists(dat_dir)){dir.create(dat_dir)}
production <- read.csv(paste0(dat_dir,"/","production_data.csv"),header=T,sep="|")
out_dir <- paste0(mainDir,"/",chapter,"/","outcomes"); if(!file.exists(out_dir)){dir.create(out_dir)}
graph_dir <- paste0(out_dir,"/","graphics"); if(!file.exists(graph_dir)){dir.create(graph_dir)}
csv_dir <- paste0(out_dir,"/","csv"); if(!file.exists(csv_dir)){dir.create(csv_dir)}
production$Treatment <- gsub("bottom","Bottom",production$Treatment); production$Treatment <- gsub("top","Top",production$Treatment)

p1a <- ggplot(data=production,aes(x=Treatment,y=FCE))+
  geom_boxplot(aes(fill=Treatment)) +
  stat_boxplot(geom ='errorbar') +
  stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, geom="point")+
  guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Average daily gain  (kg DM)")+ xlab("Treatment")+
  ggtitle("")+
  scale_fill_manual("legend", values = c("blue","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=40),axis.text.x  = element_text(size=40,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=40,colour="black"),
        legend.title=element_blank(),legend.position="none")
p1a

ggsave(paste0(graph_dir,"/","ADG",".pdf"),p1a,dpi=600,width =90,height=34.395,units = "cm",scale=1.2,limitsize = FALSE)




pH <- as.data.frame(cbind(as.factor(as.character(production$Sample_ID)),production$Phase,production$Treatment,production$X1hrpH,production$X3hrpH,production$X5hrpH,production$X25hrpH))
colnames(pH) <- c("sampleID","Phase","Treatment","1h","3h","5h","25h")
pH$sampleID <- as.factor(as.character(production$Sample_ID));pH$Treatment <- as.factor(as.character(production$Treatment));pH$Phase <- as.factor(as.character(production$Phase))
pH$ID <- NA; pH$ID <- paste0(pH$sampleID,"-",pH$Phase)
pH <- pH[,c(ncol(pH),1:(ncol(pH)-1))]
#pH$Treatment <- gsub("bottom","Bottom",pH$Treatment); pH$Treatment <- gsub("top","Top",pH$Treatment)


pH_join <- melt(pH[,-2],c("ID","Phase","Treatment"))
pH_join$variable2 <- NA
pH_join$variable2[which(pH_join$variable=="1h")] <- "1"
pH_join$variable2[which(pH_join$variable=="3h")] <- "3"
pH_join$variable2[which(pH_join$variable=="5h")] <- "5"
pH_join$variable2[which(pH_join$variable=="25h")] <- "25"
pH_join$variable <- NULL
 colnames(pH_join) [5]<- "variable"


variable <- factor(pH_join$variable, levels = c("1","3","5","25"))

pH_join$variable <- as.factor(variable)
pH_PLOT  <-  ggplot(data=pH_join,
                             aes(x=variable, y=value, fill=Treatment)) +
  geom_boxplot(size=2,alpha=0.8)+
  #   geom_point()+
  ylab("pH")+ xlab("hours")+
  ggtitle("")+   
  scale_fill_manual("legend", values = c("blue","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=40),axis.text.x  = element_text(size=40,colour="black",angle = 90, hjust = 1),
        axis.text.y  = element_text(size=40,colour="black"),
        legend.title=element_blank())# +
  
pH_PLOT

#  scale_x_continuous(limits = c(0, 25),breaks=seq(0,25,1))
                    
ggsave(paste0(graph_dir,"/","pH_behavior",".pdf"),pH_PLOT,dpi=600,width =90,height=34.395,units = "cm",scale=1.2,limitsize = FALSE)











