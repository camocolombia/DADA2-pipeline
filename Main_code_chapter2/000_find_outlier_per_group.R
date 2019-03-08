
require("phyloseq")
require("plyr")
require("dplyr")
require("ggplot2")

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}



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

#Calling production data
# production <- read.csv(paste0(production_dir,"/","production_data.csv"),header=T,sep=",")
# ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))
production <- read.csv(paste0(production_dir,"/","production_data.csv"),header=T,sep=",")
quality_file <- read.csv(paste0(seq_dir,"/","track_quality.csv"),header=T,sep="|")
#ps
colnames(quality_file)[1] <-"subject"
data_to_use <- left_join(quality_file,production,by="subject")
data_to_use <- data_to_use[-c(46:48),]



data_to_use[is_outlier(data_to_use$nonchim),]
data_to_use2 <- data_to_use %>% tibble::rownames_to_column(var="outlier") %>% group_by(Breed) %>% 
  mutate(is_outlier=ifelse(is_outlier(input), input, as.numeric(NA)))
data_to_use2$outlier[which(is.na(data_to_use2$is_outlier))] <- as.numeric(NA)

for(i in 1:nrow(data_to_use2)){
  if(is.na(data_to_use2$outlier[[i]])){
    data_to_use2$outlier[[i]] <- NA
  } else {
    data_to_use2$outlier[[i]] <- data_to_use2$subject[[i]]
    
  }
};rm(i)

coverage_boxplot<- ggplot(data=data_to_use2,aes(x=Breed,y=input,fill=Breed))+
  # geom_bar(stat="identity",show.legend = T,position = "dodge")+
  geom_boxplot(position = "dodge",outlier.colour = "black",outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T,show.legend = T) +
  stat_boxplot(geom ='errorbar') +
  
  geom_text(aes(label=outlier),na.rm=TRUE,nudge_y=0.05,size=15) +
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
coverage_boxplot


ggsave(paste0(graph_dir,"/","Outlier_samples_per_group",".pdf"),coverage_boxplot,dpi=600,width =240,height=90,units = "cm",scale=0.8,limitsize = FALSE)


final_table <-data.frame(Subject=data_to_use2$subject,Phase =data_to_use2$Phase,Breed=data_to_use2$Breed,Outlier=data_to_use2$outlier)
final_table <- final_table[which(!is.na(final_table$Outlier)),]
write.csv(final_table,paste0(seq_dir,"/","Outlier.csv"),quote = F,row.names = F)
