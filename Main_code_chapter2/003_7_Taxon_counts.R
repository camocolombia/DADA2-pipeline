

##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##COUNTS PER TAXONOMIC LEVEL                     #
##################################################



require("phyloseq")
require("matrixStats")


mainDir <- "E:/DADA2"
chapter <- "chapter_2"
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
#phylo_object <- readRDS(paste0(csv_dir,"/","Phyloseq_object",".RDS"))
#FDR_list <- readRDS(paste0(out_dir,"/","csv/","FDR_list.RDS"))


##################################################
#Function




##Object examples
# genfac <- factor(tax_table(ps_object[[i]])[, "Kingdom"])
# gentab <- t(otu_table(ps_object[[i]]))




##ORIGINAL
# gentab = apply(otu_table(GlobalPatterns), MARGIN = 2, function(x) {
#   tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
# })




gentab_function <- function(otu_table_m,genfac){
  gentab = apply(otu_table_m, MARGIN = 2, function(x) {
    tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
  })  
  return(gentab)
  
}
  
###
#calling colnames

findcolnumber <- function(df, thecolumnname){
  which(colnames(df) == thecolumnname)
}




##################################################
#Calling ps objects

ps_object <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list",".RDS"))
ps_object2 <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list_trams",".RDS"))
#ps_object3 <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list_melt",".RDS"))




#counts_per_taxonomy_rel <- lapply(1:length(ps_object),function(i){

counts_per_taxonomy_rel <- list()

for(i in 1:length(ps_object)){
  i <- i
cat(i," | ","Processing: ",sub("ps_","",names(ps_object)[i])," ||| ","\n")
sample_data(ps_object[[i]])$sampleID <- trimws(sample_data(ps_object[[i]])$sampleID)
val_dat <- sample_data(ps_object[[i]])
animalsID <- unique(phyloseq::sample_data(ps_object[[i]])$sampleID)
###factor 
genfac <-  factor(tax_table(ps_object[[i]])[,i])
val_dat_HE <- unique(trimws(val_dat[which(val_dat$Treatment=="High feed efficiency"),]$sampleID))
val_dat_LE <- unique(trimws(val_dat[which(val_dat$Treatment=="Low feed efficiency"),]$sampleID))

val_dat_HE_Sam <- unique(trimws(val_dat[which(val_dat$Treatment=="High feed efficiency"),]$Subject))
val_dat_LE_Sam <- unique(trimws(val_dat[which(val_dat$Treatment=="Low feed efficiency"),]$Subject))

animalsID_list <- list()
for(l in 1:length(animalsID)){
  animalsID_list[[l]] <- rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], sampleID==animalsID[[l]]))),
                  genfac = genfac))
 # animalsID_list[[l]] <-  animalsID_list[[l]][which(animalsID_list[[l]]>0)] <- 1
  
};rm(l)


animalsID_list <-do.call(rbind,animalsID_list)
animalsID_list <- as.data.frame(t(animalsID_list))

for(l in 1:ncol(animalsID_list)){
  for(m in 1:nrow(animalsID_list)){
    
    if(animalsID_list[m,l]>0){
      animalsID_list[m,l] <- 1
    } else {
      animalsID_list[m,l] <- 0
    }
  };rm(m)
};rm(l)

colnames(animalsID_list) <- animalsID

###he animals



x <- gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]]))),
                genfac = genfac)

x_pres <- x
for(y in 1:nrow(x_pres)){
  for(z in 1:ncol(x_pres)){
    if(x_pres[y,z]>0){
      x_pres[y,z] <- 1
    } else {
      x_pres[y,z] <- 0
    }
  };rm(z)  
};rm(y)

cat("formatting data to fit","\n")
data_f <- data.frame(level = sub("ps_","",names(ps_object)[i]),
                     taxon_names = row.names(x),
                     animal_per_sample=NA,
                     animal_per_HE=NA,
                     animal_per_LE=NA,
                     count_per_sample=NA,
                     count_per_HE=NA,
                     count_per_LE=NA,
                     count_total = rowSums(x),
                     count_HE =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="High feed efficiency"))),
                                                         genfac = genfac)),
                     sd_HE =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="High feed efficiency"))),
                                                         genfac = genfac)),
                     count_LE =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="Low feed efficiency"))),
                                                         genfac = genfac)),
                     sd_LE =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="Low feed efficiency"))),
                                                        genfac = genfac)),
                     count_LI =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Phase=="liquid"))),
                                                         genfac = genfac)),
                     sd_LI =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Phase=="liquid"))),
                                                        genfac = genfac)),
                     count_SO =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Phase=="solid"))),
                                                         genfac = genfac)),
                     sd_SO =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Phase=="solid"))),
                                                         genfac = genfac)),
                     count_HE_LI =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="High feed efficiency" & Phase=="liquid"))),
                                                            genfac = genfac)),
                     sd_HE_LI =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="High feed efficiency" & Phase=="liquid"))),
                                                            genfac = genfac)),
                     count_HE_SO =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="High feed efficiency" & Phase=="solid"))),
                                                         genfac = genfac)),
                     sd_HE_SO =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="High feed efficiency" & Phase=="solid"))),
                                                            genfac = genfac)),
                     count_LE_LI =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="Low feed efficiency" & Phase=="liquid"))),
                                                            genfac = genfac)),
                     sd_LE_LI =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="Low feed efficiency" & Phase=="liquid"))),
                                                            genfac = genfac)),
                     count_LE_SO =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]], Treatment=="Low feed efficiency" & Phase=="solid"))),
                                                            genfac = genfac)),
                     
                     status_FE = NA,
                     status_Phase = NA,
                     status_HE_PHASE = NA,
                     status_LE_PHASE = NA,
                     status_LIQUID= NA,
                     status_SOLID = NA
                     
                     )


data_f$status_FE = data_f$count_HE > data_f$count_LE
data_f$status_Phase = data_f$count_LI > data_f$count_SO
data_f$status_HE_PHASE = data_f$count_HE_LI > data_f$count_HE_SO
data_f$status_LE_PHASE = data_f$count_LE_LI > data_f$count_LE_SO
data_f$status_LIQUID= data_f$count_HE_LI > data_f$count_LE_LI
data_f$status_SOLID = data_f$count_HE_SO > data_f$count_LE_SO


# data_f$status_HE_PHASE = data_f$count_HE_LI > data_f$count_HE_SO
# data_f$status_LE_PHASE = data_f$count_LE_LI > data_f$count_LE_SO

cat("counting samples","\n")

for(j in 1:nrow(data_f)){
  #cat(j,"\n")
  data_f$status_FE[[j]] = if(data_f$status_FE[[j]]==T){data_f$status_FE[[j]]="HE"}else if(data_f$count_HE[[j]] == data_f$count_LE[[j]]){
    data_f$status_FE[[j]]="EQUAL"} else {data_f$status_FE[[j]]="LE"}
  
  data_f$status_Phase[[j]] = if(data_f$status_Phase[[j]]==T){data_f$status_Phase[[j]]="LI"} else if( data_f$count_LI[[j]] == data_f$count_SO[[j]]){
  data_f$status_Phase[[j]]="EQUAL"} else {data_f$status_Phase[[j]]="SO"}
  
  data_f$status_HE_PHASE[[j]] = if(data_f$status_HE_PHASE[[j]]==T){data_f$status_HE_PHASE[[j]]="LI"} else if(data_f$count_HE_LI[[j]] == data_f$count_HE_SO[[j]]){
  data_f$status_HE_PHASE[[j]]="EQUAL"} else {data_f$status_HE_PHASE[[j]]="SO"}
  
  data_f$status_LE_PHASE[[j]] = if(data_f$status_LE_PHASE[[j]]==T){data_f$status_LE_PHASE[[j]]="LI"} else if(data_f$count_LE_LI[[j]] == data_f$count_LE_SO[[j]]){
  data_f$status_LE_PHASE[[j]]="EQUAL"} else {data_f$status_LE_PHASE[[j]]="SO"}
  
  data_f$status_LIQUID[[j]] = if(data_f$status_LIQUID[[j]]==T){data_f$status_LIQUID[[j]]="HE"} else  if(data_f$count_HE_LI[[j]] == data_f$count_LE_LI[[j]]){
  data_f$status_LIQUID[[j]]="EQUAL"} else {data_f$status_LIQUID[[j]]="LE"}
  
  data_f$status_SOLID[[j]] = if(data_f$status_SOLID[[j]]==T){data_f$status_SOLID[[j]]="HE"} else   if(data_f$count_HE_SO[[j]] == data_f$count_LE_SO[[j]]){
    data_f$status_SOLID[[j]]="EQUAL"} else {data_f$status_SOLID[[j]]="LE"}
};rm(j)

                 #    )  data_f$status_SOLID[[j]]="EQUAL"} else {data_f$status_SOLID[[j]]="LE"}
####Count per animal

animals_count <- as.data.frame(matrix(nrow = nrow(data_f),
                                      ncol = length(rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object[[i]]))),
                                                                            genfac = genfac)))
                                        
                                      ))#length(animalsID)))

row.names(animals_count) <- row.names(data_f)

animals_HE_count <-  animals_count; animals_LE_count <-  animals_count

# animalsID <- animalsID
cat("counting taxon levels"," | ", length(animalsID)," animals","\n")





data_f$animal_per_sample <- rowSums(animalsID_list)
data_f$animal_per_HE <- rowSums(animalsID_list[,val_dat_HE])
data_f$animal_per_LE <- rowSums(animalsID_list[,val_dat_LE])
data_f$count_per_sample <- rowSums(x_pres)
data_f$count_per_HE <- rowSums(x_pres[,val_dat_HE_Sam])
data_f$count_per_LE <- rowSums(x_pres[,val_dat_LE_Sam])



rm(animals_count,
   animals_HE_count,  
   animals_LE_count)

cat(" ||| ","Processed: ",sub("ps_","",names(ps_object)[i]),"\n")

counts_per_taxonomy_rel[[i]] <- data_f
#return(data_f)

rm(animalsID)
};rm(i)#)


counts_per_taxonomy_rel <- do.call(rbind,counts_per_taxonomy_rel)
row.names(counts_per_taxonomy_rel) <- 1:nrow(counts_per_taxonomy_rel)


write.csv(counts_per_taxonomy_rel,paste0(out_dir,"/","csv/","taxa_count_number.csv"),quote = F,row.names = F)

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
##RELATIVE ABUNDANCE


counts_per_taxonomy_rel <- list()

for(i in 1:length(ps_object2)){
  i <- i
  cat(i," | ","Processing: ",sub("ps_","",names(ps_object2)[i])," ||| ","\n")
  sample_data(ps_object2[[i]])$sampleID <- trimws(sample_data(ps_object2[[i]])$sampleID)
  val_dat <- sample_data(ps_object2[[i]])
  animalsID <- unique(phyloseq::sample_data(ps_object2[[i]])$sampleID)
  ###factor 
  genfac <-  factor(tax_table(ps_object2[[i]])[,i])
  val_dat_HE <- unique(trimws(val_dat[which(val_dat$Treatment=="High feed efficiency"),]$sampleID))
  val_dat_LE <- unique(trimws(val_dat[which(val_dat$Treatment=="Low feed efficiency"),]$sampleID))
  
  val_dat_HE_Sam <- unique(trimws(val_dat[which(val_dat$Treatment=="High feed efficiency"),]$Subject))
  val_dat_LE_Sam <- unique(trimws(val_dat[which(val_dat$Treatment=="Low feed efficiency"),]$Subject))
  
  animalsID_list <- list()
  for(l in 1:length(animalsID)){
    animalsID_list[[l]] <- rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], sampleID==animalsID[[l]]))),
                                                   genfac = genfac))
    # animalsID_list[[l]] <-  animalsID_list[[l]][which(animalsID_list[[l]]>0)] <- 1
    
  };rm(l)
  
  
  animalsID_list <-do.call(rbind,animalsID_list)
  animalsID_list <- as.data.frame(t(animalsID_list))
  
  for(l in 1:ncol(animalsID_list)){
    for(m in 1:nrow(animalsID_list)){
      
      if(animalsID_list[m,l]>0){
        animalsID_list[m,l] <- 1
      } else {
        animalsID_list[m,l] <- 0
      }
    };rm(m)
  };rm(l)
  
  colnames(animalsID_list) <- animalsID
  
  ###he animals
  
  
  
  x <- gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]]))),
                       genfac = genfac)
  
  x_pres <- x
  for(y in 1:nrow(x_pres)){
    for(z in 1:ncol(x_pres)){
      if(x_pres[y,z]>0){
        x_pres[y,z] <- 1
      } else {
        x_pres[y,z] <- 0
      }
    };rm(z)  
  };rm(y)
  
  cat("formatting data to fit","\n")
  data_f <- data.frame(level = sub("ps_","",names(ps_object2)[i]),
                       taxon_names = row.names(x),
                       animal_per_sample=NA,
                       animal_per_HE=NA,
                       animal_per_LE=NA,
                       count_per_sample=NA,
                       count_per_HE=NA,
                       count_per_LE=NA,
                       count_total = rowSums(x),
                       count_HE =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="High feed efficiency"))),
                                                           genfac = genfac)),
                       sd_HE =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="High feed efficiency"))),
                                                                     genfac = genfac)),
                       count_LE =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="Low feed efficiency"))),
                                                           genfac = genfac)),
                       sd_LE =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="Low feed efficiency"))),
                                                                     genfac = genfac)),
                       count_LI =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Phase=="liquid"))),
                                                           genfac = genfac)),
                       sd_LI =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Phase=="liquid"))),
                                                                     genfac = genfac)),
                       count_SO =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Phase=="solid"))),
                                                           genfac = genfac)),
                       sd_SO =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Phase=="solid"))),
                                                                     genfac = genfac)),
                       count_HE_LI =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="High feed efficiency" & Phase=="liquid"))),
                                                              genfac = genfac)),
                       sd_HE_LI =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="High feed efficiency" & Phase=="liquid"))),
                                                                        genfac = genfac)),
                       count_HE_SO =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="High feed efficiency" & Phase=="solid"))),
                                                              genfac = genfac)),
                       sd_HE_SO =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="High feed efficiency" & Phase=="solid"))),
                                                                        genfac = genfac)),
                       count_LE_LI =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="Low feed efficiency" & Phase=="liquid"))),
                                                              genfac = genfac)),
                       sd_LE_LI =  matrixStats ::rowSds(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="Low feed efficiency" & Phase=="liquid"))),
                                                                        genfac = genfac)),
                       count_LE_SO =  rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]], Treatment=="Low feed efficiency" & Phase=="solid"))),
                                                              genfac = genfac)),
                       
                       status_FE = NA,
                       status_Phase = NA,
                       status_HE_PHASE = NA,
                       status_LE_PHASE = NA,
                       status_LIQUID= NA,
                       status_SOLID = NA
                       
  )
  
  
  data_f$status_FE = data_f$count_HE > data_f$count_LE
  data_f$status_Phase = data_f$count_LI > data_f$count_SO
  data_f$status_HE_PHASE = data_f$count_HE_LI > data_f$count_HE_SO
  data_f$status_LE_PHASE = data_f$count_LE_LI > data_f$count_LE_SO
  data_f$status_LIQUID= data_f$count_HE_LI > data_f$count_LE_LI
  data_f$status_SOLID = data_f$count_HE_SO > data_f$count_LE_SO
  
  
  # data_f$status_HE_PHASE = data_f$count_HE_LI > data_f$count_HE_SO
  # data_f$status_LE_PHASE = data_f$count_LE_LI > data_f$count_LE_SO
  
  cat("counting samples","\n")
  
  for(j in 1:nrow(data_f)){
    #cat(j,"\n")
    data_f$status_FE[[j]] = if(data_f$status_FE[[j]]==T){data_f$status_FE[[j]]="HE"}else if(data_f$count_HE[[j]] == data_f$count_LE[[j]]){
      data_f$status_FE[[j]]="EQUAL"} else {data_f$status_FE[[j]]="LE"}
    
    data_f$status_Phase[[j]] = if(data_f$status_Phase[[j]]==T){data_f$status_Phase[[j]]="LI"} else if( data_f$count_LI[[j]] == data_f$count_SO[[j]]){
      data_f$status_Phase[[j]]="EQUAL"} else {data_f$status_Phase[[j]]="SO"}
    
    data_f$status_HE_PHASE[[j]] = if(data_f$status_HE_PHASE[[j]]==T){data_f$status_HE_PHASE[[j]]="LI"} else if(data_f$count_HE_LI[[j]] == data_f$count_HE_SO[[j]]){
      data_f$status_HE_PHASE[[j]]="EQUAL"} else {data_f$status_HE_PHASE[[j]]="SO"}
    
    data_f$status_LE_PHASE[[j]] = if(data_f$status_LE_PHASE[[j]]==T){data_f$status_LE_PHASE[[j]]="LI"} else if(data_f$count_LE_LI[[j]] == data_f$count_LE_SO[[j]]){
      data_f$status_LE_PHASE[[j]]="EQUAL"} else {data_f$status_LE_PHASE[[j]]="SO"}
    
    data_f$status_LIQUID[[j]] = if(data_f$status_LIQUID[[j]]==T){data_f$status_LIQUID[[j]]="HE"} else  if(data_f$count_HE_LI[[j]] == data_f$count_LE_LI[[j]]){
      data_f$status_LIQUID[[j]]="EQUAL"} else {data_f$status_LIQUID[[j]]="LE"}
    
    data_f$status_SOLID[[j]] = if(data_f$status_SOLID[[j]]==T){data_f$status_SOLID[[j]]="HE"} else   if(data_f$count_HE_SO[[j]] == data_f$count_LE_SO[[j]]){
      data_f$status_SOLID[[j]]="EQUAL"} else {data_f$status_SOLID[[j]]="LE"}
  };rm(j)
  
  #    )  data_f$status_SOLID[[j]]="EQUAL"} else {data_f$status_SOLID[[j]]="LE"}
  ####Count per animal
  
  animals_count <- as.data.frame(matrix(nrow = nrow(data_f),
                                        ncol = length(rowSums(gentab_function(otu_table_m = t(otu_table(subset_samples(ps_object2[[i]]))),
                                                                              genfac = genfac)))
                                        
  ))#length(animalsID)))
  
  row.names(animals_count) <- row.names(data_f)
  
  animals_HE_count <-  animals_count; animals_LE_count <-  animals_count
  
  # animalsID <- animalsID
  cat("counting taxon levels"," | ", length(animalsID)," animals","\n")
  
  
  
  
  
  data_f$animal_per_sample <- rowSums(animalsID_list)
  data_f$animal_per_HE <- rowSums(animalsID_list[,val_dat_HE])
  data_f$animal_per_LE <- rowSums(animalsID_list[,val_dat_LE])
  data_f$count_per_sample <- rowSums(x_pres)
  data_f$count_per_HE <- rowSums(x_pres[,val_dat_HE_Sam])
  data_f$count_per_LE <- rowSums(x_pres[,val_dat_LE_Sam])
  
  
  
  rm(animals_count,
     animals_HE_count,  
     animals_LE_count)
  
  cat(" ||| ","Processed: ",sub("ps_","",names(ps_object2)[i]),"\n")
  
  counts_per_taxonomy_rel[[i]] <- data_f
  #return(data_f)
  
  rm(animalsID)
};rm(i)#)


counts_per_taxonomy_rel <- do.call(rbind,counts_per_taxonomy_rel)
row.names(counts_per_taxonomy_rel) <- 1:nrow(counts_per_taxonomy_rel)


write.csv(counts_per_taxonomy_rel,paste0(out_dir,"/","csv/","taxa_count_rel_number.csv"),quote = F,row.names = F)

