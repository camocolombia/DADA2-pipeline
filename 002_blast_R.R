##################################################
###Chrystian C. Sosa 2018 Chapter 1              #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##BLAST USING AVS FASTA FILES                    #
##################################################

#Calling libraries
require(DECIPHER);require(Biostrings);require(devtools);require(taxize);require(rentrez)
#install_github("mhahsler/rBLAST")
require(rBLAST)

#Defining folders
mainDir <- "E:/DADA2"
chapter <- "chapter_1"

#Defining enviroment parameters
blast_n_path <- "C:/Program Files/NCBI/blast-2.7.1+/bin"
Sys.setenv(PATH = paste(Sys.getenv("PATH"), blast_n_path , sep= .Platform$path.sep))
Sys.which("blastn")
blast_folder <- "E:/DB_DADA2/NCBI/16SMicrobial"
key_NCBI <- "ee630947ebb03dea7f2e757a53ea2ae29208"

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
script_dir <-  paste0(mainDir,"/","scripts")

#Reading AVS fasta fies
seq <- Biostrings::readDNAStringSet(paste0(seq_dir,"/","ASVs.fa"),format="fasta")
bl <- blast(db=paste0(blast_folder,"/","16SMicrobial"))

#Running BLAST
blast_match <- lapply(1:length(seq),function(i){
cat(round((i/length(seq)*100),2)," %","\n")
  cl <- predict(bl, seq[i,])
cl <- cl[which(cl$Perc.Ident>=99),]
cl <- cl[which.min(cl$E),]

if(nrow(cl)==0){
cl <- cl[1,]
cl$MATCH <-NA;cl$MATCH <- rep("<99",nrow(cl))

  } else {
  
    cl$MATCH <-NA;
    cl$MATCH <- rep(">99",nrow(cl))
  }

return(cl)

})

#Joining BLAST results in a unique file
blast_match <- do.call(rbind,blast_match)
blast_match <- blast_match[which(blast_match$MATCH==">99"),]

#Joining taxonomic information
class_match1 <- lapply(1:nrow(blast_match),function(i){
  cat(round((i/nrow(blast_match)*100),2)," % | row: ",i,"\n")
  
  ncbi_tax_id <- rentrez::entrez_summary("Nucleotide",id=blast_match$SubjectID[[i]])$taxid
  clasif <- taxize::classification(ncbi_tax_id,db="ncbi",key=key_NCBI)
  clasif2 <- as.data.frame(t(clasif[[1]][,1]))
  colnames(clasif2) <- as.character(t(clasif[[1]][,2]))
  
  clasif2 <- clasif2[,-which(names(clasif2) %in% "no rank")]

  return(clasif2)
})


blast_matrix <- as.data.frame(matrix(ncol=8,nrow=nrow(blast_match)))
colnames(blast_matrix) <- c("superkingdom","phylum","class","order","family","genus","species","status")
for(i in 1:nrow(blast_matrix)){
  cat(i,"\n")
  blast_matrix[i,1] <- if("superkingdom" %in% names(class_match1[[i]])==T) {as.character(class_match1[[i]][,"superkingdom"]) } else {NA}
  blast_matrix[i,2] <-  if("phylum" %in% names(class_match1[[i]])==T) {as.character(class_match1[[i]][,"phylum"]) } else {NA}
  blast_matrix[i,3] <-  if("class" %in% names(class_match1[[i]])==T) {as.character(class_match1[[i]][,"class"]) } else {NA}
  blast_matrix[i,4] <-  if("order" %in% names(class_match1[[i]])==T) {as.character(class_match1[[i]][,"order"]) } else {NA}
  blast_matrix[i,5] <-  if("family" %in% names(class_match1[[i]])==T) {as.character(class_match1[[i]][,"family"]) } else {NA}
  blast_matrix[i,6] <-  if("genus" %in% names(class_match1[[i]])==T) {as.character(class_match1[[i]][,"genus"]) } else {NA}
  blast_matrix[i,7] <-  if("species" %in% names(class_match1[[i]])==T) {as.character(class_match1[[i]][,"species"]) } else {NA}
  blast_matrix[i,8] <-  if(sum(is.na(blast_matrix[i,1:7]))>0){"CHECK"} else {"OK"}
};rm(i)

blast_match <- cbind(blast_match,blast_matrix)

#Saving file 
write.csv(blast_match,paste0(seq_dir,"/","16S_BLAST.csv"),row.names = F,quote = F)
