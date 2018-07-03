##################################################
###Chrystian C. Sosa 2018 Chapter 1              #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##################################################

trim_adapters_function <-  function(code,trim_folder,trim_version,miseq_path,names_seq,seq_proc_dir_up,seq_proc_dir_p){
  
  setwd(miseq_path)
  
  paired_name1 <- sub("_L001_R1_001.fastq.gz","_R1_paired.fastq.gz",names_seq[[1]])
  paired_name2 <- sub("_L001_R2_001.fastq.gz","_R2_paired.fastq.gz",names_seq[[2]])
  unpaired_name1 <- sub("_L001_R1_001.fastq.gz","_R1_unpaired.fastq.gz",names_seq[[1]])
  unpaired_name2 <- sub("_L001_R2_001.fastq.gz","_R2_unpaired.fastq.gz",names_seq[[2]])
  
   # sink(code)
    com <- paste0('java -jar',' ',paste0(trim_folder,"/",trim_version),' ')
    PE <- paste0('PE',' ',names_seq[[1]],' ',names_seq[[2]],' ')
    R1 <- paste0(paste0(seq_proc_dir_up,"/",unpaired_name1),' ',paste0(seq_proc_dir_p,"/",paired_name1),' ')
    R2 <- paste0(paste0(seq_proc_dir_up,"/",unpaired_name2),' ',paste0(seq_proc_dir_p,"/",paired_name2),' ')

    option1 <- paste0('ILLUMINACLIP:',paste0(trim_folder,"/","adapters","/","TruSeq3-PE-2.fa"),':2:30:10',' ')
    option2 <-paste0('LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36')
    
    statement <-paste0(com,PE,R1,R2,option1,option2)

    # cat('java -jar','',paste0(trim_folder,"/",trim_version), fill = T)
    # cat('PE','',names_seq[[1]],'',names_seq[[2]],'', fill = T)
    # cat(paste0(seq_proc_dir,"/",unpaired_name1),'',paste0(seq_proc_dir,"/",paired_name1),'', fill = T)
    # cat(paste0(seq_proc_dir,"/",unpaired_name2),'',paste0(seq_proc_dir,"/",paired_name2),'', fill = T)
    # 
    # cat('ILLUMINACLIP:',paste0(trim_folder,"/","adapters","/","TruSeq3-PE-2.fa"),':2:30:10','',fill = T)
    # cat(paste0('LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36'), fill = T)
    
    write.table(statement,code,sep = "",quote = F,row.names = F,col.names = F)
    #shell(code)# system2(paste0('python ', code));# shell.exec(code)
  system(statement)
  return(statement)
     print('Done...')
  
}



#x <- trim_adapters_function (code = paste0(seq_proc_dir,"/","trim.txt"),trim_folder,trim_version,miseq_path,names_seq,seq_proc_dir_up,seq_proc_dir_p)