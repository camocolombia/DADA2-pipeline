##################################################
###Chrystian C. Sosa 2018 Chapter 1              #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##Trimming sequences using Trimmomatic and bbduk #
##################################################

trim_adapters_function <-  function(program,code,program_dir,miseq_path,names_seq,seq_proc_dir_up,seq_proc_dir_p,primers,trimq,tpe,default,length_seqParameters){
  
  paired_name1 <- sub("_L001_R1_001.fastq.gz","_R1_paired.fastq.gz",names_seq[[1]])
  paired_name2 <- sub("_L001_R2_001.fastq.gz","_R2_paired.fastq.gz",names_seq[[2]])
  unpaired_name1 <- sub("_L001_R1_001.fastq.gz","_R1_unpaired.fastq.gz",names_seq[[1]])
  unpaired_name2 <- sub("_L001_R2_001.fastq.gz","_R2_unpaired.fastq.gz",names_seq[[2]])
  
  
  if(program=="trimmomatic"){
    cat("   ","\n");cat("USING TRIMOMMATIC","\n");cat("   ","\n")
    setwd(miseq_path)
    trim_folder <- paste0(program_dir,"/","Trimmomatic-0.38")
    trim_version <- "trimmomatic-0.38.jar"
 
  
   # sink(code)
    com <- paste0('java -jar',' ',paste0(trim_folder,"/",trim_version),' ')
    PE <- paste0('PE',' ',names_seq[[1]],' ',names_seq[[2]],' ')
    R1 <- paste0(paste0(seq_proc_dir_p,"/",paired_name1),' ',paste0(seq_proc_dir_up,"/",unpaired_name1),' ')
    R2 <- paste0(paste0(seq_proc_dir_p,"/",paired_name2),' ',paste0(seq_proc_dir_up,"/",unpaired_name2),' ')

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
    

    #shell(code)# system2(paste0('python ', code));# shell.exec(code)
    
  } else {
    
    cat("   ","\n");cat("USING BBDUK","\n");cat("   ","\n")
    setwd(miseq_path)
    trim_folder <- paste0(program_dir,"/","bbmap")
    trim_version <- "bbduk.sh"
    
    com <- paste0(paste0(trim_folder,"/",trim_version),' ','-Xmx1g',' ')
    in1 <- paste0('in1=',paste0(miseq_path,"/",names_seq[[1]]),' ')
    in2 <- paste0('in2=',paste0(miseq_path,"/",names_seq[[2]]),' ')
    out1 <- paste0('out1=',paste0(paste0(seq_proc_dir_p,"/",paired_name1),' '))
    out2 <- paste0('out2=',paste0(paste0(seq_proc_dir_p,"/",paired_name2),' '))
    literal <- paste0('literal=',"'",primers,"'",' ')
    ref <- paste0('ref=',trim_folder,"/","resources","/","adapters.fa",' ')
    
    
if(trimq==F){  trimq <- paste0('') } else { trimq <- paste0('qtrim=rl trimq=10 ')}
    
    if(is.null(length_seqParameters)){
      length_seqParameters <- '' } else { length_seqParameters <- paste0(length_seqParameters,' ')}
    
if(default==F){
options <- paste0('ordered=t mink=2 ktrim=l rcomp=f k=10 tbo')
} else {
  options <- paste0('ordered=t rcomp=f ktrim=l tbo')
}

if(tpe==T){tpe <- paste0(' tpe') }  else { tpe <- paste0('') }
    
  statement <- paste0(com,in1,in2,out1,out2,literal,ref,trimq,length_seqParameters,options,tpe)
  
}

  write.table(statement,code,sep = "",quote = F,row.names = F,col.names = F)
  system(statement)
  return(statement)
  cat("   ","\n");print('Done...');cat("   ","\n")
}



#x <- trim_adapters_function (code = paste0(seq_proc_dir,"/","trim.txt"),trim_folder,trim_version,miseq_path,names_seq,seq_proc_dir_up,seq_proc_dir_p)