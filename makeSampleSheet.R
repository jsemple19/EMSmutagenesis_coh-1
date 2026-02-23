library(dplyr)

if(Sys.info()['sysname']=="Darwin"){
  serverPath="/Volumes/external.data/MeisterLab"
} else if (Sys.info()['sysname']=="Linux"){
  serverPath="/mnt/external.data/MeisterLab"
} else if (Sys.info()['sysname']=="Windows"){
  serverPath="Z:/MeisterLab"
}


workDir=paste0(serverPath,"/Kalyan/EMS_sequencing_data")

fq<-read.csv("./fastqList.txt",header=F)

nrow(fq)
idx<-seq(1, nrow(fq),by=2)

df<-data.frame(patient = "tm580",
               sex = NA,
               status = 0,
               sample=NA,
               lane = 0,
               fastq_1=fq[idx,],
               fastq_2=fq[idx+1,],
               group=NA)

df$sample<-sapply(strsplit(df$fastq_1,'/'), "[[", 9)
df$status<-ifelse(grepl("EMS",df$sample),1,0)
df$group<-gsub("_CROSS","",df$sample)

write.csv(df,file=paste0(workDir,"/fileList.csv"),quote=F,row.names=F)


df<-data.frame(patient = "tm580",
               sex = NA,
               status = 0,
               sample=NA,
               lane = 0,
               cram = NA,
               crai = NA,
               group=NA,
               fastq_1=fq[idx,])

df$sample<-sapply(strsplit(df$fastq_1,'/'), "[[", 9)
df$status<-ifelse(grepl("EMS",df$sample),1,0)
df$group<-gsub("_CROSS","",df$sample)
df$cram<-paste0(workDir,"/preprocessing/markduplicates/",df$sample,"/",df$sample,".md.cram")
df$crai<-paste0(workDir,"/preprocessing/markduplicates/",df$sample,"/",df$sample,".md.cram.crai")
df$fastq_1<-NULL

write.csv(df,file=paste0(workDir,"/fileList_variant_calling.csv"),quote=F,row.names=F)


