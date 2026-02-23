library(dplyr)
library(ggplot2)
library(zoo)

options(tibble.width=Inf)
workDir<-"/Volumes/external.data/MeisterLab/Kalyan/EMS_sequencing_data"
dir.create(paste0(workDir,"/plots"), showWarnings = FALSE, recursive = TRUE)
dataDir<-paste0(workDir,"/MQsnpFilt")

csvList<-list.files(dataDir, pattern="high_quality_snp_counts.csv", full.names = T, recursive = F)

kwin=11
MQthresh=40

someGenes<-data.frame(CHROM=c("X","IV","V","III"),
                      start=c(7465193, 10561905, 3188343, 9060221),
                      end=c(7475018, 10563974, 3190115, 9071472),
                      name=c("coh-1","him-8","lag-2","lin-12"))

df_mean<-NULL
i=1
for(i in 1:length(csvList)){
  df<-read.csv(csvList[i], header=T)
  df<-df[df$MQM>=MQthresh & df$MQM>=MQthresh,]

  sampleName<-gsub(colnames(df)[7:8], pattern="tm580_(.*)_...Count", replacement="\\1")[1]

  df[,paste0(sampleName,"_coverage")]<-as.numeric(df[,7])+as.numeric(df[,8])


  df[,paste0(sampleName,"_altfreq")]<-as.numeric(df[,8])/df[,paste0(sampleName,"_coverage")]



  idx<-is.na(df[,paste0(sampleName,"_coverage")]) | df[,paste0(sampleName,"_coverage")]<10

  print(paste0("Rows with NA coverage: ", sum(idx)))
  df[idx,]

  df<-df[!idx,]
  df<-df[df$CHROM!="MtDNA",]

  tmp<-df %>% group_by(CHROM) %>%
    mutate(meanPOS=rollmean(POS,k=kwin,fill=NA,align="center"),
           meanFreqAlt=rollmean(get(paste0(sampleName,"_altfreq")),k=kwin,fill=NA,align="center"))
  tmp$sample<-sampleName
  if(is.null(df_mean)){
    df_mean<-tmp[,c("CHROM","POS","meanPOS","meanFreqAlt","sample")]
  }
  else {
    df_mean<-rbind(df_mean, tmp[,c("CHROM","POS", "meanPOS","meanFreqAlt","sample")])
  }
}

table(df_mean$sample)


df_mean$experiment<-gsub("_CROSS","",df_mean$sample)
df_mean$experiment[df_mean$experiment=="PMW1074"]<-"EMS_3"
df_mean$bg<-"bg"
df_mean$bg[grep("EMS_3",df_mean$sample)]<-"sup"
df_mean$bg[grep("_CROSS",df_mean$sample)]<-"sup"
df_mean$bg<-factor(df_mean$bg,levels=c("bg","sup"))

experiments<-unique(df_mean$experiment)

toHighlight<-list()
exp="EMS_1"
# add short list genes
shortlist<-read.table(paste0(workDir,"/TvNmode/EMSforest_Results/",exp,"_manually_selected_genes_snps_high_moderate_impact_short.tsv"), sep="\t",header=T)
shortestlist<-NULL
shortestlist<-c(shortestlist,shortlist %>% arrange(desc(.data[[paste0("minuslog10pval.",exp)]])) %>% slice(1:2) %>% pull(gene_name))
shortestlist<-c(shortestlist,shortlist %>% dplyr::filter(Annotation_Impact=="HIGH") %>% pull(gene_name))
shortestlist<-c(shortestlist,shortlist %>% arrange(PAM30) %>% slice(1:2) %>% pull(gene_name))
shortestlist<-c(shortestlist,shortlist %>% arrange(desc(.data[[paste0(exp,"_CROSS_altFreq")]])) %>% slice(1:2) %>% pull(gene_name))
colnames(shortlist)[colnames(shortlist)=="seqnames"]<-"CHROM"
toHighlight[[exp]]<-list(shortlist=shortlist,shortestlist=shortestlist)

exp="EMS_2"
# add short list genes
shortlist<-read.table(paste0(workDir,"/TvNmode/EMSforest_Results/",exp,"_manually_selected_genes_snps_high_moderate_impact_short.tsv"), sep="\t",header=T)
shortestlist<-NULL
shortestlist<-c(shortestlist,shortlist %>% arrange(desc(.data[[paste0("minuslog10pval.",exp)]])) %>% slice(1:2) %>% pull(gene_name))
shortestlist<-c(shortestlist,shortlist %>% dplyr::filter(Annotation_Impact=="HIGH") %>% pull(gene_name))
shortestlist<-c(shortestlist,shortlist %>% dplyr::filter(PAM30 <= -7)  %>% pull(gene_name))
shortestlist<-c(shortestlist,shortlist %>% dplyr::filter(is.na(.data[[paste0(exp,"_CROSS_altFreq")]]),
                .data[[paste0(exp,"_altFreq")]]>0.5) %>% pull(gene_name))
shortestlist<-c(shortestlist,shortlist %>% dplyr::filter(.data[[paste0(exp,"_CROSS_altFreq")]]<0.25) %>% pull(gene_name))
colnames(shortlist)[colnames(shortlist)=="seqnames"]<-"CHROM"
toHighlight[[exp]]<-list(shortlist=shortlist,shortestlist=shortestlist)




pdf(paste0(workDir,"/plots/MQsnpfilt_altfreq_k",kwin,"_MQ",MQthresh,"_shortlist.pdf"), width=12, height=8,onefile=T)
exp=experiments[2]
for(exp in experiments){
  df_exp<-df_mean[df_mean$experiment==exp,]
  bglabel<-unique(df_exp$sample[df_exp$bg=="bg"])
  suplabel<-unique(df_exp$sample[df_exp$bg=="sup"])
  p<-ggplot(df_exp, aes(x=meanPOS/1e6,y=meanFreqAlt,color=bg))+
    geom_line(alpha=0.8)+
    facet_wrap(~CHROM)+
    ggtitle(paste0("Alternative allele frequency: ",exp)) +
    theme_bw()+
    xlab("Genomic Position (Mb)")+
    ylab(paste0("Alt allele freq ",exp))+
    geom_hline(yintercept=0.5, linetype="dashed", color="red",alpha=0.3)+
    ylim(c(0,1))+
    scale_color_manual(values=c("grey60","orange"),labels=c(bglabel,suplabel))+
    geom_segment(aes(x=POS/1e6, xend=POS/1e6, y=0, yend=0.05), color="black",alpha=0.1)+
    geom_vline(data=someGenes, aes(xintercept=start/1e6), linetype="dotted", color="blue",alpha=0.3) +
    geom_text(data=someGenes, aes(x=(start+2000)/1e6, y=0.9, label=name), angle=90, vjust=-0.5, color="blue", size=2,alpha=0.3)
  shortlist<-toHighlight[[exp]]$shortlist
  shortestlist<-toHighlight[[exp]]$shortestlist
  p<-p +
    geom_vline(data=shortlist[shortlist$gene_name %in% shortestlist,], aes(xintercept=start/1e6), linetype="dotted", color="darkred",alpha=0.8) +
    geom_text(data=shortlist[shortlist$gene_name %in% shortestlist,], aes(x=(start+2000)/1e6, y=0.9, label=gene_name), angle=90, vjust=-0.5, color="darkred", size=2,alpha=0.8)+
    geom_vline(data=shortlist[!(shortlist$gene_name %in% shortestlist),], aes(xintercept=start/1e6), linetype="dotted", color="brown1",alpha=0.4) +
    geom_text(data=shortlist[!(shortlist$gene_name %in% shortestlist),], aes(x=(start+2000)/1e6, y=0.9, label=gene_name), angle=90, vjust=-0.5, color="brown1", size=2,alpha=0.4)

    print(p)
  #ggsave(paste0(workDir,"/plots/log10_altfreq_ratio_",sampleName,"_vs_",ctrlName,".png"), width=12, height=8)
}
dev.off()

