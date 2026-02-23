library(dplyr)
library(rtracklayer)
library(readxl)
workDir<-"/Volumes/external.data/MeisterLab/Kalyan/EMS_sequencing_data/_normalized"

gtf<-import(paste0(workDir,"/EMSForest/c_elegans.PRJNA13758.WS295.canonical_geneset.genes.gtf"))
gtf$score<-NULL
gtf$phase<-NULL

annoDir<-paste0(workDir,"/reports/snpeff/freebayes")

dirlist<-list.dirs(paste0(workDir,"/EMSforest_Results"), recursive = FALSE)

groups<-gsub(paste0(workDir,"/EMSforest_Results/"),"",dirlist)

df<-NULL
for(group in groups){
  print(group)
  volcano_table<-read.csv(paste0(workDir,"/EMSforest_Results/",group,"/volcano_table.csv"))
  volcano_table$X<-NULL
  volcano_table$group<-group
  mutanno<-read.csv(paste0(annoDir,"/",group,"/",group,".freebayes.filtered.norm.sorted_snpEff.genes.txt"), sep="\t", skip=1, header=T)
  mutanno$X.GeneName<-NULL
  mutanno$TranscriptId<-NULL
  mutanno$BioType<-NULL

  tmp<-left_join(left_join(volcano_table,data.frame(gtf),by=join_by("gene"=="gene_id")),mutanno,by=join_by("gene"=="GeneId"))
  if(is.null(df)){
    df<-tmp
  } else {
    if(sum(!(colnames(df) %in% colnames(tmp)))>0){
      missing<-colnames(df)[!(colnames(df) %in% colnames(tmp))]
      tmp[,missing]<-0
    }
    else if(sum(!(colnames(tmp) %in% colnames(df)))>0) {
      missing<-colnames(tmp)[!(colnames(tmp) %in% colnames(df))]
      df[,missing]<-0
    }
    df<-rbind(df, tmp)
  }
}

#rowSums(df[,grepl("variants_impact_",colnames(df))])
dim(df)
dim(df[df$p_value>2 & df$p_value<40,])

df<-df %>% group_by(gene)%>%
  mutate(num_groups=n())

df$num_strains<-0
ss<-df %>% filter(group %in% c("EMS_1_vs_PMW1074", "EMS_2_vs_PMW1074", "EMS_3_vs_PMW1074")) %>% group_by(gene) %>%
  mutate(num_strains=n())
df$num_strains[match(ss$gene,df$gene)]<-ss$num_strains


to_csv<-df%>%
  dplyr::filter(p_value>1.3, p_value<40,observed<3,
                gene_biotype=="protein_coding",num_groups<4,
                (variants_impact_HIGH>0 | variants_impact_MODERATE>0)) %>%
  group_by(group) %>% arrange(group,desc(variants_impact_HIGH),desc(variants_impact_MODERATE),p_value,desc(foldchange))

write.table(to_csv, paste0(workDir,"/EMSforest_Results/EMSforest_suspect_genes_0.05.csv"), sep=",", row.names=F, quote=F)



