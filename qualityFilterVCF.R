library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(vcfR)
library(ggplot2)
library(ggVennDiagram)
library(plotly)
library(tidyr)
library(htmlwidgets)
library(ggrepel)



options(tibble.width=Inf)

someGenes<-data.frame(chr=c("X","IV","V","III","III","IV","IV"),
                      start=c(7465193, 10561905, 3188343, 9060221,3612951,16646791,8620683),
                      end=c(7475018, 10563974, 3190115, 9071472,3627142,16655883,8626162),
                      name=c("coh-1","him-8","lag-2","lin-12","sel-8","set-26","lin-49"))

blacklist<-c("bath-35")

workDir<-"/Volumes/external.data/MeisterLab/Kalyan/EMS_sequencing_data/_normalized"
#workDir<-"/Volumes/external.data/MeisterLab/Kalyan/EMS_sequencing_data/"
dir.create(paste0(workDir,"/qualityFilt"), showWarnings = FALSE)
dir.create(paste0(workDir,"/plots"), showWarnings = FALSE)

inputDirs<-list.dirs(paste0(workDir,"/annotation/freebayes"), full.names = TRUE, recursive = FALSE)

groups<-basename(inputDirs)
samples<-gsub("_vs_PMW1074","",groups)



####### functions -------

#' Extract most severe annotation from VCF
#'
#' @param vcf vcfR object
#'
#' @returns data frame of most severe annotation per variant
#' @export
extractAnnotation<-function(vcf){
  # Extract variant metadata including INFO
  fix <- getFIX(vcf) %>% as.data.frame()
  # Extract ANN field only
  ann_raw <- extract_info_tidy(vcf, info_fields = "ANN")
  # Combine variant metadata + ANN
  df <- bind_cols(fix, ann_raw)
  # Expand ANN into multiple rows
  ann_expanded <- df %>% separate_rows(ANN, sep = ",") %>% # one row per annotation
    separate( ANN, into = c("Allele","Effect","Impact","Gene","Gene_ID", "Feature_Type","Feature_ID","Transcript_Biotype", "Rank","HGVS_c","HGVS_p","cDNA_pos","CDS_pos", "AA_pos","Distance","Errors_Warnings"),
              sep = "\\|", fill = "right" )
  ann_most_severe <- ann_expanded %>% mutate(Impact = factor(Impact, levels = c("HIGH","MODERATE","LOW","MODIFIER"))) %>% group_by(Key) %>% slice_min(Impact) %>% slice_head() %>% arrange(Key) %>% ungroup()
  return(ann_most_severe)
}


#' Extract allele counts per sample from VCF
#'
#' @param vcf vcfR object
#' @param sample name of sample of interest
#'
#' @returns table of allele counts, depth, alt allele frequency and genotype for sample of interest and control (PMW1074)
#' @export
extractAlleleCounts<-function(vcf, sample){
  # extract allele counts
  refCounts<-extract.gt(vcf, element = 'RO',as.numeric=T)
  altCounts<-extract.gt(vcf, element = 'AO',as.numeric=T)

  colnames(refCounts)<-gsub("tm580_","refCounts_",colnames(refCounts))
  colnames(altCounts)<-gsub("tm580_","altCounts_",colnames(altCounts))

  # create data frame of counts
  counts<-data.frame(cbind(refCounts,altCounts))
  colMeans(counts,na.rm=T)
  counts$id<-rownames(counts)
  counts$chr<-strsplit(counts$id,"_") %>% sapply("[[",1)
  counts$start<-as.numeric(strsplit(counts$id,"_") %>% sapply("[[",2))

  # calculate depth
  counts[[paste0("depth_",sample)]]<-rowSums(counts[,grep(paste0("Counts_",sample),colnames(counts)),drop = FALSE])
  counts[["depth_PMW1074"]]<-rowSums(counts[,grep("Counts_PMW1074",colnames(counts)),drop = FALSE])

  # calculate alt Freq
  counts[[paste0("altFreq_",sample)]]<-round(counts[[paste0("altCounts_",sample)]] / counts[[paste0("depth_",sample)]],2)
  counts[["altFreq_PMW1074"]]<-round(counts[["altCounts_PMW1074"]] / counts[["depth_PMW1074"]],2)

  # extract genotype
  genotype<-extract.gt(vcf, element = 'GT')
  colnames(genotype)<-gsub("tm580_","genotype_",colnames(genotype))
  counts<-cbind(counts,genotype)
  table(counts$genotype_PMW1074,counts[,paste0("genotype_",sample)])

  # extract most severe annotation
  ann<-extractAnnotation(vcf)
  counts$ann_impact<-ann$Impact
  counts$gene<-ann$Gene
  counts$gene_id<-ann$Gene_ID
  counts$HGVS_c<-ann$HGVS_c
  counts$HGVS_p<-ann$HGVS_p
  counts$hoverLabel<-paste("ID:", counts$id, "<br>Gene:", counts$gene,
                           "<br>HGVS_c:", counts$HGVS_c,
                            ifelse(counts$HGVS_p!= "",
                                paste0("<br>HGVS_p:",counts$HGVS_p),""))
  #extract impact
  info<-getINFO(vcf)
  counts$impact<-"LOW"
  counts$impact[grep("MODERATE",info)]<-"MODERATE"
  counts$impact[grep("HIGH",info)]<-"HIGH"
  counts$impact<-factor(counts$impact, levels=c("LOW","MODERATE","HIGH"))
  return(counts)
}


#' Quality filter a VCF
#'
#' @param vcf vcf object
#' @param sample name of sample of interest
#' @param minAltFreq Minimum alternative allele frequency for sample of interest
#' @param maxAltFreq Maximum alternative allele frequency for sample of interest
#'
#' @returns index for variants that pass quality filters
#' @export
qualityFilterVCF<-function(vcf, sample,minAltFreq=0.1,maxAltFreq=1){
  counts<-extractAlleleCounts(vcf, sample)
  #depth
  dp<-extract.gt(vcf, element = 'DP',as.numeric=TRUE)
  medians<-apply(dp,2,median,na.rm=T)
  highCounts<- dp[,1] > 3*medians[1] & dp[,2] > 3*medians[2]
  #Mean mapping quality of observed alternate alleles
  mqm<-extract.info(vcf, element = 'MQM',as.numeric=TRUE)
  #Mean mapping quality of observed reference alleles
  mqmr<-extract.info(vcf, element = 'MQMR',as.numeric=TRUE)
  #Strand balance probability for the reference allele
  srp<-extract.info(vcf, element = 'SRP',as.numeric=TRUE)
  #Strand balance probability for the alternate allele
  sap<-extract.info(vcf, element = 'SAP',as.numeric=TRUE)
  qr<-extract.info(vcf, element = 'QR',as.numeric=TRUE)
  qa<-extract.info(vcf, element = 'QA',as.numeric=TRUE)

  qual <- as.numeric(vcf@fix[, "QUAL"])
  qd_dp <- qual / dp # Quality by depth
  #Proportion of observed alternate alleles which are supported by properly paired read fragments"
  paired<-extract.info(vcf, element = 'PAIRED',as.numeric=TRUE)
  #Proportion of observed reference alleles which are supported by properly paired read fragments"
  pairedr<-extract.info(vcf, element = 'PAIREDR',as.numeric=TRUE)
  #Total number of alternate alleles in called genotypes
  ac<-extract.info(vcf, element = 'AC',as.numeric=TRUE) #
  type<-extract.info(vcf, element = 'TYPE')

  idx<-which(counts[,paste0("depth_",sample)]>10 &
             counts[,"depth_PMW1074"]>10 &
             counts[,"altFreq_PMW1074"]<0.00001 &
             mqm > 30 & mqmr > 30 &
             abs(mqm-mqmr) < 10 &
             rowSums(qd_dp>2)==2 &
             paired > 0.5 & pairedr > 0.5 &
             ac > 0 & ac < 3 &
             counts[,paste0("altFreq_",sample)]>=minAltFreq &
             counts[,paste0("altFreq_",sample)]<=maxAltFreq &
             type!="complex" &
             !highCounts)

  return(list(index=(1:nrow(counts) %in% idx),alleleCounts=counts))
}

#' Plot alternative allele frequency across genome
#'
#' @param ll list of two objects, an index for vcf file and a table of counts
#' @param sample name of sample of interest
#' @param landmarks optional data frame of genomic landmarks to plot
#' @returns plot of alternative allele ferquency in genome
#' @export
plotAltFreq<-function(ll,sample,landmarks=NULL){
    ss<-ll$alleleCounts[ll$index,]
    p<-ggplot(data=ss,aes(x=start/1e6,
                      y=.data[[paste0("altFreq_", sample)]],
                     text = hoverLabel))  +
      facet_wrap(~chr,scales="free_x",nrow=1) +
      geom_point(aes(shape=impact),alpha=0.4) +
      theme_classic() +
      labs(x="Genomic position (Mb)", y=paste0("Alternative allele frequency - ",sample)) +
      ggtitle(sample) + ylim(0,1)

    if(!is.null(landmarks)){
      p<-p+ geom_vline(data=landmarks, aes(xintercept=start/1e6), color="grey40", linetype="dashed",alpha=0.6,inherit.aes = FALSE) +
        geom_text(data=landmarks, aes(x=start/1e6, y=0.9, label=name), angle=90, vjust=-0.5, size=3,color="grey40",alpha=0.6,inherit.aes = FALSE)
    }

    print(table(ss$genotype_PMW1074,ss[,paste0("genotype_",sample)]))
    return(p)
  }

###############-
mutNames<-list()

# fitler VCFs -----

# EMS-1: 2 recessive genes, selecting egl suppression

# EMS-1
g1=groups[2]
print(g1)
annot_file1<-paste0(workDir,"/annotation/freebayes/",g1,"/",g1,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz")
sample1<-gsub("_vs_PMW1074","",g1)
vcf1<-read.vcfR(annot_file1, verbose = FALSE)

ll1<-qualityFilterVCF(vcf1, sample1,minAltFreq = 0.1, maxAltFreq = 1)
vcf1[ll1$index]
dir.create(paste0(workDir,"/qualityFilt/",g1), showWarnings = FALSE)
write.vcf(vcf1[ll1$index],file=paste0(workDir,"/qualityFilt/",g1,"/",g1,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz"))

fix<-getFIX(vcf1[ll1$index])
mutNames[[g1]]<-paste0(fix[,"CHROM"],"_",fix[,"POS"])
length(mutNames[[g1]])
sum(duplicated(mutNames[[g1]]))


#EMS-1 backcrossed, egl suppressor selected
g1a=groups[1]
print(g1a)
annot_file1a<-paste0(workDir,"/annotation/freebayes/",g1a,"/",g1a,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz")
sample1a<-gsub("_vs_PMW1074","",g1a)
vcf1a<-read.vcfR(annot_file1a, verbose = FALSE)
ll1a<-qualityFilterVCF(vcf1a, sample1a,minAltFreq = 0.1, maxAltFreq = 1)
dir.create(paste0(workDir,"/qualityFilt/",g1a), showWarnings = FALSE)
vcf1a[ll1a$index]
write.vcf(vcf1a[ll1a$index],file=paste0(workDir,"/qualityFilt/",g1a,"/",g1a,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz"))
fix<-getFIX(vcf1a[ll1a$index])
mutNames[[g1a]]<-paste0(fix[,"CHROM"],"_",fix[,"POS"])
sum(duplicated(mutNames[[g1a]]))
length(mutNames[[g1a]])



# EMS-2: 1 dominant gene, selecting egl
g2=groups[4]
print(g2)
annot_file2<-paste0(workDir,"/annotation/freebayes/",g2,"/",g2,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz")
sample2<-gsub("_vs_PMW1074","",g2)
vcf2<-read.vcfR(annot_file2, verbose = FALSE)
ll2<-qualityFilterVCF(vcf2, sample2,minAltFreq = 0.1, maxAltFreq = 1)
vcf2[ll2$index]
dir.create(paste0(workDir,"/qualityFilt/",g2), showWarnings = FALSE)
write.vcf(vcf2[ll2$index],file=paste0(workDir,"/qualityFilt/",g2,"/",g2,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz"))
fix<-getFIX(vcf2[ll2$index])
mutNames[[g2]]<-paste0(fix[,"CHROM"],"_",fix[,"POS"])
length(mutNames[[g2]])
sum(duplicated(mutNames[[g2]]))


g2a=groups[3]
print(g2a)
annot_file2a<-paste0(workDir,"/annotation/freebayes/",g2a,"/",g2a,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz")
sample2a<-gsub("_vs_PMW1074","",g2a)
vcf2a<-read.vcfR(annot_file2a, verbose = FALSE)
ll2a<-qualityFilterVCF(vcf2a, sample2a,minAltFreq = 0, maxAltFreq = 1)
vcf2a[ll2a$index]
dir.create(paste0(workDir,"/qualityFilt/",g2a), showWarnings = FALSE)
write.vcf(vcf2a[ll2a$index],file=paste0(workDir,"/qualityFilt/",g2a,"/",g2a,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz"))
fix<-getFIX(vcf2a[ll2a$index])
mutNames[[g2a]]<-paste0(fix[,"CHROM"],"_",fix[,"POS"])
length(mutNames[[g2a]])
sum(duplicated(mutNames[[g2a]]))



# EMS-3: ?

g3=groups[5]
print(g3)
annot_file3<-paste0(workDir,"/annotation/freebayes/",g3,"/",g3,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz")
sample3<-gsub("_vs_PMW1074","",g3)
vcf3<-read.vcfR(annot_file3, verbose = FALSE)
ll3<-qualityFilterVCF(vcf3, sample3,minAltFreq = 0.2, maxAltFreq = 1)
vcf3[ll3$index]
dir.create(paste0(workDir,"/qualityFilt/",g3), showWarnings = FALSE)
write.vcf(vcf3[ll3$index],file=paste0(workDir,"/qualityFilt/",g3,"/",g3,".freebayes.filtered.norm.sorted_snpEff.ann.vcf.gz"))
fix<-getFIX(vcf3[ll3$index])
mutNames[[g3]]<-paste0(fix[,"CHROM"],"_",fix[,"POS"])
length(mutNames[[g3]])
sum(duplicated(mutNames[[g3]]))



# plots -----
notUnique<-names(table(unlist(mutNames[c(g1,g2,g3)]))[table(unlist(mutNames[c(g1,g2,g3)]))>1])

### EMS-1 & CROSS plot -----
common<-intersect(mutNames[[g1]],mutNames[[g1a]])
length(common)

wierdAlleleFreq<-c("set-26", "tatn-1")
highImpact<-c("thoc-3", "C11E4.6", "Y53G8AM.4")
interesting<-c("lin-49","nhr-38")

p1<-plotAltFreq(ll1,sample1,landmarks=someGenes)
ll1$alleleCounts$group<-NA
ll1$alleleCounts[ll1$alleleCounts$id %in% common,"group"]<-paste0(sample1, "& CROSS")
ll1$alleleCounts[ll1$alleleCounts$id %in% notUnique,"group"]<-"Background"
ll1$alleleCounts[ll1$alleleCounts$gene %in% wierdAlleleFreq & ll1$alleleCounts$impact!="LOW","group"]<-"wierdAlleleFreq"
ll1$alleleCounts[ll1$alleleCounts$gene %in% highImpact & ll1$alleleCounts$impact=="HIGH","group"]<-"highImpact"
ll1$alleleCounts$group<-factor(ll1$alleleCounts$group, levels=c("Background",paste0(sample1, "& CROSS"),"wierdAlleleFreq","highImpact"))
p1<-p1 + geom_point(data=ll1$alleleCounts[!is.na(ll1$alleleCounts$group),], aes(x=start/1e6,y=.data[[paste0("altFreq_", sample1)]],color=group,shape=impact), size=0.8,alpha=0.7) +
  #geom_point(data=ll1$alleleCounts[ll1$alleleCounts$id %in% notUnique,], aes(x=start/1e6,y=.data[[paste0("altFreq_", sample1)]],color=group), size=0.2)+
  scale_color_manual(values=c("red","green","orange","purple"))
p1

notLow1a<-ll1a$alleleCounts %>% filter(id %in% common, impact != "LOW", altFreq_EMS_1_CROSS>0.6, altFreq_PMW1074<0.00001)
notLow1<-ll1$alleleCounts %>% filter(id %in% notLow1a$id, impact != "LOW", altFreq_EMS_1>0.6, altFreq_PMW1074<0.00001)

wierdAlleleFreqTbl<-ll1$alleleCounts %>% filter(gene %in% wierdAlleleFreq, impact != "LOW")
wierdAlleleFreqTbl<-left_join(wierdAlleleFreqTbl, ll1a$alleleCounts[,c("id","altFreq_EMS_1_CROSS")], by="id")

highImpactTbl<-ll1$alleleCounts %>% filter(gene %in% highImpact, impact == "HIGH")
highImpactTbl<-left_join(highImpactTbl, ll1a$alleleCounts[,c("id","altFreq_EMS_1_CROSS")], by="id")

interestingTbl<-ll1$alleleCounts %>% filter(gene %in% interesting, impact != "LOW")
interestingTbl<-left_join(interestingTbl, ll1a$alleleCounts[,c("id","altFreq_EMS_1_CROSS")], by="id")


topHits<-left_join(notLow1,notLow1a[,c("id","altFreq_EMS_1_CROSS")], by="id")

toIgnore<-c("B0491.t2")
topHits<-topHits[!(topHits$gene %in% toIgnore),]
toLabel<-rbind(topHits, wierdAlleleFreqTbl,highImpactTbl,interestingTbl)
toLabel<-toLabel[!duplicated(toLabel$id),]
toLabel<-toLabel[!(toLabel$gene %in% toIgnore) & !(toLabel$id %in% notUnique),]

p1<-p1 + geom_text_repel( data = toLabel, aes(x = start/1e6, y = .data[[paste0("altFreq_",sample1)]], label = gene), vjust = -0.5, size=2,
                          max.overlaps=Inf,  min.segment.length = 0, segment.color = "grey20", segment.size = 0.2, segment.alpha = 0.5, force = 2)

saveWidget(ggplotly(p1), file = paste0(workDir,"/plots/",g1,"_altFreq_plot.html"))


p1a<-plotAltFreq(ll1a,sample1a,landmarks=someGenes)
ll1a$group<-NA
ll1a$alleleCounts[ll1a$alleleCounts$id %in% common,"group"]<-paste0(sample1, "& CROSS")
ll1a$alleleCounts[ll1a$alleleCounts$id %in% notUnique,"group"]<-"Background"
ll1a$alleleCounts[ll1a$alleleCounts$gene %in% wierdAlleleFreq & ll1a$alleleCounts$impact!="LOW","group"]<-"wierdAlleleFreq"
ll1a$alleleCounts[ll1a$alleleCounts$gene %in% highImpact & ll1a$alleleCounts$impact=="HIGH","group"]<-"highImpact"
ll1a$alleleCounts$group<-factor(ll1a$alleleCounts$group, levels=c("Background",paste0(sample1, "& CROSS"),"wierdAlleleFreq","highImpact"))
p1a<-p1a + geom_point(data=ll1a$alleleCounts[!is.na(ll1a$alleleCounts$group),], aes(x=start/1e6,y=.data[[paste0("altFreq_", sample1a)]],color=group,shape=impact), size=0.8,alpha=0.7)+
  #geom_point(data=ll1a$alleleCounts[ll1a$alleleCounts$id %in% notUnique,], aes(x=start/1e6,y=.data[[paste0("altFreq_", sample1a)]],color=group), size=0.2) +
  scale_color_manual(values=c("red","green","orange","purple"))

p1a<-p1a + geom_text_repel( data = toLabel, aes(x = start/1e6, y = .data[[paste0("altFreq_",sample1a)]], label = gene), vjust = -0.5, size=2,
                            max.overlaps=Inf,  min.segment.length = 0, segment.color = "grey20", segment.size = 0.2, segment.alpha = 0.5, force = 2)
p1a

saveWidget(ggplotly(p1a), file = paste0(workDir,"/plots/",g1a,"_altFreq_plot.html"))

p<-ggpubr::ggarrange(p1,p1a,ncol=1,nrow=2)

ggsave(filename = paste0(workDir,"/plots/",g1,"_altFreq_plot.png"), plot = p, width = 10, height = 6)
topHits$hoverLabel<-NULL
write.csv(topHits, file = paste0(workDir,"/plots/",g1,"_topHitsByFreq.csv"), row.names = FALSE)
toLabel$hoverLabel<-NULL
write.csv(toLabel, file = paste0(workDir,"/plots/",g1,"_topHits.csv"), row.names = FALSE)



### EMS-2 & CROSS plot -----
common<-intersect(mutNames[[g2]],mutNames[[g2a]])
length(common)

highImpact<-c("Y49G5B.1","lips-6", "srw-89", "cyp-35A1", "cyp-33D1", "nlp-47", "C32D5.11")
interesting<-c("mdt-10", "hda-11", "top-1", "chd-1")

p2<-plotAltFreq(ll2,sample2,landmarks=someGenes[1:5,])
ll2$group<-NA
ll2$alleleCounts[ll2$alleleCounts$id %in% common,"group"]<-paste0(sample2, "& CROSS")
ll2$alleleCounts[ll2$alleleCounts$id %in% notUnique,"group"]<-"Background"
ll2$alleleCounts[ll2$alleleCounts$gene %in% highImpact & ll2$alleleCounts$impact=="HIGH","group"]<-"highImpact"
ll2$alleleCounts[ll2$alleleCounts$gene %in% interesting & ll2$alleleCounts$impact!="LOW","group"]<-"chromatin"
ll2$alleleCounts$group<-factor(ll2$alleleCounts$group, levels=c("Background",paste0(sample2, "& CROSS"),"highImpact","chromatin"))
p2<-p2 + geom_point(data=ll2$alleleCounts[!is.na(ll2$alleleCounts$group),], aes(x=start/1e6,y=.data[[paste0("altFreq_", sample2)]],color=group,shape=impact), size=0.8,alpha=0.7)+
#geom_point(data=ll2$alleleCounts[ll2$alleleCounts$id %in% notUnique,], aes(x=start/1e6,y=.data[[paste0("altFreq_", sample2)]],color=group), size=0.2)+
  scale_color_manual(values=c("red","green","purple","navyblue"))

notLow2a<-ll2a$alleleCounts[ll2a$index,] %>% filter(altFreq_EMS_2_CROSS>0.1) # remove any allele that is detected at high freq in cross
notLow2<-ll2$alleleCounts[ll2$index,] %>% filter(impact != "LOW", altFreq_EMS_2>0.2, altFreq_PMW1074<0.00001)

highImpactTbl<-ll2$alleleCounts %>% filter(gene %in% highImpact, impact == "HIGH")
highImpactTbl<-left_join(highImpactTbl, ll2a$alleleCounts[,c("id","altFreq_EMS_2_CROSS")], by="id")

interestingTbl<-ll2$alleleCounts %>% filter(gene %in% interesting, impact != "LOW")
interestingTbl<-left_join(interestingTbl, ll2a$alleleCounts[,c("id","altFreq_EMS_2_CROSS")], by="id")


topHits<-notLow2[!(notLow2$id %in% notLow2a$id),]
topHits<-left_join(topHits,notLow2a[,c("id","altFreq_EMS_2_CROSS")], by="id")
topHits

toLabel<-rbind(topHits,highImpactTbl,interestingTbl)
toLabel<-toLabel[!duplicated(toLabel$id),]
toLabel<-toLabel[!(toLabel$id %in% notUnique),]

p2<-p2 + geom_text_repel( data = toLabel, aes(x = start/1e6, y = .data[[paste0("altFreq_",sample2)]], label = gene), vjust = -0.5, size=2,
                          max.overlaps=Inf,  min.segment.length = 0, segment.color = "grey20", segment.size = 0.2, segment.alpha = 0.5, force = 2)



saveWidget(ggplotly(p2), file = paste0(workDir,"/plots/",g2,"_altFreq_plot.html"))


p2a<-plotAltFreq(ll2a,sample2a,landmarks=someGenes[1:5,])
ll2a$group<-NA
ll2a$alleleCounts[ll2a$alleleCounts$id %in% common,"group"]<-paste0(sample2, "& CROSS")
ll2a$alleleCounts[ll2a$alleleCounts$id %in% notUnique,"group"]<-"Background"
ll2a$alleleCounts[ll2a$alleleCounts$gene %in% highImpact & ll2a$alleleCounts$impact=="HIGH","group"]<-"highImpact"
ll2a$alleleCounts[ll2a$alleleCounts$gene %in% interesting & ll2a$alleleCounts$impact!="LOW","group"]<-"chromatin"
ll2a$alleleCounts$group<-factor(ll2a$alleleCounts$group, levels=c("Background",paste0(sample2, "& CROSS"),"highImpact","chromatin"))

p2a<-p2a + geom_point(data=ll2a$alleleCounts[!is.na(ll2a$alleleCounts$group),], aes(x=start/1e6,y=.data[[paste0("altFreq_", sample2a)]],color=group,shape=impact), size=0.8, alpha=0.7)+
#geom_point(data=ll2a$alleleCounts[ll2a$alleleCounts$id %in% notUnique,], aes(x=start/1e6,y=.data[[paste0("altFreq_", sample2a)]],color=group), size=0.2) +
  scale_color_manual(values=c("red","green","purple","navyblue"))

p2a<-p2a + geom_text_repel( data = toLabel[!is.na(toLabel$altFreq_EMS_2_CROSS),], aes(x = start/1e6, y = .data[[paste0("altFreq_",sample2a)]], label = gene), vjust = -0.5, size=2,
                          max.overlaps=Inf,  min.segment.length = 0, segment.color = "grey20", segment.size = 0.2, segment.alpha = 0.5, force = 2)

saveWidget(ggplotly(p2a), file = paste0(workDir,"/plots/",g2a,"_altFreq_plot.html"))

p<-ggpubr::ggarrange(p2,p2a,ncol=1,nrow=2)
ggsave(filename = paste0(workDir,"/plots/",g2,"_altFreq_plot.png"), plot = p, width = 10, height = 6)
topHits$hoverLabel<-NULL
write.csv(topHits, file = paste0(workDir,"/plots/",g2,"_topHitsByFreq.csv"), row.names = FALSE)
toLabel$hoverLabel<-NULL
write.csv(toLabel, file = paste0(workDir,"/plots/",g2,"_topHits.csv"), row.names = FALSE)


### EMS-3 Low penetrance. no backcross ------
topHits<-ll3$alleleCounts %>% filter(impact != "LOW", altFreq_EMS_3>0.8,
                           !(id %in% notUnique), altFreq_PMW1074<0.00001)
dim(topHits) #56
emsforest<-read.csv(paste0(workDir,"/EMSforest_Results/EMSforest_suspect_genes_0.05.csv"))
emsforest<-emsforest[emsforest$group=="EMS_3_vs_PMW1074",]
dim(emsforest) #53

emsforestvar<-ll3$alleleCounts %>% filter(gene_id %in% emsforest$gene, impact!="LOW", altFreq_PMW1074<0.00001)
dim(emsforestvar) # 53

p3<-plotAltFreq(ll3,sample3,landmarks=someGenes[1:4,])
ll3$group<-NA
ll3$alleleCounts[ll3$alleleCounts$id %in% notUnique,"group"]<-"Background"
ll3$alleleCounts[ll3$alleleCounts$id %in% topHits$id,"group"]<-"highAltFreq"
ll3$alleleCounts[ll3$alleleCounts$id %in% emsforestvar$id,"group"]<-"EMSforest"
ll3$alleleCounts[ll3$alleleCounts$id %in% topHits$id & ll3$alleleCounts$id %in% emsforestvar$id,"group"] <-"Freq&forest"

ll3$alleleCounts$group<-factor(ll3$alleleCounts$group, levels=c("Background","highAltFreq","EMSforest","Freq&forest"))

table(ll3$alleleCounts$group)

p3<-p3 + geom_point(data=ll3$alleleCounts[!is.na(ll3$alleleCounts$group),], aes(x=start/1e6,y=.data[[paste0("altFreq_", sample3)]], color=group, shape=impact), size=0.8, alpha=0.7)+
  scale_color_manual(values=c("red","pink3","lightblue","purple"))
p3


toLabel<-ll3$alleleCounts[!is.na(ll3$alleleCounts$group) & ll3$alleleCounts$group!="Background",]
toLabel<-toLabel[!duplicated(toLabel$id),]

p3<-p3 + geom_text_repel( data = toLabel, aes(x = start/1e6, y = .data[[paste0("altFreq_",sample3)]], label = gene, color=group), vjust = -0.5, size=2,
                          max.overlaps=Inf,  min.segment.length = 0, segment.color = "grey20", segment.size = 0.2, segment.alpha = 0.5, force = 2)

saveWidget(ggplotly(p3), file = paste0(workDir,"/plots/",g3,"_altFreq_plot.html"))

p3
ggsave(filename = paste0(workDir,"/plots/",g3,"_altFreq_plot.png"), plot = p3, width = 10, height = 3)

topHits$hoverLabel<-NULL
write.csv(topHits, file = paste0(workDir,"/plots/",g3,"_topHitsByFreq.csv"), row.names = FALSE)

toLabel$hoverLabel<-NULL #60
write.csv(toLabel, file = paste0(workDir,"/plots/",g3,"_topHits.csv"), row.names = FALSE)
