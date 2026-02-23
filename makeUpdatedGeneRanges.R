library(rtracklayer)

gtfFile<-"/Volumes/external.data/MeisterLab/publicData/genomes/WS295/c_elegans.PRJNA13758.WS295.canonical_geneset.gtf"

gtf <- import(gtfFile)
gtf<-gtf[gtf$type=="gene"]
mcols(gtf)<-mcols(gtf)[,c("source","type","gene_id","gene_biotype","gene_name")]
seqlevels(gtf)<-c("I","II","III","IV","V","X","MtDNA")
gtf<-sort(gtf)
export(gtf,"/Volumes/external.data/MeisterLab/Kalyan/EMS_sequencing_data/TvNmode/EMSForest/c_elegans.PRJNA13758.WS295.canonical_geneset.genes.gtf")

gtfNoStrand<-gtf
strand(gtfNoStrand)<-"*"
gtfNoStrand<-sort(gtfNoStrand)
gene_ranges<-data.frame(gene_id=gtfNoStrand$gene_id,
                        chr=seqnames(gtfNoStrand),
                        start=start(gtfNoStrand),
                        end=end(gtfNoStrand))

write.table(gene_ranges,"/Volumes/external.data/MeisterLab/Kalyan/EMS_sequencing_data/TvNmode/EMSForest/gene_range_WS295.csv",row.names=F,quote=F,col.names=F,sep=",")


