# EMSmutagenesis_coh-1

Scripts for analysing sequencing data from EMS mutagenesis strains and their cross progeny from coh-1 mutant (tm580) suppressor screen performed by Kalyan Ghadage.

## nf-sarek detection of variants

Used nf-core [sarek](https://nf-co.re/sarek/3.7.1) pipeline to detect variants. 

Used `makeSampleSheet.R` script to prepare input sample sheet. Designed it in tumour vs normal mode, with the orignal unmutagenised strain (PMW1074) as the normal tissue and the mutagenised cloned (original and backcrossed) as tumour samples. 

Used `nf_sarek_submit.sh` script to perform the analysis. Important to use normalization parameter to avoid too much polymorphism. Aligned with bwa-mem2, called variants with freebayes and annotated with snpEff. 

## Quality filtering of variants

A lot of variants called by freebayes simply arise from bad sequencing quality/alignment, or are common to all samples. The script `qualityFilterVCF.R` performs stringent quality filtering and plots frequency of alternative alleles along the genome. 

## EMS forest

The quality filtered variants (~200-400 per sample) were tested with the [EMS forest](https://github.com/young55775/EMSForest-A-Machine-Learning-Enhanced-EMS-Mutagenesis-Probability-Map) algorithm.

The `makeUpdateGeneRanges.R` was used to create gene ranges from WS295 gtf file.

The `EMSforest_submit.sh` script was used to run the EMS forest model on the quality filtered SNPs

The `post_ems_forest.R` script was used to filter EMS forest results to only keep mutations in coding regions that had high or moderate impact and had a significant -log10(p-value).

The `qualityFilterVCF.R` script incorporated some of these genes in the plots. (manually chosen after comparing the frequency plots and EMS forest results).





