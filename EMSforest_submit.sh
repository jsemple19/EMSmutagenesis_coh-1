#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE EMSforest

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

WORK_DIR=/mnt/external.data/MeisterLab/Kalyan/EMS_sequencing_data/_normalized
MODEL_FILE=$WORK_DIR/EMSForest/model_1_2023_3_2.pkl
VCF_DIRS=( EMS_1_CROSS_vs_PMW1074 EMS_1_vs_PMW1074 EMS_2_CROSS_vs_PMW1074 EMS_2_vs_PMW1074 EMS_3_vs_PMW1074 )
GENE_RANGES=$WORK_DIR/EMSForest/gene_range_WS295.csv
THRESHOLD=20
OUTDIR=$WORK_DIR/EMSforest_Results
#BACKGROUND=$WORK_DIR/annotation/freebayes/PMW1074/PMW1074.freebayes.filtered.norm.sorted_snpEff.ann.vcf
#gunzip -k ${BACKGROUND}.gz


for VCF_DIR in ${VCF_DIRS[@]}
do
	echo $VCF_DIR
	mkdir -p $OUTDIR/$VCF_DIR
	VCF_FILE=${WORK_DIR}/qualityFilt/${VCF_DIR}/${VCF_DIR}.freebayes.filtered.norm.sorted_snpEff.ann.vcf
	gunzip -k ${VCF_FILE}.gz
	python $WORK_DIR/EMSForest/compare.py --model $MODEL_FILE --data $WORK_DIR/qualityFilt/$VCF_DIR --ref $GENE_RANGES --threshold $THRESHOLD --out $OUTDIR/$VCF_DIR 
	#--background $BACKGROUND
done
