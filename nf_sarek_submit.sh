#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

genomeVer=WS295
genomeDir=/mnt/external.data/MeisterLab/publicData/genomes/${genomeVer}
genomeFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.genomic.fa
#gtfFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.canonical_geneset.gtf

WORK_DIR=/mnt/external.data/MeisterLab/Kalyan/EMS_sequencing_data
CONFIG_FILE=/mnt/external.data/MeisterLab/nf-core/unibe_izb.config


nextflow run nf-core/sarek -profile singularity -r 3.7.1 --input ${WORK_DIR}/csv/markduplicates_no_table.csv --step variant_calling --multiqc_title multiqc_sarek --outdir $WORK_DIR -c $CONFIG_FILE --aligner bwa-mem2 --genome WBcel235 --no_intervals --skip_tools baserecalibrator --tools freebayes,mpileup,snpeff --length_required 50 --join_germline --freebayes_filter 30 --normalize_vcfs 

