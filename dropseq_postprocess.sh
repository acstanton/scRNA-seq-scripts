#!/bin/bash

# usage: module load java; module load picard
#        sbatch dropseq_postprocess.sh <staraligned.sam> <reference.fasta> <unaligned_mc_tagged_polyA_filtered.bam> <path/to/dropseq.jar> <annotations.refFlat>

# for postprocessing of aligned scRNA-seq data. takes sam file from STAR output and returns bamfile with reads
# corresponding to exons tagged for dge extraction. will accept annotations in either .gtf or .refFlat format,
# but refFlat has proven more manipulatable upstream.

# unless troubleshooting is necessary, delete all .tmp.bam files after completion.

#SBATCH -p short
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-04:00
#SBATCH --mem=15G
#SBATCH -o postprocess_%j.out
#SBATCH -e postprocess_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

java -jar $PICARD/picard-2.8.0.jar SortSam \
    I=$1 \
    O=$1_sorted.tmp.bam \
    SO=queryname

java -jar $PICARD/picard-2.8.0.jar MergeBamAlignment \
    REFERENCE_SEQUENCE=$2 \
    UNMAPPED_BAM=$3 \
    ALIGNED_BAM=$1_sorted.tmp.bam \
    OUTPUT=$1_merged.tmp.bam \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    PAIRED_RUN=false

java -jar $4 TagReadWithGeneExon \
    I=$1_merged.tmp.bam \
    O=$1_star_gene_exon_tagged.bam \
    ANNOTATIONS_FILE=$5 \
    TAG=GE

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
