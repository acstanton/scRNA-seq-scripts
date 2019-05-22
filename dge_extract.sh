#!/bin/bash

# usage: module load java
#        sbatch dge_extract.sh </path/to/dropseq.jar> <gene_exon_tagged.bam> <outfile prefix> <number of cells to use>

# creates a gzipped matrix of transcript counts (gene x cell) from the final cleaned bamfile output by Drop-Seq postprocessing.
# input bamfile is the output of dropseq_postprocess.sh, and it should consist of cell, UMI, and exon-tagged reads. 

#SBATCH -p short
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem=15G
#SBATCH -t 0-01:00
#SBATCH -o dge_%j.out
#SBATCH -e dge_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

java -jar $1 DigitalExpression \
    I=$2 \
    O=$3.dge.txt.gz \
    SUMMARY=$3.dge.summary.txt \
    NUM_CORE_BARCODES=$4

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
