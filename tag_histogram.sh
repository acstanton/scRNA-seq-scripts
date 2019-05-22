#!/bin/bash

# usage: module load java
#        sbatch tag_histogram.sh </path/to/dropseq.jar> <out_gene_exon_tagged.bam> <outfile.txt.gz>

# creates a histogram of number of reads overlapping an exon in each cell

#SBATCH -p short
#SBATCH -c 1
#SBATCH -t 0-00:10
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o taghist_%j.out
#SBATCH -e taghist_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

java -jar $1 BAMTagHistogram \
    I=$2 \
    O=$3 \
    TAG=XC

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
