#!/bin/bash

# usage: module load java; module load picard
#        sbatch picard_samfastq.sh <file.bam> <output.fastq>

# converts a bamfile to a fastq file

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 0-03:00
#SBATCH --mem=20G
#SBATCH -o samfastq_%j.out
#SBATCH -e samfastq_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

java -jar $PICARD/picard-2.8.0.jar SamToFastq \
    I=$1 \
    FASTQ=$2

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
