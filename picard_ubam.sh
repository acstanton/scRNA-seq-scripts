#!/bin/bash


# usage: module load java; module load picard
#        sbatch picard_ubam.sh <file1.fastq> <file2.fastq> <outfile.bam> <sample name>

# creates an unaligned bamfile from paired fastq files (i.e. from paired-end sequencing). read 1 typically labeled
# eg file_1_sequence.fastq, etc.

#SBATCH -p short
#SBATCH -t 0-01:00
#SBATCH -o picard_ubam_%j.out
#SBATCH -e picard_ubam_%j.err
#SBATCH --mem=20G
#SBATCH -c 2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

java -jar $PICARD/picard-2.8.0.jar FastqToSam \
    F1= $1 \
    F2= $2 \
    O= $3 \
    SM= $4

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
