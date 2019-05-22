#!/bin/bash

# usage: module load star
#        sbatch staralign.sh <path/to/indexedgenomeDIR> <file.fastq> <outfileprefix>

# aligns a fastq file to an indexed genome using STAR. outputs a sam file along with some quality control
# summary files.

#SBATCH -p short
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --mem=40G
#SBATCH -t 0-02:00
#SBATCH -o staralign_%j.out
#SBATCH -e staralign_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

STAR --runThreadN 4 --genomeDir $1 --readFilesIn $2 --outFileNamePrefix $3

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
