#!/bin/bash

# usage: module load java; module load picard
#        sbatch createreference.sh <reference.fasta> <outfile.dict>

# creates a sequence dictionary from a reference fasta

#SBATCH -p short
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:05
#SBATCH -o createref_%j.out
#SBATCH -e createref_%j.err
#SBATCH --mem=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary \
    R=$1 \
    O=$2

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
