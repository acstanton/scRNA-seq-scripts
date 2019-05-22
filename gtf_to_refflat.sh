#!/bin/bash

# usage: module load java
#        sbatch gtf_to_refflat.sh <path/to/dropseq.jar> <annotations.gtf> <sequence.dict> <outfile.refFlat>

# converts a gtf annotation file to refFlat format

#SBATCH -p short
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:10
#SBATCH -o gtf_to_refflat_%j.out
#SBATCH -e gtf_to_refflat_%j.err
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

java -jar $1 ConvertToRefFlat \
    ANNOTATIONS_FILE=$2 \
    SEQUENCE_DICTIONARY=$3 \
    O=$4

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID

