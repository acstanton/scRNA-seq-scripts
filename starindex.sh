#!/bin/bash

 
# usage: module load star
#        sbatch starindex.sh </path/to/genomedir> </path/to/refgenome.fasta> </path/to/refgenome.gtf> 

# creates a STAR-indexed genome directory called </path/to/genomedir>. note that this user input is NOT the location of
# the reference fasta file, but the name of a directory to be created by STAR.

#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 0-03:00
#SBATCH -p short
#SBATCH --mem=40G
#SBATCH -o STAR_%j.out
#SBATCH -e STAR_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $1 --genomeFastaFiles $2 --sjdbGTFfile $3 --sjdbOverhang 100

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
