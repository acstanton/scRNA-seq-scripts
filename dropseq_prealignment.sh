#!/bin/bash

# usage: module load java; module load picard
#        sbatch dropseq_prealignment.sh <path/to/dropseq.jar> <bamfile.bam> <adapter sequence> 

# takes an unaligned bamfile of single cell sequencing reads and performs the pre-alignment processing
# necessary for STAR. This includes (in order) tagging with cell barcode, tagging with molecular (UMI)
# barcode, filtering out low quality reads, trimming remaining adapter sequences, trimming remaining
# polyA tails, and converting to a fastq file.

# unless troubleshooting is necessary, delete all files ending in .tmp.bam after completion.

#SBATCH -p short
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH -t 0-05:00
#SBATCH -o prealign_%j.out
#SBATCH -e prealign_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=acstanton@g.harvard.edu

java -jar $1 TagBamWithReadSequenceExtended \
    INPUT=$2 \
    OUTPUT=$2_tagged_Cell.tmp.bam \
    SUMMARY=$2_tagged_Cellular.bam_summary.txt \
    BASE_RANGE=1-12 \
    BASE_QUALITY=10 \
    BARCODED_READ=1 \
    DISCARD_READ=False \
    TAG_NAME=XC \
    NUM_BASES_BELOW_QUALITY=1

java -jar $1 TagBamWithReadSequenceExtended \
    INPUT=$2_tagged_Cell.tmp.bam \
    OUTPUT=$2_tagged_CellMolecular.tmp.bam \
    SUMMARY=$2_tagged_CellMolecular.bam_summary.txt \
    BASE_RANGE=13-20 \
    BASE_QUALITY=10 \
    BARCODED_READ=1 \
    DISCARD_READ=True \
    TAG_NAME=XM \
    NUM_BASES_BELOW_QUALITY=1

java -jar $1 FilterBAM \
    TAG_REJECT=XQ \
    INPUT=$2_tagged_CellMolecular.tmp.bam \
    OUTPUT=$2_tagged_filtered.tmp.bam

java -jar $1 TrimStartingSequence \
    INPUT=$2_tagged_filtered.tmp.bam \
    OUTPUT=$2_tagged_trimmed_smart.tmp.bam \
    OUTPUT_SUMMARY=$2_adapter_trimming_report.txt \
    SEQUENCE=$3 \
    MISMATCHES=0 \
    NUM_BASES=5

java -jar $1 PolyATrimmer \
    INPUT=$2_tagged_trimmed_smart.tmp.bam \
    OUTPUT=$2_mc_tagged_polyA_filtered.bam \
    OUTPUT_SUMMARY=$2_polyA_trimming_report.txt \
    MISMATCHES=0 \
    NUM_BASES=7

java -jar $PICARD/picard-2.8.0.jar SamToFastq \
    INPUT=$2_mc_tagged_polyA_filtered.bam \
    FASTQ=$2_mc_tagged_polyA_filtered.fastq

sleep 5
sacct --format=JobID,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID
