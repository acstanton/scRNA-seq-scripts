-----------------------
This file provides information on the alignment and processing of scRNA-seq data. Many of the tools used in this pipeline are derived from the Drop-Seq computational pipeline developed by the McCarroll lab. 

Significant coding background is NOT necessary for the implementation of this pipeline. However, a basic understanding of bash scripting and knowledge of how to navigate in the command line environment will be helpful.

It should be possible to condense this pipeline by using pipes or reducing it to a single script, which would avoid creating a large number of temporary files. However, this implementation does not do so in order to allow the user to troubleshoot at several stages during the process and to improve resource management.

This pipeline has been optimized for use on O2, a linux-based high-performance computing cluster at Harvard Medical School. Use on other systems or in other environments may require therefore require some tweaking of settings. In particular, O2 uses the SLURM job scheduler. If this pipeline is employed locally or on a cluster using LSF or another scheduler, the appropriate changes will need to be made.

Allie Stanton, 2018
acstanton@g.harvard.edu
-----------------------

Contents:

1. Included scripts and implementation overview
2. Necessary packages
3. Necessary files and metadata
4. Processing and alignment
5. DGE extraction guide

-----------------------

1. Included scripts and implementation overview

This pipeline contains several scripts used to create the metadata files used in the actual data processing. They are:
     - starindex.sh || create a STAR-indexed genome
     - createreference.sh || create a .dict file from reference fasta file
     - gtf_to_refflat.sh || convert a .gtf annotations file to a .refFlat file
     - gb2gtf.py || convert a genbank file to .gtf format

The actual pipeline, starting with paired fastq files, is implemented as follows:
    1. picard_ubam.sh || create unaligned bamfile from paired fastq files
    2. dropseq_prealignment.sh || tag reads with cell and molecular barcodes and filter reads
    3. staralign.sh || align scRNA-seq data to reference genome with STAR
    4. dropseq_postprocess.sh || takes tagged and filtered bamfile and STAR-aligned sam file and outputs tagged aligned bamfile

Finally, two scripts are included to help the user interpret the results and extract a gene expression matrix.
     - tag_histogram.sh || create a histogram of readcounts that can be used to select cells for extraction
     - dge_extract.sh || create a digital gene expression (dge) matrix from the tagged aligned bamfile

-----------------------

2. Necessary packages

java
	used by picard and dropseq tools

python 2.7; biopython, numpy, scipy
	used by gb2gtf.py

STAR
	for indexing and alignment. version 2.5.4b

picard
	version 2.8.0, used for file management

dropseq tools
	a homebrew package that can be downloaded from the broad institute github page. version 1.13

-----------------------

3. Necessary files and metadata

a. Input fastq files

The input scRNA-seq data will be in fastq format (if it is instead in Illumina basecalls, convert it first to fastq format with Illumina's bcl2fastq). It should be paired-end data, where each pair consists of a 50bp read (read 2) corresponding to actual cellular mRNA and a 20 bp read (read 1) that corresponds to the cell and molecular barcodes. Bases 1-12 of read 1 are the cell barcode. They are unique to each bead (and therefore each cell) and all of the oligomers on a single bead have the same cell barcode. Bases 13-20 of read 1 are the molecular barcode, identifying one mRNA molecule.

b. Fasta reference genome

A reference genome in .fasta format is required for STAR indexing and alignment. For the human genome, this can be downloaded as a pair with the .gtf annotations file from ENSEMBL or other sources.

c. Reference annotations

Reference annotations in .gtf format is required for STAR. For the human genome, this can be downloaded as a pair with the .fasta file from ENSEMBL or other sources. For other genomes, if a .gtf file is not readily available, one can be created with the included gb2gtf.py script from a genbank file. If the .gtf file doesn't contain labeled exons (eg if it is from a virus), it will be useful to convert the .gtf file to a .refFlat file for dropseq_postprocess.sh, as refFlat is more amenable to manual editing so that appropriate viral "exons" can be added. This can be done with create_refflat.sh, which itself requires a .dict file, which can be created with createreference.sh.

d. STAR-indexed reference genome

A STAR-indexed reference genome is required for STAR alignment. It can be created with starindex.sh.

-----------------------

4. Processing and alignment

a. Convert paired fastq files to unaligned bamfile with picard_ubam.sh.

This script will take the paired fastq input files and convert them to an unaligned bam format file. 

b. Perform pre-alignment processing steps on the unaligned bamfile with dropseq_prealignment.sh.

This script performs several preprocessing steps that prepare the data for alignment with STAR. First, the program extracts the cell barcode--the first 12 bases of read 1--and tags the read with it. Then, it extracts the remaining 8 bases of read 1--the molecular (UMI) barcode--and tags the read with it. This will later allow us to associate each read with a particular transcript from a particular cell. After both barcodes are extracted, read 1 is dropped, which means that the output bam consists of tagged unpaired reads. Additionally, during these steps, if the quality of a base pair in a read is too low, the read is dropped. This parameter can be changed in the script. The script also trims leftover adapter starting sequences (provided by the user) and polyA tails. Summary text files are output from these steps. Finally, the bam file is converted to fastq format. The user should delete .tmp.bam files produced by this script provided that additional troubleshooting is not necessary. Both the final output bamfile and the fastq file are required for later steps.

c. Align the tagged and filtered fastq file to the reference genome with staralign.sh.

This script implements STAR to align the tagged and filtered fastq file to the reference genome. The default settings work reasonably well in my experience, but they can be tailored as necessary. Consult the STAR manual for more information.

d. Perform postprocessing of the STAR-aligned sam file to tag reads that overlap an exon with dropseq_postprocess.sh.

This script performs postprocessing on the aligned RNAseq data. First, it sorts the sam file output by staralign.sh and outputs it as a bamfile in order to guarantee that the alignment is sorted in queryname order. Then, it merges the aligned bamfile with the unaligned but tagged and filtered bamfile output by dropseq_prealignment.sh. This recovers the tags that were lost during alignment. Finally, the script takes a .gtf or .refFlat annotation file and uses it to tag the reads in the merged bamfile that overlap an exon. Either a .gtf or .refFlat format is accepted; however, I had difficulty making the .gtf for a ZIKV genome compatible with the program as the ZIKV genome does not have annotated exons (the viral genome is unsplices, so the concept of an exon does not apply in this context). As a workaround, I manually created a refFlat file for the ZIKV genome, labeling the CDS as an exon. The refFlat format was used simply because I found it more easy to manipulate than gtf format.

-----------------------

5. DGE extraction guide

The dropseq_postprocess.sh script outputs a tagged and cleaned bamfile, from which we want to extract a matrix of gene expression that shows the number of transcripts of each gene in each cell. We do not want to extract every cell into this matrix, as many cells will be low quality. So first we want to figure out what cells to extract. I recommend doing this using the tag_histogram.sh script, which takes the tagged and cleaned bamfile as input and outputs a text file containing the number of reads per cell barcode. This file can be exported into R, and inspection of a cumulative distribution plot can allow the user to select the "knee" of the distribution. Implement this in R as follows:

    a=read.table("exp_readcounts.txt.gz", header=F, stringsAsFactors=F)
    x=cumsum(a$V1)
    x=x/max(x)
    plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=c(1,5000))

Then run dge_extract.sh, using the number of cells at the knee of the distribution as input for how many cells to extract.
