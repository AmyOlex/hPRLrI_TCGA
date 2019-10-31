#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

# Script that extracts chromosome 5 from each of the input files and writes the reads to a new BAM file.

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to extract chromosome 5 from each of the input files and writes the reads to a new BAM file. Only read pairs with at least one mate mapped to chr5 are selected. The BAM file is then converted back into a paired FastQ file for re-alignment to the transcriptome and read counting done by Salmon.

REQUIRED ARGUMENTS:
        -f	Input file with list of BAM files to process.  *Must include full path to files.*

OPTIONAL ARGUMENTS:
        -d              Directory of where results should be saved.  Sub-folder 01_chr5 will be created here.

-f INPUT File:
Input file with list of BAM files to process.  *Must include full path to files.*

-d RESULTS Directory:
Directory where results folder should be created.

EXAMPLE USAGE:
   >> ~/01_extract_chr5.sh -f /path/to/filelist.txt -d /path/to/results/dir/

EOF
}


## Parsing input arguments

FILE=
DIR=
while getopts “hf:d:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         f)
             FILE=$OPTARG
             ;;
         d)
             DIR=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

## Handle empty arguments

if [[ -z $FILE ]]; then usage; exit 1; fi
if [[ -z $DIR ]]; then usage; exit 1; fi

cd $DIR

INPUT=$FILE
OUTDIR=$DIR/01_chr5fastq

mkdir -p $OUTDIR

echo "Starting Chr5 Extraction "`date`

cat $INPUT | while read bam
do
  	echo "Starting extraction of $bam on "`date`
        prefix=$(basename $bam .bam)
        
	echo "Indexing bam $prefix..."
	samtools index $bam
	
	echo "Getting all mapped reads where both mates are mapped..."
	samtools view -bh -f 1 -F 12 $bam chr5 > $OUTDIR/$prefix.paired.chr5.bam
	
	echo "Getting reads where R1 is unmapped, R2 is mapped in a primary alignment..."
	samtools view -bh -f 4 -F 264 $bam chr5 > $OUTDIR/$prefix.matemapped.chr5.bam
 
	echo "Getting reads where R1 is mapped, and R2 is unmapped..."
	samtools view -bh -f 8 -F 260 $bam chr5 > $OUTDIR/$prefix.mateunmapped.chr5.bam

	echo "Merge filtered bam files for $prefix..."
	samtools merge $OUTDIR/$prefix.merged.chr5.bam  $OUTDIR/$prefix.paired.chr5.bam $OUTDIR/$prefix.matemapped.chr5.bam $OUTDIR/$prefix.mateunmapped.chr5.bam

	echo "Sort merged file..."
	samtools sort -n $OUTDIR/$prefix.merged.chr5.bam $OUTDIR/$prefix.merged.chr5.sorted
	samtools index $OUTDIR/$prefix.merged.chr5.sorted.bam
	bedtools bamtofastq -i $OUTDIR/$prefix.merged.chr5.sorted.bam -fq $OUTDIR/$prefix.chr5.R1.fq -fq2 $OUTDIR/$prefix.chr5.R2.fq
	
	echo "Deleting temporary files..."
	rm $OUTDIR/$prefix.paired.chr5.bam
	rm $OUTDIR/$prefix.matemapped.chr5.bam
	rm $OUTDIR/$prefix.mateunmapped.chr5.bam
        
        echo "Finished extraction of $prefix on "`date`

done

echo "Completed Chr5 Extraction on "`date`


