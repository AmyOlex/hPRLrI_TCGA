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
FQ_OUT=$DIR/chr5fastq
BAM_OUT=$DIR/chr5bam

mkdir -p $FQ_OUT
mkdir -p $BAM_OUT

echo "Starting Chr5 Extraction "`date`

cat $INPUT | while read bam
do
  	echo "Starting extraction of $bam on "`date`
        prefix=$(basename $bam .bam)
        
	if [[ -e $bam.bai ]]
	then
		echo "index already exists...moving to extraction."
	else
		echo "Indexing bam $prefix..."
		samtools index $bam
	fi	

	echo "Getting all mapped reads where both mates are mapped..."
	samtools view -bh -f 1 -F 12 $bam chr5 > $BAM_OUT/$prefix.paired.chr5.bam
	
	echo "Getting reads where R1 is unmapped, R2 is mapped in a primary alignment..."
	samtools view -bh -f 4 -F 264 $bam chr5 > $BAM_OUT/$prefix.matemapped.chr5.bam
 
	echo "Getting reads where R1 is mapped, and R2 is unmapped..."
	samtools view -bh -f 8 -F 260 $bam chr5 > $BAM_OUT/$prefix.mateunmapped.chr5.bam

	echo "Merge filtered bam files for $prefix..."
	samtools merge $BAM_OUT/$prefix.merged.coordSorted.chr5.bam  $BAM_OUT/$prefix.paired.chr5.bam $BAM_OUT/$prefix.matemapped.chr5.bam $BAM_OUT/$prefix.mateunmapped.chr5.bam

	echo "Sort merged file on read names..."
	samtools sort -n $BAM_OUT/$prefix.merged.coordSorted.chr5.bam $BAM_OUT/$prefix.merged.nameSorted.chr5
	samtools index $BAM_OUT/$prefix.merged.coordSorted.chr5.bam
	bedtools bamtofastq -i $BAM_OUT/$prefix.merged.nameSorted.chr5.bam -fq $FQ_OUT/$prefix.chr5.R1.fq -fq2 $FQ_OUT/$prefix.chr5.R2.fq
	
	echo "Deleting temporary files..."
	rm $BAM_OUT/$prefix.paired.chr5.bam
	rm $BAM_OUT/$prefix.matemapped.chr5.bam
	rm $BAM_OUT/$prefix.mateunmapped.chr5.bam
        
        echo "Finished extraction of $prefix on "`date`

done

echo "Completed Chr5 Extraction on "`date`


