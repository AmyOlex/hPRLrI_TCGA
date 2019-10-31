#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run Salmon on a list of BAM files aligned to the transcriptome.

REQUIRED ARGUMENTS:
   	-f      Input file with list of BAM files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 04_Salmon will be created here.

-f INPUT File:
Input file with list of BAM files to process.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/04_Salmon.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/

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


INPUT=$FILE
OUTDIR=$DIR/04_Salmon

mkdir -p $OUTDIR

echo "Starting Salmon on "`date`

cat $INPUT | while read file
do
	echo "Starting read counting of $file on "`date` 
	prefix=$(basename $file .Aligned.toTranscriptome.out.bam)
	echo "outputting results to $OUTDIR/$prefix"
	salmon quant -t /data/refGenomes/human/GDC-GRCh38/GRCh38.d1.vd1.transcriptome.fa -l IU -a $file -o $OUTDIR/$prefix
	mv $OUTDIR/$prefix/quant.sf $OUTDIR/$prefix.quant.sf
	echo $prefix
	echo "Finished read counting of $file on "`date`
done

echo "Completed Salmon on "`date`
