#!/bin/bash

# Script that extracts chromosome 5 from each of the input files and writes the reads to a new BAM file.

cat $1 | while read bam
do
	
	#samtools index ../$bam
	# Get all mapped reads where both mates are mapped
	samtools view -bh -f 1 -F 12 ../$bam chr5 > $bam.paired.chr5.bam
	
	# Get reads where R1 is unmapped, R2 is mapped in a primary alignment
	samtools view -bh -f 4 -F 264 ../$bam chr5 > $bam.matemapped.chr5.bam
 
	# Get reads where R1 is mapped, and R2 is unmapped
	samtools view -bh -f 8 -F 260 ../$bam chr5 > $bam.mateunmapped.chr5.bam

	# merge bam files
	samtools merge $bam.merged.chr5.bam  $bam.paired.chr5.bam $bam.matemapped.chr5.bam $bam.mateunmapped.chr5.bam

	# sort merged file
	samtools sort -n $bam.merged.chr5.bam $bam.merged.chr5.sorted
	samtools index $bam.merged.chr5.sorted.bam
	bedtools bamtofastq -i $bam.merged.chr5.sorted.bam -fq	$bam.chr5.R1.fq -fq2	$bam.chr5.R2.fq
done


