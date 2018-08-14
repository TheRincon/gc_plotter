#!/bin/bash

# index for order
samtools faidx $1

# get sizes of chromosomes
cut -f 1,2 $1.fai > $1.sizes

# make windows from bins
bedtools makewindows -g $1.sizes -w $2 > $1_$2bps.bed

# get GC
bedtools nuc -fi $1 -bed $1_$2bps.bed > $1_nuc.txt

# gwak and print only the bases we need in nice way
gawk -v w=$2 'BEGIN{FS="\t"; OFS="\t"} { if (FNR>1) {print $1,$2,$3,"GC_pc"w"bps",$5} }' $1_nuc.txt > $1_$2bps.igv

# horizontal ideograms   ->  Rscript gcplot.r test.fasta_1000bps.igv 1000
if [ $3 = "d" ]
then
	Rscript gcplot.r $1_$2bps.igv $2
fi

# circular plot   ->  Rscript cirplot.r test.fasta_10000bps.igv 10000 jpg 0.5
if [ $3 = "c" ]
then
	Rscript cirplot.r $1_$2bps.igv $2 $4 $5
fi

