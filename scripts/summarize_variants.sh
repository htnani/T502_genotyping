#!/bin/bash

#### 
#workingdir=/path/to/your/T502_genotyping/directory
workingdir=/N/dc2/scratch/rtraborn/T502_genotyping
refDir=fasta
outDir=/N/dc2/scratch/rtraborn/T502_genotyping/vcf_out
bamDir=test_alignments
EC_annot=E_cloacae_annotation.gff
ECL_annot=E_coli_annotation.gff
KP_annot=K_pneumoniae_annotation.gff 
nThreads=24

gatk='java -jar /N/soft/rhel7/gatk/3.8/GenomeAnalysisTK.jar'
picard='java -jar /N/soft/rhel7/picard/2.14.0/picard.jar'

module load bedtools

cd $workingdir

echo ""
echo "=============================================================================="
echo " Summarizing mutations by gene across samples"
echo "=============================================================================="
echo ""

cd $bamDir

for VCF in `ls KP-*.vcf`; do

bedtools intersect -c -a ../annotation/$KP_annot -b $VCF > $(basename $VCF .vcf)_intersect.gff3

done

#for BM in `ls ECL-*.vcf`; do

#bedtools intersect -c -a ../annotation/$ECL_annot -b $VCF > $(basename $VCF .vcf)_intersect.gff3

#done

for VCF in `ls EC-*.vcf`; do

bedtools intersect -c -a ../annotation/$EC_annot -b $VCF > $(basename $VCF .vcf)_intersect.gff3

done

exit
