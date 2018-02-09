#!/bin/bash

#PBS -N WGS_genotyping_variant_calling_T502
#PBS -k o
#PBS -l nodes=1:ppn=16,vmem=40gb
#PBS -l walltime=6:00:00
#PBS -m abe

workingdir=/N/dc2/scratch/rtraborn/T502_genotyping
outDir=/N/dc2/scratch/rtraborn/T502_genotyping/vcf_out
bamDir=alignments
EC_ref=E_cloacae_genomic.fasta
ECL_ref=E_coli_genome.fasta
KP_ref=K_pneumoniae_genome.fasta
nThreads=16

module load samtools
module load bcftools

cd $workingdir

        echo ""
	echo "=============================================================================="
	echo "Step 0 Generating fasta index for references ..."
	echo "=============================================================================="
	echo ""

	samtools faidx fasta/$EC_ref
	samtools faidx fasta/$ECL_ref
	samtools faidx fasta/$KP_ref

        echo ""
	echo "=============================================================================="
	echo "Step 1 Generating BCF file from .bam ..."
	echo "=============================================================================="
	echo ""

cd ${bamDir}
	for file1 in `ls KP-*.bam`; do
	    echo "samtools index $file1"
	    samtools index $file1
	    echo "samtools mpileup -g -f $reference $file1"
	    samtools mpileup -g -f ${workingdir}/fasta/$KP_ref $file1 \
		      > $(basename $file1 .bam).raw.bcf

	    echo "Step 2 Generating VCF file"
	    echo "=============================================================================="
	    echo ""
	    echo "bcftools view -bvcg $(basename $file1 .bam).bcf > $(basename $file1 .bam).bcf"
	    bcftools view -bvcg $(basename $file1 .bam).raw.bcf > $outDir/$(basename $file1 .bam).bcf

	    echo "VCF file created"
	    echo ""
	    done

	for file1 in `ls ECL-*.bam`; do
	    echo "samtools index $file1"
	    samtools index $file1
	    echo "samtools mpileup -g -f $reference $file1"
	    samtools mpileup -g -f ${workingdir}/fasta/$ECL_ref $file1 \
		      > $(basename $file1 .bam).raw.bcf

	    echo "Step 2 Generating VCF file"
	    echo "=============================================================================="
	    echo ""
	    echo "bcftools view -bvcg $(basename $file1 .bam).bcf > $(basename $file1 .bam).bcf"
	    bcftools view -bvcg $(basename $file1 .bam).raw.bcf > $outDir/$(basename $file1 .bam).bcf

	    echo "VCF file created"
	    echo ""
	    done

	for file1 in `ls EC-*.bam`; do
	    echo "samtools index $file1"
	    samtools index $file1
	    echo "samtools mpileup -g -f $reference $file1"
	    samtools mpileup -g -f ${workingdir}/fasta/$EC_ref $file1 \
		      > $(basename $file1 .bam).raw.bcf

	    echo "Step 2 Generating VCF file"
	    echo "=============================================================================="
	    echo ""
	    echo "bcftools view -bvcg $(basename $file1 .bam).bcf > $(basename $file1 .bam).bcf"
	    bcftools view -bvcg $(basename $file1 .bam).raw.bcf > $outDir/$(basename $file1 .bam).bcf

	    echo "VCF file created"
	    echo ""
	    done


echo ""
echo "=============================================================================="
echo " Finished generating files"
echo "=============================================================================="
echo ""

exit
