#!/bin/bash

#PBS -N WGS_gatk_variant_calling_T502
#PBS -k o
#PBS -l nodes=1:ppn=24,vmem=40gb
#PBS -l walltime=18:00:00
#PBS -m abe

#### 
#workingdir=/path/to/your/T502_genotyping/directory
workingdir=/N/dc2/scratch/rtraborn/T502_genotyping
refDir=fasta
outDir=/N/dc2/scratch/rtraborn/T502_genotyping/vcf_out
bamDir=alignments
EC_ref=E_cloacae_genomic.fasta
ECL_ref=E_coli_genome.fasta
KP_ref=K_pneumoniae_genome.fasta
nThreads=24

gatk='java -jar /N/soft/rhel7/gatk/3.8/GenomeAnalysisTK.jar'
picard='java -jar /N/soft/rhel7/picard/2.14.0/picard.jar'

module load samtools
module load picard
module load gatk
module load r

cd $workingdir


        echo ""
	echo "=============================================================================="
	echo "Step 0 Sorting each bam file  ..."
	echo "=============================================================================="
	echo ""

cd ${bamDir}
        for BM in *.bam; do
	    echo "$samtools sort -@ 8 -o $(basename $BM .bam).sorted.bam $BM"
	          samtools sort -@ 8 -o $(basename $BM .bam).sorted.bam $BM
	    done


        echo ""
	echo "=============================================================================="
	echo "Step 0 Adding a read group to each bam file  ..."
	echo "=============================================================================="
	echo ""

for BM in *.sorted.bam; do

$picard AddOrReplaceReadGroups \
      I=$BM \
      O=$(basename $BM .sorted.bam).rh.bam \
      RGID=4 \
      RGLB=$(basename $BM .sorted.bam) \
      RGPL=Illumina \
      RGPU=unit1 \
      RGSM=20

done
        echo ""
	echo "=============================================================================="
	echo "Step 1 Performing picard markDuplicates on each bam file  ..."
	echo "=============================================================================="
	echo ""

cd ..

echo "Processing KP files"

cd ${bamDir}
	for file1 in `ls KP-*.rh.bam`; do
	    echo "$picard MarkDuplicates $file1"
	    $picard MarkDuplicates \
		I=$file1 \
		O=$(basename $file1 .rh.bam).markedDup.bam \
		M=$(basename $file1 .rh.bam)_metrics.txt
	    done

echo "Processing EC files"

	for file1 in `ls EC-*.rh.bam`; do
	    echo "$picard MarkDuplicate $file1"
	    $picard MarkDuplicates \
		I=$file1 \
		O=$(basename $file1 .rh.bam).markedDup.bam \
		M=$(basename $file1 .rh.bam)_metrics.txt
	    done

echo "Processing ECL files"

	for file1 in `ls ECL-*.rh.bam`; do
	    echo "$picard MarkDuplicate $file1"
	    $picard MarkDuplicates \
		I=$file1 \
		O=$(basename $file1 .rh.bam).markedDup.bam \
		M=$(basename $file1 .rh.bam)_metrics.txt
	done

        echo ""
	echo "=============================================================================="
	echo "Step 2 Generating VCF files from BAMs using GATK ..."
	echo "=============================================================================="
	echo ""

###################################################################
#   Prepping reference genome fasta file for GATK
echo "Prepping reference genome fasta file for GATK....."
cd ../${refDir}

# Create sequence dictionary using Picard Tools.
# the following command produces a SAM-style header file describing the contents of our fasta file.
$picard CreateSequenceDictionary \
    REFERENCE=$KP_ref \
    OUTPUT=$(basename $KP_ref .fasta).dict

$picard CreateSequenceDictionary \
    REFERENCE=$EC_ref \
    OUTPUT=$(basename $EC_ref .fasta).dict

$picard CreateSequenceDictionary \
    REFERENCE=$ECL_ref \
    OUTPUT=$(basename $ECL_ref .fasta).dict

echo "created sequence dictionary for the reference genomes."

echo "indexing the reference genomes...."

# Create the fasta index file.
# The index file describes byte offset in the fasta file for each contig. It is a text file with one record
# per line for each of the fasta contigs. Each record is of the type -
# contig, size, location, basePerLine, bytesPerLine
samtools faidx $KP_ref
samtools faidx $ECL_ref
samtools faidx $EC_ref

echo "Reference genomes are now ready for GATK."

cd ../${bamDir}

###############################################################
## Summary Statistics

for BM in `ls KP-*.markedDup.bam`; do
        $picard MeanQualityByCycle \
	    INPUT=$BM \
	    CHART_OUTPUT=$(basename $BM .markedDup.bam)_mean_quality_by_cycle.pdf \
	    OUTPUT=$(basename $BM .markedDup.bam)_read_quality_by_cycle.txt \
	    REFERENCE_SEQUENCE=../${refDir}/$KP_ref

        $picard QualityScoreDistribution \
	    INPUT=$BM \
	    CHART_OUTPUT=$(basename $BM .markedDup.bam)_mean_quality_overall.pdf \
	    OUTPUT=$(basename $BM .markedDup.bam)_read_quality_overall.txt \
	    REFERENCE_SEQUENCE=../${refDir}/$KP_ref

	$picard CollectWgsMetrics \
	    INPUT=$BM \
	    OUTPUT=$(basename $BM .markedDup.bam)_stats_picard.txt \
	    REFERENCE_SEQUENCE=../${refDir}/$KP_ref \
	    MINIMUM_MAPPING_QUALITY=20 \
	    MINIMUM_BASE_QUALITY=20
done

for BM in `ls EC-*.markedDup.bam`; do
       $picard MeanQualityByCycle \
 	    INPUT=$BM \
	    CHART_OUTPUT=$(basename $BM .markedDup.bam)_mean_quality_by_cycle.pdf \
	    OUTPUT=$(basename $BM .markedDup.bam)_read_quality_by_cycle.txt \
	    REFERENCE_SEQUENCE=../${refDir}/$EC_ref

       $picard QualityScoreDistribution \
	    INPUT=$BM \
	    CHART_OUTPUT=$(basename $BM .markedDup.bam)_mean_quality_overall.pdf \
	    OUTPUT=$(basename $BM .markedDup.bam)_read_quality_overall.txt \
	    REFERENCE_SEQUENCE=../${refDir}/$EC_ref

	$picard CollectWgsMetrics \
	    INPUT=$BM \
	    OUTPUT=$(basename $BM .markedDup.bam)_stats_picard.txt \
	    REFERENCE_SEQUENCE=../${refDir}/$EC_ref \
	    MINIMUM_MAPPING_QUALITY=20 \
	    MINIMUM_BASE_QUALITY=20
done

for BM in `ls ECL-*.markedDup.bam`; do
        $picard MeanQualityByCycle \
	    INPUT=$BM \
	    CHART_OUTPUT=$(basename $BM .markedDup.bam)_mean_quality_by_cycle.pdf \
	    OUTPUT=$(basename $BM .markedDup.bam)_read_quality_by_cycle.txt \
	    REFERENCE_SEQUENCE=../${refDir}/$ECL_ref

        $picard QualityScoreDistribution \
	    INPUT=$BM \
	    CHART_OUTPUT=$(basename $BM .markedDup.bam)_mean_quality_overall.pdf \
	    OUTPUT=$(basename $BM .markedDup.bam)_read_quality_overall.txt \
	    REFERENCE_SEQUENCE=../${refDir}/$ECL_ref

	$picard CollectWgsMetrics \
	    INPUT=$BM \
	    OUTPUT=$(basename $BM .markedDup.bam)_stats_picard.txt \
	    REFERENCE_SEQUENCE=../${refDir}/$ECL_ref \
	    MINIMUM_MAPPING_QUALITY=20 \
	    MINIMUM_BASE_QUALITY=20
done

#############################################################
## GATK Data Pre-Processing

# Step 1 - Local realignment around indels.
# Create a target list of intervals to be realigned.

echo "Creating a target list of intervals to be realigned...."

for BM in `ls KP-*.markedDup.bam`; do

samtools index $BM

$gatk \
-T RealignerTargetCreator \
-R ../${refDir}/$KP_ref \
-I $BM \
-o $(basename $BM .markedDup.bam)_target_intervals.list

# do the local realignment.
echo "local realignment..."

$gatk \
-T IndelRealigner \
-R ../${refDir}/$KP_ref \
-I $BM \
-targetIntervals $(basename $BM .markedDup.bam)_target_intervals.list \
-o $(basename $BM .markedDup.bam)_realigned.bam

done 

for BM in `ls EC-*.markedDup.bam`; do

samtools index $BM

$gatk \
-T RealignerTargetCreator \
-R ../${refDir}/$EC_ref \
-I $BM \
-o $(basename $BM .markedDup.bam)_target_intervals.list

# do the local realignment.
echo "local realignment..."

$gatk \
-T IndelRealigner \
-R ../${refDir}/$EC_ref \
-I $BM \
-targetIntervals $(basename $BM .markedDup.bam)_target_intervals.list \
-o $(basename $BM .markedDup.bam)_realigned.bam

done 

for BM in `ls ECL-*.markedDup.bam`; do

samtools index $BM

$gatk \
-T RealignerTargetCreator \
-R ../${refDir}/$ECL_ref \
-I $BM \
-o $(basename $BM .markedDup.bam)_target_intervals.list

# do the local realignment.
echo "local realignment..."

$gatk \
-T IndelRealigner \
-R ../${refDir}/$ECL_ref \
-I $BM \
-targetIntervals $(basename $BM .markedDup.bam)_target_intervals.list \
-o $(basename $BM .markedDup.bam)_realigned.bam

done 

###########################################################################
# GATK Variant Calling -  HaplotypeCaller
# Set -nct, outmode, emit_thresh, call_threh,

outmode="EMIT_ALL_CONFIDENT_SITES"
hetrate=0.03	#Popgen heterozygosity rate (that is, for any two random chrom in pop, what is rate of mismatch).               Human is ~0.01, so up maize to ~0.03
minBaseScore=20	#Minimum Phred base score to count a base (20 = 0.01 error, 30=0.001 error, etc)

for BM in `ls KP-*.markedDup.bam`; do

echo "indexing the realigned bam file..."

# Create a new index file.
samtools index $(basename $BM .markedDup.bam)_realigned.bam

echo "calling variants...."

$gatk \
-T HaplotypeCaller \
-R ../${refDir}/$KP_ref \
-I $(basename $BM .markedDup.bam)_realigned.bam \
--emitRefConfidence GVCF \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-hets $hetrate \
-mbq $minBaseScore \
-out_mode $outmode \
-nct 1 \
-ploidy 1 \
-o $(basename $BM .markedDup.bam)_output.snps.indels.g.vcf

done

for BM in `ls EC-*.markedDup.bam`; do

echo "indexing the realigned bam file..."

# Create a new index file.
samtools index $(basename $BM .markedDup.bam)_realigned.bam

echo "calling variants...."

$gatk \
-T HaplotypeCaller \
-R ../${refDir}/$EC_ref \
-I $(basename $BM .markedDup.bam)_realigned.bam \
--emitRefConfidence GVCF \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-hets $hetrate \
-mbq $minBaseScore \
-out_mode $outmode \
-nct 1 \
-nt $nThreads \
-ploidy 1 \
-o $(basename $BM .markedDup.bam)_output.snps.indels.g.vcf

done

for BM in `ls ECL-*.markedDup.bam`; do

echo "indexing the realigned bam file..."

# Create a new index file.
samtools index $(basename $BM .markedDup.bam)_realigned.bam

echo "calling variants...."

$gatk \
-T HaplotypeCaller \
-R ../${refDir}/$ECL_ref \
-I $(basename $BM .markedDup.bam)_realigned.bam \
--emitRefConfidence GVCF \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-hets $hetrate \
-mbq $minBaseScore \
-out_mode $outmode \
-nct 1 \
-nct $nThreads \
-ploidy 1 \
-o $(basename $BM .markedDup.bam)_output.snps.indels.g.vcf

done

echo ""
echo "=============================================================================="
echo " Finished generating files"
echo "=============================================================================="
echo ""

exit
