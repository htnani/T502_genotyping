#!/bin/bash

#PBS -N WGS_genotyping_align_T502
#PBS -k o
#PBS -l nodes=1:ppn=16,vmem=40gb
#PBS -l walltime=16:00:00
#PBS -m abe

module load bwa
module load samtools

fileDir=/N/dc2/scratch/rtraborn/T502_fastqs/Genotyping_Ecoli_new
####### Before running the script, please enter path to desired output directory, below ####
WD=/N/dc2/scratch/rtraborn/T502_genotyping/
#WD=<enter/path/to/T502_genotyping/>
fqDir=fastqs
genomedir=${WD}/fasta
KPgenome=${genomedir}/K_pneumoniae_genome.fasta
ECIgenome=${genomedir}/E_coli_genome.fasta
ECgenome=${genomedir}/E_cloacae_genomic.fasta
nThreads=16


echo "Starting job"

echo "Making symbolic links to fastq files."

cd $WD

if [ ! -d "$fqDir" ] ; then
    mkdir $fqDir
  fi

cd $fqDir

ln -s ${fileDir}/GSF1659-EC-7_S17_R1_001.fastq EC-7_S17_R1.fastq
ln -s ${fileDir}/GSF1659-EC-7_S17_R2_001.fastq EC-7_S17_R2.fastq
ln -s ${fileDir}/GSF1659-ECL-1_S15_R1_001.fastq ECL-1_S15_R1.fastq
ln -s ${fileDir}/GSF1659-ECL-1_S15_R2_001.fastq ECL-1_S15_R2.fastq
ln -s ${fileDir}/GSF1659-ECL-4_S16_R1_001.fastq ECL-4_S16_R1.fastq
ln -s ${fileDir}/GSF1659-ECL-4_S16_R2_001.fastq ECL-4_S16_R2.fastq
ln -s ${fileDir}/GSF1659-KP-177_S18_R1_001.fastq KP-177_S18_R1.fastq
ln -s ${fileDir}/GSF1659-KP-177_S18_R2_001.fastq KP-177_S18_R2.fastq
ln -s ${fileDir}/GSF1659-KP-49_S9_R1_001.fastq KP-49_S9_R1.fastq
ln -s ${fileDir}/GSF1659-KP-49_S9_R2_001.fastq KP-49_S9_R2.fastq
ln -s ${fileDir}/GSF1659-KP-80_S10_R1_001.fastq KP-80_S10_R1.fastq
ln -s ${fileDir}/GSF1659-KP-80_S10_R2_001.fastq KP-80_S10_R2_001.fastq
ln -s ${fileDir}/GSF1659-KP-83_S11_R1_001.fastq KP-83_S11_R1.fastq
ln -s ${fileDir}/GSF1659-KP-83_S11_R2_001.fastq KP-83_S11_R2.fastq
ln -s ${fileDir}/GSF1659-KP-84_S12_R1_001.fastq KP-84_S12_R1.fastq
ln -s ${fileDir}/GSF1659-KP-84_S12_R2_001.fastq KP-84_S12_R2.fastq
ln -s ${fileDir}/GSF1659-KP-86_S14_R1_001.fastq KP-86_S14_R1.fastq
ln -s ${fileDir}/GSF1659-KP-86_S14_R2_001.fastq KP-86_S14_R2.fastq
ln -s ${fileDir}/GSF1659-KP-85_S13_R1_001.fastq KP-85_S13_R1.fastq
ln -s ${fileDir}/GSF1659-KP-85_S13_R2_001.fastq KP-85_S13_R2.fastq

cd $WD

#echo "Retrieving the appropriate assembly files"

#source 0README

echo "Performing alignment on the bacterial files"

cd $genomedir
bwa index $KPgenome

cd ${WD}/${fqDir}
  for file1 in `ls KP-*_R1.fastq`; do
    echo $file1
    file2=$(basename $file1 _R1.fastq)_R2.fastq
    echo "bwa mem $KPgenome -t $nThreads $file1 $file2 | samtools view > $(basename $file1 _R1.fastq).bam"
    bwa mem -t $nThreads $KPgenome $file1 $file2 | samtools view -h > $(basename $file1 _R1.fastq).bam
  done

cd $genomedir
bwa index $ECIgenome

cd ${WD}/${fqDir}
  for file1 in `ls ECL-*_R1.fastq`; do
    file2=$(basename $file1 _R1.fastq)_R2.fastq
    echo "bwa mem -t $nThreads $KPgenome $file1 $file2 | samtools view > $(basename $file1 _R1.fastq).bam"
    bwa mem -t $nThreads $ECLgenome $file1 $file2 | samtools view -h > $(basename $file1 _R1.fastq).bam
  done

cd $genomedir
bwa index $ECgenome

cd ${WD}/${fqDir}
  for file1 in `ls EC-*_R1.fastq`; do
    file2=$(basename $file1 _R1.fastq)_R2.fastq
    echo "bwa mem -t $nThreads $KPgenome $file1 $file2 | samtools view > $(basename $file1 _R1.fastq).bam"
    bwa mem -t $nThreads $ECgenome $file1 $file2 | samtools view -bh > $(basename $file1 _R1.fastq).bam
  done

  cd ${WD}
  if [ ! -d "alignments" ] ; then
    mkdir alignments
  fi
  cd alignments
  for file1 in ${WD}/$fqDir/*.bam ; do
    echo $file1
    file2=$(basename $file1)
    echo $file2
    ln -s $file1 ./$file2
  done
  cd ..

  echo ""
  echo " Done with step 3 (read mapping)."
  echo ""
  echo "================================================================================"

exit


