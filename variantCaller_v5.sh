#!/bin/bash

echo -e "-----------------------------------------------------------------------------"
echo -e "                                                                             "
echo -e "This script will perform the variant calling operation step by step !"
echo -e "This program is developed for the analysis of the lab made dataset....       "
echo -e "In this program I am taking some reference from Abhishek..............       "
echo -e "                                                                             "
echo -e "-----------------------------------------------------------------------------"



######################################################################################
#            This portion of code is used for alignment and fastQC
#####################################################################################
#Location of the tools
bwaMemLoc=/data/sata_data/workshop/wsu28/packages/bwa/bwa
fastqc=/data/sata_data/workshop/wsu28/packages/FastQC/fastqc
samtools=/data/sata_data/workshop/wsu28/packages/samtools1/samtools/samtools
bamtools=/data/sata_data/workshop/wsu28/anaconda3/envs/pdas/bin/bamtools
picard=/data/sata_data/workshop/wsu28/packages/picard/picard.jar
gatk=/data/sata_data/workshop/wsu28/packages/gatk/gatk
bbduk=/data/sata_data/workshop/wsu28/anaconda3/bin/bbduk.sh


#Reference genome location
refFasta=/data/sata_data/workshop/wsu28/reference/hg38.fa

#Delete the folder if any
# OutPutPath="/data/sata_data/workshop/wsu28/vcResults"
# cd OutPutPath
# rm -r 14810642_v5
# echo "There is no such folder named 14810642_v5"

##############################################################################################
#                                 Location of the raw reads
##############################################################################################


longReadsL1R1=/data/sata_data/workshop/wsu28/rawdata/DATA/14810642/14810642_TruSeqCD_BHW2KKDSXX_L1_R1.fastq.gz
longReadsL1R2=/data/sata_data/workshop/wsu28/rawdata/DATA/14810642/14810642_TruSeqCD_BHW2KKDSXX_L1_R2.fastq.gz

longReadsL2R1=/data/sata_data/workshop/wsu28/rawdata/DATA/14810642/14810642_TruSeqCD_BHW2KKDSXX_L2_R1.fastq.gz
longReadsL2R2=/data/sata_data/workshop/wsu28/rawdata/DATA/14810642/14810642_TruSeqCD_BHW2KKDSXX_L2_R2.fastq.gz

longReadsL3R1=/data/sata_data/workshop/wsu28/rawdata/DATA/14810642/14810642_TruSeqCD_BHW2KKDSXX_L3_R1.fastq.gz
longReadsL3R2=/data/sata_data/workshop/wsu28/rawdata/DATA/14810642/14810642_TruSeqCD_BHW2KKDSXX_L3_R2.fastq.gz

longReadsL4R1=/data/sata_data/workshop/wsu28/rawdata/DATA/14810642/14810642_TruSeqCD_BHW2KKDSXX_L4_R1.fastq.gz
longReadsL4R2=/data/sata_data/workshop/wsu28/rawdata/DATA/14810642/14810642_TruSeqCD_BHW2KKDSXX_L4_R2.fastq.gz


echo "-------------------------------------------"
echo "All the sample files are taken as an Input"
echo "-------------------------------------------"



######################################################################
#                            Cleaning the reads
######################################################################
#The output file will be assembled under a directory
OutPutPath="/data/sata_data/workshop/wsu28/vcResults/"
samplelabel="14810642_v5"

#Frist create Sample output folder
if [ ! -d $OutPutPath$samplelabel ]; then
  mkdir -p $OutPutPath$samplelabel;
fi

#Make a FastQc folder to save the FastQc results
FQoutDir="1_FastQC"
if [ ! -d $OutPutPath$samplelabel/$FQoutDir ]; then
  mkdir -p $OutPutPath$samplelabel/$FQoutDir;
fi

#cd /data/sata_data/workshop/wsu28/vcResults/FastQC
OutPutDir=$OutPutPath$samplelabel/$FQoutDir


# #cleaning reads
# $bbduk in=$longReadsL1R1 \
# in2=$longReadsL1R2 \
# out=$OutPutDir/Cleaned_TruSeqCD_L1_R1.fastq.gz \
# out2=$OutPutDir/Cleaned_TruSeqCD_L1_R2.fastq.gz \
# minavgquality=20


echo "clean read is compleated"



##################################################################################
#			        	Mapping
##################################################################################

#Make a MapFileContainer folder to save the mapped sam files results
MapFileContainer="2_MapFileContainer"
if [ ! -d $OutPutPath$samplelabel/$MapFileContainer ]; then
  mkdir -p $OutPutPath$samplelabel/$MapFileContainer;
fi

#Now change the directory.
OutPutDir1=$OutPutPath$samplelabel/$MapFileContainer


fql1r1=/data/sata_data/workshop/wsu28/vcResults/14810642_v5/1_FastQC/Cleaned_TruSeqCD_L1_R1.fastq.gz
fql1r2=/data/sata_data/workshop/wsu28/vcResults/14810642_v5/1_FastQC/Cleaned_TruSeqCD_L1_R2.fastq.gz


#Alinging reads with the reference

cd $OutPutDir1
echo "BWA running"
#$bwaMemLoc mem -t 20 -M -R "@RG\tID:group1\tSM:sample*\tPL:Illumina\tLB:lib1\tPU:unit1" $refFasta $fql1r1 $fql1r2 > mappedL1R1R2.sam



#Converting sam to bam
echo "Converting Sam to bam"
$samtools view -Sb mappedL1R1R2.sam > mappedL1R1R2.bam


 #samtools sort :
 echo "Sorted bam file production"
 $samtools sort mappedL1R1R2.bam -o mappedL1R1R2.sorted.bam -O bam -@ 20


 #Command to remove duplicate reads
 echo "MarkDuplicates is running"
$gatk MarkDuplicates -I mappedL1R1R2.sorted.bam -O dedup_mappedL1R1R2.bam -M marked_dup.txt
