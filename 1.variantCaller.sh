#!/bin/bash

#only for LONG READS mapping
#USAGE: runMapper.sh minimap2 ref.fa reads.fa/fq 48 ont default/strict

echo -e "-----------------------------------------------------------------------------"
echo -e "                                                                             "
echo -e "This is a bash scrip to map your long reads and make it visualization ready !"
echo -e "                                                                             "
echo -e "-----------------------------------------------------------------------------"

#################################################################################################################
#For transfer the program file from server to local machine
#scp wsu28@10.10.4.201:/data/sata_data/workshop/wsu28/variantCaller/variantCaller.sh /home/icgc-main/Desktop/asdf

#For transfer the program file from local machine to server
#scp /home/icgc-main/Desktop/asdf/variantCaller.sh wsu28@10.10.4.201:/data/sata_data/workshop/wsu28/variantCaller
#################################################################################################################

######################################################################################
#            This portion of code is used for alignment and fastQC
#####################################################################################
#Location of the tools
bwaMemLoc=/data/sata_data/workshop/wsu28/packages/bwa/bwa
fastqc=/data/sata_data/workshop/wsu28/packages/FastQC/fastqc
samtools=/data/sata_data/workshop/wsu28/packages/samtools1/samtools/samtools
picard=/data/sata_data/workshop/wsu28/packages/picard/picard.jar
gatk=/data/sata_data/workshop/wsu28/packages/gatk/gatk


#Reference genome location
refFasta=/data/sata_data/workshop/wsu28/reference/hg38.fa


#For concatinating the fasta files used the code given below.
#cat U0a_CGATGT_L001_R1_00*.fastq > U0a_CGATGT_L001_R1_Marged.fastq
#cat U0a_CGATGT_L001_R2_00*.fastq > U0a_CGATGT_L001_R2_Marged.fastq


#Location of the raw reads
longReadsL1R1S1=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_001.fastq
longReadsL1R1S2=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_002.fastq
longReadsL1R1S3=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_003.fastq
longReadsL1R1S4=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_004.fastq
longReadsL1R1S5=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_005.fastq

longReadsL1R2S1=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_001.fastq
longReadsL1R2S2=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_002.fastq
longReadsL1R2S3=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_003.fastq
longReadsL1R2S4=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_004.fastq
longReadsL1R2S5=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_005.fastq


######################################################################
#                               FastQC
######################################################################
#The output file will be assembled under a directory
OutPutPath="/data/sata_data/workshop/wsu28/vcResults/"
samplelabel="Sample_U0a"
FQoutDir="FastQC"
MapFileContainer="MapFileContainer"


#Frist create Sample output folder
if [ ! -d $OutPutPath$samplelabel ]; then
  mkdir -p $OutPutPath$samplelabel;
fi

#Make a FastQc folder to save the FastQc results
if [ ! -d $OutPutPath$samplelabel$FQoutDir ]; then
  mkdir -p $OutPutPath$samplelabel/$FQoutDir;
fi

#Now change the directory.
#cd /data/sata_data/workshop/wsu28/vcResults/FastQC
#OutPutDir=$OutPutPath$samplelabel/$FQoutDir


#Doing the FastQC
#$fastqc $longReadsL1R1S1 -o $OutPutDir
#$fastqc $longReadsL1R1S2 -o $OutPutDir
#$fastqc $longReadsL1R1S3 -o $OutPutDir
#$fastqc $longReadsL1R1S4 -o $OutPutDir
#$fastqc $longReadsL1R1S5 -o $OutPutDir

#$fastqc $longReadsL1R2S1 -o $OutPutDir
#$fastqc $longReadsL1R2S2 -o $OutPutDir
#$fastqc $longReadsL1R2S3 -o $OutPutDir
#$fastqc $longReadsL1R2S4 -o $OutPutDir
#$fastqc $longReadsL1R2S5 -o $OutPutDir



#logPaths=("vv" "$longReadsL1R1S1" "$longReadsL1R1S2" "$longReadsL1R1S3" "$longReadsL1R1S4" "$longReadsL1R1S5")
#logPaths1=("vv" "$longReadsL1R2S1" "$longReadsL1R2S2" "$longReadsL1R2S3" "$longReadsL1R2S4" "$longReadsL1R2S5")


#for i in {1..6}
#do
  # spl=$(longReadsL1R1S${i})
#   $fastqc ${logPaths[$i]} -o $OutPutDir
#   $fastqc ${logPaths1[$i]} -o $OutPutDir
#   echo "Processing complete: $i times"
#done


###################################################################################
#              Three type of indexing using different tools
###################################################################################

#1 BWA Indexing
#$bwaMemLoc index $refFasta


#2 Index by using samtools
#$samtools faidx $refFasta


#3 Index by using picard
#java -jar $picard CreateSequenceDictionary \
#REFERENCE=$refFasta \
#OUTPUT=hg38.dict


##################################################################################
#			        	Mapping
##################################################################################

#Make a MapFileContainer folder to save the mapped sam files results
if [ ! -d $OutPutPath$samplelabel$FQoutDir ]; then
  mkdir -p $OutPutPath$samplelabel/$MapFileContainer;
fi

#Now change the directory.
OutPutDir1=$OutPutPath$samplelabel/$MapFileContainer



#Output File name
fileName=U0a_CGATGT_Default_marged

#Number of thread used
thread=20


#accuracy == "default"
#cd $OutPutDir1
#$bwaMemLoc mem $refFasta $longReadsL1R1S1 $longReadsL1R2S1 -t $thread > SplL1R1R2S1.out.sam
#$bwaMemLoc mem $refFasta $longReadsL1R1S2 $longReadsL1R2S2 -t $thread > SplL1R1R2S2.out.sam
#$bwaMemLoc mem $refFasta $longReadsL1R1S3 $longReadsL1R2S3 -t $thread > SplL1R1R2S3.out.sam
#$bwaMemLoc mem $refFasta $longReadsL1R1S4 $longReadsL1R2S4 -t $thread > SplL1R1R2S4.out.sam
#$bwaMemLoc mem $refFasta $longReadsL1R1S5 $longReadsL1R2S5 -t $thread > SplL1R1R2S5.out.sam





#accuracy == "strict"
#$bwaMemLoc mem -k 16 -W 60 $refFasta $longReads1 $longReads2 -t $thread > $fileName.out.sam



############################################################################################
#                        AddOrReplaceReadGroups (Picard)
############################################################################################


#logPaths=("vv" "SplL1R1R2S1.out.sam" "SplL1R1R2S2.out.sam" "SplL1R1R2S3.out.sam" "SplL1R1R2S4.out.sam" "SplL1R1R2S5.out.sam")
#logPaths1=("vv" "SplL1R1R2S1.bam" "SplL1R1R2S2.bam" "SplL1R1R2S3.bam" "SplL1R1R2S4.bam" "SplL1R1R2S5.bam")

#for i in {1..6}
#  do
    #Input file
#    input=$OutPutDir1/${logPaths[$i]}
#    echo "----------------------------------"
#    echo $input
#    echo "Input file name and location"
#    #Output file
#    output=$OutPutDir1/${logPaths1[$i]}
#    echo "----------------------------------"
#    echo $output
#    echo "Output file name and location"
#
#    java -jar $picard AddOrReplaceReadGroups \
#          I=$input \
#          O=$output \
#          RGID=4 \
#          RGLB=lib1 \
#          RGPL=ILLUMINA \
#          RGPU=unit1 \
#          RGSM=20

#          echo "Processing complete: $i times"

#done


#############################################################################################
#                 Base Quality Score Recalibration (Create table)
#############################################################################################
#Make a BQSRcreateTable folder to save the mapped sam files results
BQSRcreateTable="BQSRcreateTable"
if [ ! -d $OutPutPath$samplelabel$FQoutDir ]; then
  mkdir -p $OutPutPath$samplelabel/$BQSRcreateTable;
fi

#New Directory location
OutPutDir2=$OutPutPath$samplelabel/$BQSRcreateTable


#Before using the vcf file you would required to make a tab index vcf file
#tabix -p vcf 1000G_phase1.snps.high_confidence.hg38.vcf.gz

goldRndelsVCF=/data/sata_data/workshop/wsu28/vcf1/1000G_phase1.snps.high_confidence.hg38.vcf.gz
setOfSitesToMask=/data/sata_data/workshop/wsu28/vcf1/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

#Input file names for this function
logPaths1=("vv" "SplL1R1R2S1.bam" "SplL1R1R2S2.bam" "SplL1R1R2S3.bam" "SplL1R1R2S4.bam" "SplL1R1R2S5.bam")


#Location of input files (This is the input file)
dedupReads=$OutPutDir1/SplL1R1R2S1.bam

#Output File name
#targetIntervalsList=/data/sata_data/workshop/wsu28/markedDuplicates/recal_data.table.txt

#Change the directory for the output
cd OutPutDir2
#Now using the base recalibration
$gatk BaseRecalibrator \
   -I $dedupReads \
   -R $refFasta \
   --known-sites $goldRndelsVCF \
   --known-sites $setOfSitesToMask \
   -O recal_data.table.txt
   >$mysamplebase-recalibrate.out 2>$mysamplebase-recalibrate.err




#Base Quality Score Recalibration (apply)

$gatk PrintReads \
    -R $refFasta \
    -I $dedupReads \
    -BQSR $mysamplebase"-recal-table.txt" \
    -o  $mysamplebase"-GATK.bam" \
    > $mysamplebase-recalbam.out 2> $mysamplebase-recalbam.err
