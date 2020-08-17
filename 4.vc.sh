#!/bin/bash


#Location of the tools
bwaMemLoc=/data/sata_data/workshop/wsu28/packages/bwa/bwa
samtools=/data/sata_data/workshop/wsu28/packages/samtools1/samtools/samtools
bcftools=/data/sata_data/workshop/wsu28/packages/samtools1/bcftools/bcftools


#Mosquito Reference genome location
refFasta=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/test1/finalRefGen.fa



#Location of the raw reads
Read1=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simu_genome204/Read1.fq
Read2=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simu_genome204/Read2.fq
Read3=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simu_genome204/Read3.fq
Read4=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simu_genome204/Read4.fq





#Where the result will be stored.
presentwd=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/test204
cd $presentwd

#$bwaMemLoc index $refFasta
#echo "Index completed"

#Now mapping all those files
$bwaMemLoc mem -t 20 -M -R "@RG\tID:group1\tSM:sample1\tPL:Illumina\tLB:lib1\tPU:unit1" $refFasta $Read1 $Read2 > mappedR1R2.sam
$bwaMemLoc mem -t 20 -M -R "@RG\tID:group2\tSM:sample1\tPL:Illumina\tLB:lib1\tPU:unit1" $refFasta $Read3 $Read4 > mappedR3R4.sam
echo "Mapping completed"

#Convert SAM to BAM
$samtools view -S -b mappedR1R2.sam > mappedR1R2.bam -@ 20
$samtools view -S -b mappedR3R4.sam > mappedR3R4.bam -@ 20
echo "Sam to Bam compleated"

#Sort BAM files
$samtools sort mappedR1R2.bam > mappedR1R2_sorted.bam
$samtools sort mappedR1R2.bam > mappedR3R4_sorted.bam
echo "Sort compleated"


#Now gather all the BAM files.
$samtools merge finalBamFile.bam mappedR*R*_sorted.bam


#Samtools indexing is there.
#$samtools faidx $refFasta
#echo "Samtool index completed"


#####################################################################################
#$samtools mpileup -g -f $refFasta mappedR1R2_sorted.bam > mappedR1R2_sorted_raw.bcf
#echo "mpileup cpmpleated"

#$bcftools call -O b -vc mappedR1R2_sorted_raw.bcf > mappedR1R2_sorted_var.bcf
#echo "Call SNP compleated"
#####################################################################################


#Use this one insted of using upper two line.
$bcftools mpileup -f $refFasta finalBamFile.bam | $bcftools call -mv -Ob -o mappedR1234_sorted_var.bcf
echo "Mpileup and Call SNP compleated"


varFilter=/data/sata_data/workshop/wsu28/packages/samtools1/bcftools/misc
$bcftools view mappedR1234_sorted_var.bcf | $varFilter/vcfutils.pl varFilter -Q > MappedR1234.var-final.vcf
echo "Final vcf compleated"



######################################################################################
#$bcftools view -v snps Result_GCF_genomic.refseq2simseq.SNP.vcf | grep -v "^#" | wc -l
#281346

#Filtering variants
#$bcftools view -i '%QUAL>=20' mappedR1R2_sorted_var1.bcf
