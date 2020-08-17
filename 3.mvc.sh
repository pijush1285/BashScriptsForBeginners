#!/bin/bash


#Location of the tools
bwaMemLoc=/data/sata_data/workshop/wsu28/packages/bwa/bwa
samtools=/data/sata_data/workshop/wsu28/packages/samtools1/samtools/samtools
bcftools=/data/sata_data/workshop/wsu28/packages/samtools1/bcftools/bcftools


#Mosquito Reference genome location
refFasta=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/ncbi-genomes-2020-06-18/GCF_genomic.fa


#Location of the raw reads
#simuRead=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/Result_GCF_genomic.simseq.genome.fa
simuRead1=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simulated_reads/GCF_GPI_1000SNP.read1.fq
simuRead2=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simulated_reads/GCF_GPI_1000SNP.read2.fq


#Location to save the output results
result=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/mapped


#1) Align reads to reference (using BWA)

#1. Index the reference (genome) sequence
#$bwaMemLoc index $refFasta

#2. Perform the alignment
#Number of thread used
thread=20

cd $result
#$bwaMemLoc mem $refFasta $simuRead1 $simuRead2 -t $thread > simu.out.sam


#accuracy == "strict"
#$bwaMemLoc mem -k 16 -W 60 $refFasta $simuRead1 $simuRead2 -t $thread > simu_Strick.out.sam



#1. Convert SAM to BAM for sorting
#$samtools view -S -b $result/simu_Strick.out.sam > simu_out.bam


#2. Sort BAM for SNP calling
#cd $result¢
#$samtools sort $result/simu_out.bam > simu_sorted.bam


#2) Call SNPs (using SAMtools)
#1. Index the genome assembly (again!)
locaofref=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/ncbi-genomes-2020-06-18
#cd $locaofref
#$samtools faidx $refFasta
#echo "samtools faidx completed."


#2. Run 'mpileup' to generate VCF format¢
cd $result
#$samtools mpileup -g -f $refFasta $result/simu_sorted.bam > file.mpileup.bcf
#echo "mileup completed."



#3. Call SNPs...
#$bcftools view -bcg $result/file.mpileup.bcf > my-var.bcf

#$bcftools mpileup -f $refFasta $result/simu_sorted.bam | $bcftools call -mv -Ob -o calls.bcf
#echo "Call SNPs completed."



#3) Filter SNPs
#1. Filter SNPs
varFilter=/data/sata_data/workshop/wsu28/packages/samtools1/bcftools/misc
$bcftools view calls.bcf | $varFilter/vcfutils.pl varFilter -p > my.var-final.vcf


##########################################################################
#$samtools mpileup -E -uf $refFasta $result/simu_sorted.bam > file.mpileup
#echo "mileup completed."


#$bcftools view -cg file.mpileup > file.vcf
#echo "Call SNPs completed.


#$bcftools view -bcg file.mpileup > file.bcf
#echo "2nd Call SNPs completed.


#samtools mpileup -u -l chr1.snps -f hg19.fa sample1_filter_sort.bam | bcftools view -v -c - > out.vcf


#$bcftools view -i '%QUAL>=20' calls.bcf > file.bcf
#$bcftools query -f '%MyAnnotation\n' calls.bcf | my-plotting-program
#$bcftools filter -i'%QUAL>20' calls.bcf | $bcftools stats | grep TSTV
echo "Filter completed."
