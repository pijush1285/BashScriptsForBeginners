read -p " enter your Reference:" Re1
read -p " enter the Read 1:" r1
read -p " enter the Read 2:" r2


Re1=finalRefGen.fa
r1=NRead3.fq
r2=NRead4.fq



#Location of the tools
bwa=/data/sata_data/workshop/wsu28/packages/bwa/bwa
samtools=/data/sata_data/workshop/wsu28/packages/samtools1/samtools/samtools
bcftools=/data/sata_data/workshop/wsu28/packages/samtools1/bcftools/bcftools
picard=/data/sata_data/workshop/wsu28/packages/picard/build/libs/picard.jar
gatk=/data/sata_data/workshop/wsu28/packages/gatk/gatk





#Indexing
$bwa index "$Re1"


#ALignment 
$bwa mem -M -t 20 "$Re1" "$r1" "$r2" > "$r1"_"$r2"_bwa.sam


#Bam_File #Running
$samtools view -S -b "$r1"_"$r2"_bwa.sam > "$r1"_"$r2"_bwa.bam -@ 20


#Sort_Bam 
$samtools sort "$r1"_"$r2"_bwa.bam > "$r1"_"$r2"_Sort_bwa.bam


#Marked Duplicates  
java -jar $picard MarkDuplicates \
      I="$r1"_"$r2"_Sort_bwa.bam \
      O="$r1"_"$r2"_marked_duplicates.bam \
      M="$r1"_"$r2"_marked_dup_metrics.txt


#Add replace Group   
java -jar $picard AddOrReplaceReadGroups \
       I="$r1"_"$r2"_marked_duplicates.bam \
       O="$r1"_"$r2"_RG_bwa.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=20


#Now creating the index file. I am not sure this code is needed or not.
#java -jar $picard BuildBamIndex \
#      -I "$r1"_"$r2"_RG_bwa.bam \
#      -O "$r1"_"$r2"_RG_bwa.bai \
#      >"$r1"_"$r2"_RG_bwa_bamindex.out 2>"$r1"_"$r2"_RG_bwa_bamindex.err


#Creating the index file   
jAdd replace Groupava -jar $picard BuildBamIndex \
      I="$r1"_"$r2"_RG_bwa.bam


#$samtools faidx "$Re1"


#Index file generate of the reference file.
java -jar $picard CreateSequenceDictionary \ 
      R="$Re1" \ 
      O=finalRefGen.dict


#Variant call  #Running   
$gatk --java-options "-Xmx10g" HaplotypeCaller -R "$Re1" -I "$r1"_"$r2"_RG_bwa.bam -O "$r1"_BWA_GATK_"$r2".vcf.gz


#Merge VCF files. This code is working here.
java -jar $picard MergeVcfs I=NRead1.fq_BWA_GATK_NRead2.fq.vcf.gz I=NRead3.fq_BWA_GATK_NRead4.fq.vcf.gz O=Final_merged_variants.vcf.gz

#The code removes the duplicates.
bcftools norm -D Final_merged_variants.vcf.gz -o Final_merged_variants1.vcf.gz


#################################################################################################################


#variant_Sepration_Indel_SNV
#vcftools --vcf output/"$r1"/BWA_GATK_"$r1".vcf --remove-indels --recode --recode-INFO-all --out output/"$r1"/BWA_GATK_SNP_"$r1".vcf
#vcftools --vcf output/"$r1"/BWA_GATK_"$r1".vcf --keep-only-indels  --recode --recode-INFO-all --out output/"$r1"/BWA_GATK_Indels_"$r1".vcf


#################################################################################################################



bcftools merge NRead*.fq_BWA_GATK_NRead*.fq.vcf.gz -Oz -o Final_merged1234_variants.vcf.gz


#Total no of SNP
bcftools view -v snps Final_merged_variants1.vcf | grep -v "^#" | wc -l
#948172
bcftools view -v snps NRead1.fq_BWA_GATK_NRead2.fq.vcf.gz | grep -v "^#" | wc -l
#944962
bcftools view -v snps NRead3.fq_BWA_GATK_NRead4.fq.vcf.gz | grep -v "^#" | wc -l
#945012



#This code is also required.
bcftools view -v snps Final_merged_variants1.vcf | perl -lane 'if (/^#/) { print } elsif (length($F[3]) == 1) { print }' | bgzip > Final_merged_variants1.vcf.gz


# index
tabix -p vcf Final_merged_variants1.vcf.gz

# calculate Jaccard index
bedtools jaccard -a snpOriginal.vcf.gz -b Final_merged_variants1.vcf.gz

#intersection	union	jaccard	n_intersections
#948053	1000119	0.94794	947307

