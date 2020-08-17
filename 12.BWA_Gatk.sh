read -p " enter your Reference:" Re1
read -p " enter the Read 1:" r1
read -p " enter the Read 2:" r2


Re1=finalRefGen.fa
r1=NRead1.fq
r2=NRead2.fq



#Location of the tools
bwa=/data/sata_data/workshop/wsu28/packages/bwa/bwa
samtools=/data/sata_data/workshop/wsu28/packages/samtools1/samtools/samtools
bcftools=/data/sata_data/workshop/wsu28/packages/samtools1/bcftools/bcftools
picard=/data/sata_data/workshop/wsu28/packages/picard/build/libs/picard.jar
gatk=/data/sata_data/workshop/wsu28/packages/gatk/gatk




#indexing
$bwa index "$Re1"
#ALignment
$bwa mem -M -t 20 "$Re1" "$r1" "$r2" > "$r1"_"$r2"_bwa.sam 
#Bam_File
$samtools view -bS "$r1"_"$r2"_bwa.sam -o "$r1"_"$r2"_bwa.bam   
#Sort_Bam  #Running . This sort code is not working. 
java -Xmx10g -jar $picard SortSam VALIDATION_STRINGENCY=SILENT I="$r1"_"$r2"_bwa.bam O="$r1"_"$r2"_Sort_bwa.bam SORT_ORDER=coordinate
#Pcr_Duplicates
java -Xmx10g -jar $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT I="$r1"_"$r2"_Sort_bwa.bam O="$r1"_"$r2"_PCR_bwa.bam REMOVE_DUPLICATES=true M="$r1"_"$r2"_pcr_bwa.metrics
#ID_Addition
java -Xmx10g -jar $picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT I="$r1"_"$r2"_PCR_bwa.bam O="$r1"_"$r2"_RG_bwa.bam SO=coordinate RGID=SRR"$r1" RGLB=SRR"$r2" RGPL=illumina RGPU=SRR"$r1" RGSM=SRR"$r2" CREATE_INDEX=true
#Variant_calling
$gatk --java-options "-Xmx10g" HaplotypeCaller -R "$Re1" -I "$r1"_"$r2"_RG_bwa.bam -O "$r1"_BWA_GATK_"$r2".vcf.gz

#variant_Sepration_Indel_SNV
#vcftools --vcf output/"$r1"/BWA_GATK_"$r1".vcf --remove-indels --recode --recode-INFO-all --out output/"$r1"/BWA_GATK_SNP_"$r1".vcf
#vcftools --vcf output/"$r1"/BWA_GATK_"$r1".vcf --keep-only-indels  --recode --recode-INFO-all --out output/"$r1"/BWA_GATK_Indels_"$r1".vcf



##################################################################################################################
#Sort SAM
#java -Dpicard.useLegacyParser=false -jar $picard SortSam -I "$r1"_"$r2"_bwa.sam -SORT_ORDER coordinate -O sorted.sam
#java -Dpicard.useLegacyParser=false -jar $picard SortSam -I "$r1"_"$r2"_bwa.bam -SORT_ORDER coordinate -O "$r1"_"$r2"_Sort_bwa.bam -TMP_DIR /data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/test204/working_temp




################################################################################################

#Indexing
$bwa index "$Re1"


#ALignment
$bwa mem -M -t 20 "$Re1" "$r1" "$r2" > "$r1"_"$r2"_bwa.sam


#Bam_File
$samtools view -S -b "$r1"_"$r2"_bwa.sam > "$r1"_"$r2"_bwa.bam -@ 20


#Sort_Bam  
$samtools sort "$r1"_"$r2"_bwa.bam > "$r1"_"$r2"_Sort_bwa.bam


#Marked Duplicates
java -jar $picard MarkDuplicates \
      I="$r1"_"$r2"_Sort_bwa.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt


#Add replace Group
java -jar $picard AddOrReplaceReadGroups \
       I=marked_duplicates.bam \
       O="$r1"_"$r2"_RG_bwa.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=20


#Now creating the index file. I am not sure this code is needed or not.
java -jar $picard BuildBamIndex \
      -I "$r1"_"$r2"_RG_bwa.bam \
      -O "$r1"_"$r2"_RG_bwa.bai \
      >"$r1"_"$r2"_RG_bwa_bamindex.out 2>"$r1"_"$r2"_RG_bwa_bamindex.err


#java -jar $picard BuildBamIndex \
#      I="$r1"_"$r2"_RG_bwa.bam

#$samtools faidx "$Re1"


#Index file generate of the reference file.
java -jar $picard CreateSequenceDictionary \ 
      R="$Re1" \ 
      O=finalRefGen.dict


#Variant call
$gatk --java-options "-Xmx10g" HaplotypeCaller -R "$Re1" -I "$r1"_"$r2"_RG_bwa.bam -O "$r1"_BWA_GATK_"$r2".vcf.gz


#################################################################################################################





$gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $Re1 \
   -I "$r1"_"$r2"_RG_bwa.bam \
   -O "$r1"_BWA_GATK_"$r2".vcf.gz \
   -ERC GVCF \
   -G Standard \
   -G AS_Standard






java -jar $picard SortSam -I "$r1"_"$r2"_bwa.bam -O "$r1"_"$r2"_Sort_bwa.bam -SORT_ORDER coordinate

java -jar $picard SortSam \
      I="$r1"_"$r2"_bwa.bam \
      O="$r1"_"$r2"_Sort_bwa.bam \
      SORT_ORDER=coordinate \
      TMP_DIR=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/test204/working_temp



"java -jar $PICARD_TOOLS_DIR/picard.jar MarkDuplicates INPUT=output.bam 
OUTPUT=output.marked.bam METRICS_FILE=metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/short/a32/working_temp"



#SAM to BAM
$samtools view -S -b sorted.sam > mappedR1R2_sorted.bam -@ 20

#Marked Duplicate
java -jar picard.jar MarkDuplicates \
      I=input.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt



