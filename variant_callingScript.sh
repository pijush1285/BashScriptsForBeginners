

# Downloading data
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz


# downloading and installing BWA tools. Can skip if it is already installed

cd
curl -L https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2/download > bwa-0.7.15.tar.bz2
tar xjvf bwa-0.7.15.tar.bz2
cd bwa-0.7.15
make

sudo cp bwa /usr/local/bin

echo 'export PATH=$PATH:/usr/local/bin' >> ~/.bashrc
source ~/.bashrc

# downloading picard
wget https://github.com/broadinstitute/picard/releases/download/2.21.9/picard.jar

# downloading genome files
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz


gunzip Homo_sapiens.GRCh38.dna.alt.fa.gz

# downloading vcf file
wget ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz

gunzip 1000GENOMES-phase_3.vcf.gz

#cleaning reads
bbduk.sh in=NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
in2=NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
out=Cleaned_NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
out2=Cleaned_NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
minavgquality=20

#Indexing genome using bwa
echo "BWA indexing"
bwa index -a bwtsw Homo_sapiens.GRCh38.dna.alt.fa


#Alinging reads with the reference
echo "BWA running"
bwa mem -t 25 -M -R "@RG\tID:group1\tSM:sample1\tPL:Illumina\tLB:lib1\tPU:unit1"
/bwa_build/Homo_sapiens/NCBI/build37.2/Sequence/BWAIndex/Homo_sapiens.GRCh38.dna.alt.fa
Cleaned_NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
Cleaned_NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
> NIST7035_TAAGGCGA_L001.sam

#Converting sam to bam
echo "Converting Sam to bam"
samtools view -Sb NIST7035_TAAGGCGA_L001.sam > NIST7035_TAAGGCGA_L001.bam

#Sorting bam file
echo "Sorting bam"
samtools sort NIST7035_TAAGGCGA_L001.bam -o NIST7035_TAAGGCGA_L001.sorted.bam -O bam -@ 25


#Command to remove duplicate reads
echo "MarkDuplicates is running"
./gatk MarkDuplicates
-I NIST7035_TAAGGCGA_L001.sorted.bam
-O dedup_NIST7035_TAAGGCGA_L001.bam
-M marked_dup.txt

echo "Genome directory file is generating"
java -jar picard.jar CreateSequenceDictionary
R=Homo_sapiens.GRCh38.dna.alt.fa
O=Homo_sapiens.GRCh38.dna.alt.dict

#indexing vcf files
./gatk IndexFeatureFile -I 1000GENOMES-phase_3.vcf


#Indexing reference genome using samtools
samtools faidx Homo_sapiens.GRCh38.dna.alt.fa

# reordering BAM files
echo "Bam file reordering"
java -jar picard.jar ReorderSam
I=dedup_NIST7035_TAAGGCGA_L001.bam
O=reordered_dedup_NIST7035_TAAGGCGA_L001.bam
R=Homo_sapiens.GRCh38.dna.alt.fa
CREATE_INDEX=TRUE


echo "Base recalibration is running"
./gatk BaseRecalibrator
--input reordered_dedup_NIST7035_TAAGGCGA_L001.bam
--reference genome.fa
--known-sites 1000GENOMES-phase_3.vcf
--output recal1_data.table

echo "Appling BQSR"
./gatk ApplyBQSR -R genome.fa
-I reordered_dedup_NIST7035_TAAGGCGA_L001.bam
--bqsr-recal-file recal1_data.table
-O reordered_dedup_NIST7035_TAAGGCGA_L001_AppBQSR.bam

echo "HaplotypeCaller is running"
./gatk HaplotypeCaller -R genome.fa
-I reordered_dedup_NIST7035_TAAGGCGA_L001_AppBQSR.bam
-O output.vcf.gz
-ERC GVCF
