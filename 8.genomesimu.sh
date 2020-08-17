
###################################################################################
#
#This program will able to simulate the genome inserting n number od SNP randomely.
#The name of the package is simuG used for this simulation operaion.
#
###################################################################################


pl=/data/sata_data/workshop/wsu28/mosquito/simuG


rg=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/ncbi-genomes-2020-06-18



perl $pl/simuG.pl \
     -refseq /data/sata_data/workshop/wsu28/mosquito/mosquitoRG/GSE113256_AGWG.HiC.fasta.gz \
     -snp_count 1000 \
     -prefix output_prefix



perl $pl/simuG.pl \
     -refseq $rg/finalRefGen.fa \
     -snp_count 1000000 \
     -prefix output_1000000snps

