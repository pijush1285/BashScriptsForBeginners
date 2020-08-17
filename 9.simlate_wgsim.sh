#!/bin/bash


wgsim=/data/sata_data/workshop/wsu28/packages/wgsim/wgsim
wgs=/data/sata_data/workshop/wsu28/packages/wgsim


#For compilation
gcc -g -O2 -Wall -o $wgsim $wgs/wgsim.c -lz -lm

#Address of the reference genome
#refg=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/ncbi-genomes-2020-06-18/GCF_genomic.fa
#refg=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simu_genome1/output_prefix.simseq.genome.fa
refg=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simu_genome204/output_1000000snps.simseq.genome.fa


#Location of the output
#result=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simu_genome1
result=/data/sata_data/workshop/wsu28/mosquito/mosquitoRG/simulated_data/simu_genome204


#Change the location of the output
cd $result


#Doing the simulation
#$wgsim -1151 -2151 -d500 -r0 -e0 -N10000 -R0 -X0 $refg GCF_GPI.read1.fq GCF_GPI.read2.fq
#$wgsim -1151 -2151 -d500 -r0 -e0 -N1000000 -R0 -X0 $refg newRead1.fq newRead2.fq

#Generating 4 files containing 138 m reads in each file.
#$wgsim -1151 -2151 -d500 -r0 -e0 -N138000000 -R0 -X0 $refg NRead1.fq NRead2.fq
$wgsim -1151 -2151 -d500 -r0 -e0 -N138000000 -R0 -X0 $refg NRead3.fq NRead4.fq


