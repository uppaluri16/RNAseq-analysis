#!/bin/bash
module load Trimmomatic
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 -phred33 /scratch/bhavani/GBCF_Data/Adipocyte/220705_Yoon_Adipocyte_Pool2_RNAseq/Undetermined_S0_R1_001.fastq.gz /scratch/bhavani/GBCF_Data/Adipocyte/220705_Yoon_Adipocyte_Pool2_RNAseq/Undetermined_S0_R2_001.fastq.gz Undetermined_S0_pForward.fq.gz Undetermined_S0_uForward.fq.gz Undetermined_S0_pReverse.fq.gz Undetermined_S0_uReverse.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
