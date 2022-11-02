#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=1024
#SBATCH --output=mousedata_%J_stdout.txt
#SBATCH --error=mousedata_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=mousedata
#SBATCH --mail-user=lakshmibhavaniuppaluri-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/bhavani/RNAseq-analysis/
#

# BASH script to retrieve genomic data for project 220705
# usage: bash genomicData_220705.sh outputsPath
# usage ex: bash genomicData_220705.sh /scratch365/ebrooks5/yoon_July2022/220705_Yoon_Adipocyte_Pool2_RNAseq

# retrieve outputs absolute path
outputsPath="$1"

# move to the outputs directory
cd $outputsPath

# retrieve the genome assembly
wget https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# retrieve the gene annotations
wget https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz


