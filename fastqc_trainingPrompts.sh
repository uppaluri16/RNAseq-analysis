#!/bin/bash
## job script header
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=1024
#SBATCH --output=fastqc_%J_stdout.txt
#SBATCH --error=fastqc_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=fastqc
#SBATCH --mail-user=lakshmibhavaniuppaluri-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/bhavani/RNAseq-analysis
#

# script to perform fastqc quality control of paired end reads
# usage: sbatch fastqc_trainingPrompts.sh readPath outputsPath
# usage Ex: sbatch fastqc_trainingPrompts.sh /scratch/bhavani/GBCF_Data/Adipocyte/220705_Yoon_Adipocyte_Pool2_RNAseq /scratch/bhavani/project1


# required software for OSCER servers
module load FastQC

# retrieve paired reads absolute path for alignment
readPath="$1"

# retrieve analysis outputs absolute path
outputsPath="$2"

# make the new directory for project analysis files
mkdir $outputsPath

# name of a new directory for outputs of this analysis stage
qcOut=$outputsPath"/qc"

# make the outputs directory
mkdir $qcOut

# move to the outputs directory
cd $qcOut

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do

	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')

	# set paired file name
	f2=$curSample"_R2_001.fastq.gz"

	# print status message
	echo "Processing $curSample"

	# perform QC on the first paired end read for the current sample
	fastqc $f1 -o $qcOut --extract

	# perform QC on the second paired end read for the current sample
	fastqc $f2 -o $qcOut --extract

# end loop
done

# print final status message
echo "Analysis Complete"
