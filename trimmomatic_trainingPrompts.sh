#!/bin/bash
## job script header that requests 8 threads
#
#SBATCH --partition=normal
#SBATCH --ntasks=8
#SBATCH --mem=1024
#SBATCH --output=trim_%J_stdout.txt
#SBATCH --error=trim_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=trim
#SBATCH --mail-user=lakshmibhavaniuppaluri-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/bhavani/RNAseq-analysis
#

# script to perform trimmomatic trimming of paired end reads
# usage: sbatch trimmomatic_trainingPrompts.sh readPath outputsPath
# usage Ex: sbatch trimmomatic_trainingPrompts.sh /scratch/bhavani/GBCF_Data/Adipocyte/220705_Yoon_Adipocyte_Pool2_RNAseq /scratch/bhavani/project1/adipocyte
# usage Ex: sbatch trimmomatic_trainingPrompts.sh /afs/crc.nd.edu/group/genomics/R2D2/220707_Yoon_Jurkat_Pool1_RNAseq /scratch365/ebrooks5/yoon_July2022/220707_Yoon_Jurkat_Pool1_RNAseq

# required software for OSCER
module load Trimmomatic

# retrieve paired reads absolute path for alignment
readPath="$1"

# retrieve analysis outputs absolute path
outputsPath="$2"

# retrieve adapter absolute path for alignment
adapterPath="/opt/oscer/software/Trimmomatic/0.39-Java-1.8.0_191/adapters/TruSeq3-PE.fa:2:30:10"

# make a new directory for project analysis files
mkdir $outputsPath

# name of a new directory for outputs of this analysis stage
trimOut=$outputsPath"/trimmed"

# make the new outputs directory
mkdir $trimOut

# move to the outputs directory
cd $trimOut

# loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*_R1_001.fastq.gz; do

	# trim extension from current file name
	curSample=$(echo $f1 | sed 's/_R._001\.fastq\.gz//')

	# set paired file name
	f2=$curSample"_R2_001.fastq.gz"

	# print status message
	echo "Processing $curSample"

	# trim to sample tag
	sampleTag=$(basename $f1 | sed 's/_R._001\.fastq\.gz//')

	# determine phred score for trimming
	if grep -iF "Illumina 1.5" $outputsPath"/qc/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=64
	elif grep -iF "Illumina 1.9" $outputsPath"/qc/"$sampleTag"_R1_001_fastqc/fastqc_data.txt"; then
		score=33
	else
		echo "ERROR: Illumina encoding not found... exiting"
		exit 1
	fi
	
	# perform adapter trimming on paired reads using 8 threads
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 -phred"$score" $f1 $f2 $sampleTag"_pForward.fq.gz" $sampleTag"_uForward.fq.gz" $sampleTag"_pReverse.fq.gz" $sampleTag"_uReverse.fq.gz" ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# end loop
done

# print final status message
echo "Analysis Complete"
