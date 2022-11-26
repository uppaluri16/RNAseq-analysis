#!/bin/bash
## job script header that requests 1 thread
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=1024
#SBATCH --output=count_%J_stdout.txt
#SBATCH --error=count_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=count
#SBATCH --mail-user=lakshmibhavaniuppaluri-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/bhavani/RNAseq-analysis
#
# script to run htseq-count on trimmed, aligned, and name sorted paired end reads
# usage: sbatch htseq_trainingPrompts.sh inputsFile outputsFile genomeFeatures
# usage Ex: sbatch htseq_trainingPrompts.sh /scratch/bhavani/Adipocyte/sorted /scratch/bhavani/Adipocyte /scratch/bhavani/Adipocyte/Mus_musculus.GRCm39.108.gtf
# usage Ex: sbatch htseq_trainingPrompts.sh /YOUR/PATH/Jurkat/sorted /YOUR/PATH/Jurkat /YOUR/PATH/Jurkat/....gtf

# required software for OSCER
module load HTSeq/0.6.1p1-intel-2016a-Python-2.7.10
pip install pysam

# retrieve paired reads absolute path for alignment
inputsPath="$1"

# retrieve analysis outputs absolute path
outputsPath="$2"

# retrieve genome features absolute path
genomeFile="$3"

# make a new directory for project analysis files
anOut=$outputsPath"/counted"

# make the new outputs directory
mkdir $anOut

# move to the outputs directory
cd $anOut

# loop through all sorted forward and reverse paired reads and store the file locations in an array
for f1 in "$inputsPath"/*/; do

	# name of aligned file
	curAlignedSample=$f1"accepted_hits.bam"

	# trim file path from current folder name
	curSampleNoPath=$(basename $f1)

	# create directory for current sample outputs
	mkdir $curSampleNoPath
	
	# output status message
	echo "Sample $curSampleNoPath is being name counted..."

	# count reads using htseq-count
	python -m htseq-count -f bam -s no -t gene $curAlignedSample $genomeFile > $curSampleNoPath"/counts.txt"

	# output status message
	echo "Sample $curSampleNoPath has been name counted!"

# end loop
done

#clean up
#rm -r $inputsPath

# print status message
echo "Analysis complete!"
