#!/bin/bash
## job script header that requests 8 threads
#
#SBATCH --partition=normal
#SBATCH --ntasks=8
#SBATCH --mem=1024
#SBATCH --output=sort_%J_stdout.txt
#SBATCH --error=sort_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=sort
#SBATCH --mail-user=lakshmibhavaniuppaluri-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/bhavani/RNAseq-analysis
#
# script to perform samtools sorting of trimmed and aligned paired end reads
# usage: sbatch samtools_trainingPrompts.sh inputsFile outputsFile
# usage Ex: sbatch samtools_trainingPrompts.sh /scratch/bhavani/project1/adipocyte/aligned /scratch/bhavani/project1/adipocyte
# usage Ex: sbatch samtools_trainingPrompts.sh /YOUR/PATH/jurkat/aligned /YOUR/PATH/jurkat

# required software for OSCER
module load SAMtools

# retrieve paired reads absolute path for alignment
inputsPath="$1"

# retrieve analysis outputs absolute path
outputsPath="$2"

# make a new directory for project analysis files
anOut=$outputsPath"/sorted"

# make the new outputs directory
mkdir $anOut

# move to the outputs directory
cd $anOut

# loop through all reads and sort sam/bam files for input to samtools
for f1 in "$inputsPath"/*/; do

	# name of aligned file
	curAlignedSample=$f1"accepted_hits.bam"

	# trim file path from current folder name
	curSampleNoPath=$(basename $f1)

	# create directory for current sample outputs
	mkdir $curSampleNoPath

	# print status message
	echo "Sample $curSampleNoPath is being name sorted..."

	# run samtools to prepare mapped reads for sorting using 8 threads
	samtools sort -@ 8 -n -o $curSampleNoPath"/sortedName.bam" -T "/tmp/"$curSampleNoPath".sortedName.bam" $curAlignedSample

	# run samtools fixmate to update paired-end flags for singletons
	samtools fixmate -m $curSampleNoPath"/sortedName.bam" $curSampleNoPath"/accepted_hits.bam"

	# remove the current sortedName.bam file
	rm "$curSampleNoPath"/sortedName.bam

	# remove duplicate reads using samtools markdup
	samtools markdup -r $curSampleNoPath"/accepted_hits.bam" $curSampleNoPath"/noDups.bam"

	# print status message
	echo "Sample $curSampleNoPath has been name sorted!"	

# end loop
done

# clean up and remove extra files
rm -r $inputsPath

# print final status message
echo "Analysis complete!"


