#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=8
#SBATCH --mem=1024
#SBATCH --output=align_%J_stdout.txt
#SBATCH --error=align_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=align
#SBATCH --mail-user=lakshmibhavaniuppaluri-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/bhavani/RNAseq-analysis
#

# script to perform trimmomatic trimming of paired end reads
# usage: sbatch hisat2_projects.sh readPath outputsPath
# usage Ex: sbatch hisat2_projects.sh /scratch/bhavani/GBCF_Data/Adipocyte/220705_Yoon_Adipocyte_Pool2_RNAseq /scratch/bhavani/project1/adipocyte /scratch/bhavani/project1/adipocyte/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# required software for OSCER
module load HISAT2
module load SAMtools

# retrieve paired reads absolute path for alignment
inputsPath="$1"

# retrieve analysis outputs absolute path
outputsPath="$2"

#retrieve genome
ref="$3"

# make a new directory for project analysis files
mkdir $outputsPath

# name of a new directory for outputs of this analysis stage
anOut=$outputsPath"/alined"

# make the new outputs directory
mkdir $anOut

# move to the outputs directory
cd $anOut

#Create build output directory for Hisat reference
buildOut="build"
mkdir $buildOut
#Trim path and file extension from build file
refNoPath=$(basename $ref)
refNoEx=$(echo $refNoPath | sed 's/\.fasta//' | sed 's/\.fa//')
#Copy genome build fasta file to hisat2 build folder
cp "$ref" "$buildOut"/"$refNoEx"
#Trim file extension
refNoPath=$(echo $refNoEx | sed 's/\.fa//g')
#Begin hisat2 build
echo "Beginning hisat2 build... "
hisat2-build -p 8 -f "$buildOut"/"$refNoEx" "$buildOut"/"$refNoPath"
echo "hisat2 build complete!"

#Loop through all forward and reverse paired reads and run Hisat2 on each pair
# using 8 threads and samtools to convert output sam files to bam
for f1 in "$inputsPath"/*pForward.fq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/.pForward\.fq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/.pForward\.fq\.gz//')
	#Create directory for current sample outputs
	mkdir "$curSampleNoPath"
	#Print status message
	echo "Processing $curSampleNoPath"
	#Run hisat2 with default settings
	hisat2 -p 8 -q -x $buildOut"/"$refNoEx -1 $f1 -2 $curSample"_pReverse.fq.gz" -S $curSampleNoPath"/accepted_hits.sam" --summary-file $curSampleNoPath"/alignedSummary.txt"
	#Convert output sam files to bam format for downstream analysis
	samtools view -@ 8 -bS $curSampleNoPath"/accepted_hits.sam" > $curSampleNoPath"/accepted_hits.bam"
	#Remove the now converted .sam file
	rm $curSampleNoPath"/accepted_hits.sam"
	#Print status message
	echo "Processed!"
done

#Clean up
rm -r "build"
rm -r $inputsPath

#Print status message
echo "Analysis complete!"
