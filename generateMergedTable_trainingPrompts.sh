#!/bin/bash
## job script header that requests 1 thread
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=1024
#SBATCH --output=merge_%J_stdout.txt
#SBATCH --error=merge_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=merge
#SBATCH --mail-user=lakshmibhavaniuppaluri-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/bhavani/RNAseq-analysis
#
# script to generate guide file and merge gene counts using the merge_tables.py script
# usage: bash generateMergedTable_trainingPrompts.sh inputsFile outputsFile
# usage ex: bash generateMergedTable_trainingPrompts.sh /scratch/bhavani/project1/Adipocyte/counted /scratch/bhavani/project1/Adipocyte
# usage ex: bash generateMergedTable_trainingPrompts.sh /YOUR/PATH/Jurkat/counted /YOUR/PATH/Jurkat

# retrieve paired reads absolute path for alignment
inputsPath="$1"

# retrieve analysis outputs absolute path
outputsPath="$2"

# make a new directory for project analysis files
anOut=$outputsPath"/merged"

# make the new outputs directory
mkdir $anOut

# move to the outputs directory
cd $anOut

# name tmp guide file for merging
tmpGuide=$inputsPath"/tmp_guideFile.txt"

# loop through all counted paired reads
for f1 in "$inputsPath"/*/; do

	# trim extension from current file
	currTag=$(echo $f1 | sed 's/.$//')

	# create a guide file
	echo $f1"counts.txt $currTag" >> $tmpGuide

# end loop
done

# move to location of merge_tagles.py script in the util director
cd ../util

# merge gene counts using the generated guide file and merge_tables.py
python merge_tables.py $tmpGuide

# trim / from the inputs path
path=$(echo $inputsPath | sed "s/\///g")

# remove file paths from sample tags
cat merged_counts.txt | sed 's/\///g' | sed "s/$path//g" > $inputsPath"/"$projectDir"_merged_counts.txt"

# clean up sample names for downstream analysis with R
header=$(head -1 $inputsPath"/"$projectDir"_merged_counts.txt" | tr '\t' '\n' | sed "s/\-/\_/g" | sed "s/Jurkat.*\_S//g" | sed "s/3T3.*\_S//g" | sed "s/Undetermined\_S//g" | awk '{$1 = sprintf("S%02d", $1); print}' | tr '\n' '\t' | cut -d$'\t' -f2-)

# add the header to the final re-formatted merged counts file
echo -e "gene\t"$header | sed 's/\t*$//' > $inputsPath"/"$projectDir"_merged_counts_formatted.txt"

# add the count data to the final re-formatted merged counts file
tail -n +2 $inputsPath"/"$projectDir"_merged_counts.txt" >> $inputsPath"/"$projectDir"_merged_counts_formatted.txt"

# clean up and remove the tmp guide file
rm $tmpGuide

# clean up and remove the un-formatted merged counts file
#rm merged_counts.txt
