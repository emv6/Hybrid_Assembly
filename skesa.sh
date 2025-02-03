#!/bin/bash -e
#SBATCH --cpus-per-task=20 --mem 50Gb --time 166:00:00 -J EV_SKESA_MULTI

module purge
module load SKESA/2.4.0-gimpi-2022a_saute.1.3.0_1

INPUT=Illumina.txt #txt file filename in column 1, R1 in column 2 and R2 in column 3

#Output Directory for contigs
mkdir -p ~/SKESA_Contigs/
output_dir=~/SKESA_Contigs

#Ensure the INPUT variable is set
if [ -z $INPUT ]; then echo "Error: INPUT variable is not set!"
exit 1
fi

#Check the INPUT file exists
if [ ! -f $INPUT ]; then echo "Error: File $INPUT not found!"
exit 1
fi

#Read the input file line by line
while IFS=$'\t' read -r isolate_name r1_path r2_path;
do
if [ -f $r1_path ] && [ -f $r2_path ]; then

echo "Starting processing for isolate: $isolate_name"

#Run Skesa
echo Running: skesa --reads \"$r1_path\" \"$r2_path\" --contigs_out \"$output_dir/${isolate_name}_contigs.fasta\ ""

skesa --reads $r1_path $r2_path --contigs_out $output_dir/${isolate_name}_contigs.fasta

if [ $? -eq 0 ]; then
echo "Finished processing for isolate: $isolate_name"
else
"Error: SKESA failed for isolate: $isolate_name"
fi
else
echo "Warning: Files for isolate $isolate_name not found, skipping....."
fi
done < $INPUT

echo "Contigs generation completed. Check the $output_dir directory"
