#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G              
#SBATCH --time=03:00:00
#SBATCH --job-name=SRA_to_FASTQ

#1) Load the sratoolkit
module load sra-toolkit

mkdir -p fastq_files

#Loop through the folders prefetched
while read -r id; do
  #Check if the folder or file exists before splitting
  if [ -d "$id" ] || [ -f "$id.sra" ]; then
    echo "Splitting $id..."
    #Fasterqdump has to go through the folder
    fasterq-dump "./$id" --split-3 --threads 8 --outdir ./fastq_files --progress
  else
    echo "Skipping $id - not downloaded yet."
  fi
done < sra_list.txt
