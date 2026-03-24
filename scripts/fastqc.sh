#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=FastQC_Metagenomics

# 1. Load the FastQC module
module load fastqc/0.12.1 multiqc/1.14

# 2. Create a folder for the fastqc reports
mkdir -p fastqc_reports

# 3. Run FastQC on all .fastq files in the output folder
# -threads to 8
fastqc ./fastq_files/*.fastq -o ./fastqc_reports/ -t 8

echo "FastQC complete. Reports are in ./fastqc_reports/"

# 4. Load the MultiQC module
module load multiqc/1.14

# 5. Create a repo for the aggregated reports
mkdir -p multiqc_final

# 6. Aggregate all fastqc reports into one
conda activate genome_assembly
multiqc ./fastqc_reports/ -o ./multiqc_final/

echo "MultiQC complete."