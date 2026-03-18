#!/bin/bash
#SBATCH --account=def-itobias
#SBATCH --cpus-per-task=32      # Kraken needs more CPUs
#SBATCH --mem=64G               # Needed to load 50GB database
#SBATCH --time=12:00:00
#SBATCH --job-name=Kraken2_taxonomy

module load kraken2/2.1.6 bracken/3.0

#Using standard database path on Narval
DB="/scratch/ssaab/kraken_db"

#Make directory for Kraken output
mkdir -p kraken_results

while read -r id; do
  echo "Classifying $id..."

  #1. Run Kraken2
  #--threads = 32
  #--minimum-base-quality = 20
  #--db = $DB
  #--paired because paired end reads
  #--confidence = 0.15 to avoid false positives
  #--report kraken_results/${id}.report \
  #--output - \
  #--use-names to use scientific names instead of just taxon ID
  kraken2 --db $DB \
          --threads 32 \
          --paired \
          --confidence 0.15 \
          --minimum-base-quality 20 \
          --use-names \
          --report kraken_results/${id}.report \
          --output - \
          fastq_files/${id}_1.fastq fastq_files/${id}_2.fastq

  #2. Bracken run (v3.0)
  #-d = Database repo
  #-i = The .report file from Kraken2
  #-o = Bracken output for R
  #-w = Bracken style report (not needed)
  #-r = 150 (multiqc showed 150bp reads)
  #-l = S (species level)
  bracken -d $DB \
          -i kraken_results/${id}.report \
          -o kraken_results/${id}.bracken \
          -r 150 \
          -l S

done < sra_list.txt
