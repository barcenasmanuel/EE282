#!/usr/bin/env bash
mamba activate ee282
# Navigate to the project directory
cd ~/myrepos/ee282
# Ensure the data directories exist
mkdir -p data/raw
mkdir -p data/processed
mkdir -p output
# Download file and checksum
gtf_url="https://ftp.flybase.net/releases/current/dmel_r6.60/gtf/dmel-all-r6.60.gtf.gz"
gtf_checksum_url="https://ftp.flybase.net/releases/current/dmel_r6.60/gtf/md5sum.txt"
wget -q -O data/raw/dmel-all-r6.60.gtf.gz $gtf_url
wget -q -O data/raw/md5sum_gtf.txt $gtf_checksum_url
# perform checksum 
cd data/raw
grep "dmel-all-r6.60.gtf.gz" md5sum_gtf.txt | md5sum -c -
cd ../..
# Process GTF file for feature counts
gunzip -c data/raw/dmel-all-r6.60.gtf.gz \
| bioawk -c gff '{print $3}' \
| sort | uniq -c | sort -nr > data/processed/feature_counts.txt

# Count genes per primary chromosome arm
declare -a arms=("2L" "2R" "3L" "3R" "4" "X" "Y")
for arm in "${arms[@]}"; do
  count=$(gunzip -c data/raw/dmel-all-r6.60.gtf.gz \
  | bioawk -c gff '$3 == "gene"' \
  | awk -v arm="$arm" '$1 == arm' \
  | wc -l)
  echo "$arm: $count"
done > data/processed/genes_per_primary_chromosome.txt

# Compile GTF report and save to output
{
  echo "GTF Annotation Summary:"
  echo "Feature Counts:"
  cat data/processed/feature_counts.txt
  echo ""
  echo "Total Number of Genes per Primary Chromosome Arm:"
  cat data/processed/genes_per_primary_chromosome.txt
} > output/gtf_output.txt