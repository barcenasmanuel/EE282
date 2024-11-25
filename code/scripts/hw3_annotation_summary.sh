#!/usr/bin/env bash
mamba activate ee282
# Navigate to the project directory
cd ~/myrepos/ee282
# Ensure the data directories exist
mkdir -p data/raw
mkdir -p data/processed
mkdir -p output
# Download file and checksum
chromosome_url="https://ftp.flybase.net/releases/current/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz"
chromosome_checksum_url="https://ftp.flybase.net/releases/current/dmel_r6.60/fasta/md5sum.txt"
wget -q -O data/raw/dmel-all-chromosome-r6.60.fasta.gz $chromosome_url
wget -q -O data/raw/md5sum_chromosome.txt $chromosome_checksum_url
# Perform the checksum 
cd data/raw
grep "dmel-all-chromosome-r6.60.fasta.gz" md5sum_chromosome.txt | md5sum -c -
cd ../..
# Process chromosome file
gunzip -c data/raw/dmel-all-chromosome-r6.60.fasta.gz \
| faSize -detailed /dev/stdin \
| tee data/processed/chromosome_summary.txt \
| cut -f 2 > data/processed/chromosome_sizes.txt

# Extract chromosome statistics
total_nucleotides=$(awk '{sum += $2} END {print sum}' data/processed/chromosome_summary.txt)
total_ns=$(gunzip -c data/raw/dmel-all-chromosome-r6.60.fasta.gz | grep -o "N" | wc -l)
total_sequences=$(wc -l < data/processed/chromosome_summary.txt)

# Compile chromosome report and save to output
{
  echo "Chromosome Summary:"
  echo "Total nucleotides: $total_nucleotides"
  echo "Total Ns: $total_ns"
  echo "Total sequences: $total_sequences"
} > output/chromosome_output.txt