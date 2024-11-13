# Homework 3 - EE282 - Manuel Barcenas

## Genome Summary

### Description

In this task, I process chromosome data to generate a summary of genomic features such as total nucleotides, total Ns, and total sequences.

### Implementation

- **Script Location**: `code/scripts/hw3_genome_summary.sh`
- **Key Commands**:
  - Use `faSize` to calculate detailed genome statistics.
  - Use `grep` and `wc` to count specific nucleotides.

### Handling Ns Counting Challenge

One of the challenges I faced was accurately counting the total number of 'Ns' in the fasta file. Initially, I attempted to use `awk`, but it proved inefficient for this task due to the multiline nature of fasta sequences. `awk` processes input line-by-line, which can lead to inaccuracies when sequences are split across multiple lines.

To overcome this issue, I used `grep` to search for 'Ns'. `grep` pattern matches across files, making it ideal for counting occurrences of 'Ns' regardless of line breaks. This approach ensured that I obtained an accurate count of 'Ns' in the dataset.
### Code Snippet
```
{
# Example of a command used in the script
gunzip -c data/raw/dmel-all-chromosome-r6.60.fasta.gz | grep -o "N" | wc -l
}
```

### Results

- **Total Nucleotides**: 143726002
- **Total Ns**: 1154851
- **Total Sequences**: 1870

## Annotation Summary of GTF file:

### Description

This task involves processing GTF files to summarize annotation features and count genes per primary chromosome arm.

### Implementation

- **Script Location**: `code/scripts/hw3_annotation_summary.sh`
- **Key Commands**:
  - Use `bioawk` to parse GTF files.
  - Filter and count genes on primary chromosome arms.

### Handling Chromosome Arm Filtering Challenges

During the analysis of the GTF file, I encountered challenges in filtering out the correct chromosome arms. The GTF file included not only the primary chromosome arms (2L, 2R, 3L, 3R, 4, X, Y) but also additional scaffolds and contigs with numerical identifiers.

Our goal was to focus on the primary chromosome arms, excluding these additional entries. Initially, attempts to filter using basic string matching led to the inclusion of these unwanted portions. To address this, I implemented a precise filtering step using `awk` to match only the specified main chromosome arms.

### Code Snippet
```
{
gunzip -c data/raw/dmel-all-r6.60.gtf.gz | bioawk -c gff '\$3 == "gene"' | awk '\$1 ~ /^(2L|2R|3L|3R|4|X|Y)$/' | wc -l
}
```

### Results
- **Genes per Chromosome Arm**:
  - **2L**: 3508
  - **2R**: 3649
  - **3L**: 3481
  - **3R**: 4226
  - **4**: 114
  - **X**: 2704
  - **Y**: 113

## Running the Scripts

### Prerequisites

Before running the scripts, ensure that you have the necessary tools installed and your environment is set up correctly. This includes:

- **Conda Environment**: Make sure the `ee282` conda environment is activated.
- **Srun is used before runing scripts**

1. **Navigate to the Scripts Directory**:

```
cd ~/myrepos/ee282/code/scripts
./hw3_genome_summary.sh
./hw3_annotation_summary.sh
```
2. **Corresponding scripts**:


    - hw3_genome_summary.sh: This script processes the chromosome data and outputs a summary of genomic features.
    
    - hw3_annotation_summary.sh: This script processes the GTF file to summarize annotation features and counts genes per primary chromosome arm.
## Full codes below for convenience:

### All chromosome file
```
{
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
}
```

### GTF file

```
{
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
}
```
