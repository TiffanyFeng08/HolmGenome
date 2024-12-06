#!/bin/bash
#SBATCH --account=def-jronho
#SBATCH --time=5:00:00
#SBATCH --job-name=PlamidSpades
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=zhangbin.cai@mail.mcgill.ca
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

module load StdEnv/2023
module load spades/3.15.4

# Define input and output directories
input_dir="/lustre04/scratch/zhangbin/Daryna/raw"
output_dir="/lustre04/scratch/zhangbin/Daryna/plasmidspades_output"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Load SPAdes (if not already in the PATH)
# module load spades  (Uncomment if using a module system)

# Loop over all paired-end read files (R1 and R2)
for r1 in "$input_dir"/*_R1_001.fastq.gz; do
    # Find the matching R2 file
    r2="${r1/_R1_001/_R2_001}"

    # Extract sample name (assuming filenames like 10000467_S39_L001_R1_001.fastq.gz)
    sample_name=$(basename "$r1" | cut -d'_' -f1)

    # Create a directory for the current sample in the output directory
    sample_output_dir="$output_dir/$sample_name"
    mkdir -p "$sample_output_dir"

    # Run PlasmidSPAdes
    spades.py --plasmid \
          -1 "$r1" \
          -2 "$r2" \
          -o "$sample_output_dir" \
          -t 16 \
          -m 128

    # Optionally: Check exit status and log it
    if [ $? -eq 0 ]; then
        echo "PlasmidSPAdes successfully completed for $sample_name"
    else
        echo "PlasmidSPAdes failed for $sample_name" >&2
    fi
done
