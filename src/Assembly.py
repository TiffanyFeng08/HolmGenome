import os
import subprocess
import shutil

# Define directories
base_dir = "path/to/base"
assembly_dir = os.path.join(base_dir, "Assembly")
trimmed_data_path = os.path.join(base_dir, 'Trim_data')
all_contigs_dir = os.path.join(assembly_dir, "contigs", "all_contigs")
filtered_contigs_dir = os.path.join(assembly_dir, "contigs", "filtered_contigs")
quast_dir = os.path.join(assembly_dir, "Quast")
pre_quast_dir = os.path.join(quast_dir, "pre_Quast")
filtered_quast_dir = os.path.join(quast_dir, "filtered_Quast")
multi_quast_dir = os.path.join(quast_dir, "multiQuast")
coverage_dir = os.path.join(assembly_dir, "Coverage")
multi_coverage_dir = os.path.join(coverage_dir, "multiCoverage")

# Ensure necessary directories exist
os.makedirs(all_contigs_dir, exist_ok=True)
os.makedirs(filtered_contigs_dir, exist_ok=True)
os.makedirs(pre_quast_dir, exist_ok=True)
os.makedirs(filtered_quast_dir, exist_ok=True)
os.makedirs(multi_quast_dir, exist_ok=True)
os.makedirs(coverage_dir, exist_ok=True)
os.makedirs(multi_coverage_dir, exist_ok=True)

def run_assembly(input_dir, output_dir):
    """
    Run the assembly process on the input data and store the results in the output directory.
    """
    for file in os.listdir(input_dir):
        if file.endswith("_R1_001_paired.fastq.gz"):
            sample_name = file.split("_")[0]
            # Example assembly command, replace with actual assembly tool command
            subprocess.run([
                "spades.py",
                "--pe1-1", os.path.join(input_dir, f"{sample_name}_R1_001_paired.fastq.gz"),
                "--pe1-2", os.path.join(input_dir, f"{sample_name}_R2_001_paired.fastq.gz"),
                "-o", os.path.join(output_dir, sample_name)
            ])
            # Move the resulting contigs to the all_contigs_dir
            shutil.move(os.path.join(output_dir, sample_name, "contigs.fasta"), os.path.join(all_contigs_dir, f"{sample_name}.fasta"))

# Run the assembly process
run_assembly(trimmed_data_path, all_contigs_dir)

# Process each file in all_contigs_dir
for file in os.listdir(all_contigs_dir):
    if file.endswith("contigs.fasta"):
        sample_name = file.split(".")[0]
        # Reformat contigs
        subprocess.run(["reformat.sh", "in=" + os.path.join(all_contigs_dir, file), "out=" + os.path.join(filtered_contigs_dir, f"{sample_name}_filtered_contigs.fasta"), "minlength=1000"])

# Use QUAST to test quality of all contigs
subprocess.run(["quast", os.path.join(all_contigs_dir, "*.fasta"), "-o", pre_quast_dir])

# Use QUAST to test quality of filtered contigs
subprocess.run(["quast", os.path.join(filtered_contigs_dir, "*.fasta"), "-o", filtered_quast_dir])

# Generate multiQC report for QUAST results
subprocess.run(["multiqc", quast_dir, "-o", multi_quast_dir])

# Run bbmap to calculate coverage
for file in os.listdir(trimmed_data_path):
    if file.endswith("_R1_001_paired.fastq.gz"):
        sample_name = file[:-len("_R1_001_paired.fastq.gz")]
        sample_coverage_dir = os.path.join(coverage_dir, sample_name)
        os.makedirs(sample_coverage_dir, exist_ok=True)
        subprocess.run([
            "bbmap.sh",
            f"in1={os.path.join(trimmed_data_path, sample_name + '_R1_001_paired.fastq.gz')}",
            f"in2={os.path.join(trimmed_data_path, sample_name + '_R2_001_paired.fastq.gz')}",
            f"ref={os.path.join(filtered_contigs_dir, sample_name + '_filtered_contigs.fasta')}",
            f"covstats={os.path.join(sample_coverage_dir, sample_name + '_covstats.txt')}",
            f"covhist={os.path.join(sample_coverage_dir, sample_name + '_covhist.tsv')}",
            f"basecov={os.path.join(sample_coverage_dir, sample_name + '_basecov.txt')}",
            "usejni=t",
            f"out={os.path.join(sample_coverage_dir, sample_name + '_mapped.bam')}"
        ])
# Generate multiQC report for BBMap coverage results
subprocess.run(["multiqc", coverage_dir, "-o", multi_coverage_dir])