import os
import subprocess
import shutil
import argparse
import glob
import logging
import sys

def setup_logging():
    logging.basicConfig(filename='assembly.log', filemode='a', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def setup_directories(directories):
    for directory in directories:
        os.makedirs(directory, exist_ok=True)

def run_subprocess(cmd, log_file):
    logging.info(f'Running command: {" ".join(cmd)}')
    with open(log_file, 'a') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        logging.error(f'Command failed: {" ".join(cmd)}')
        sys.exit(f'Error: Command failed. Check {log_file} for details.')

def run_assembly(trimmed_data_path, assembly_dir, spades_path='spades.py'):
    # Find all R1 paired files in the trimmed_data_path
    r1_paired_files = glob.glob(os.path.join(trimmed_data_path, '*_R1_paired.fastq*'))
    for r1_paired_file in r1_paired_files:
        base_name = os.path.basename(r1_paired_file)
        sample_name = base_name.replace('_R1_paired.fastq.gz', '').replace('_R1_paired.fastq', '')
        r2_paired_file = os.path.join(trimmed_data_path, f"{sample_name}_R2_paired.fastq.gz")
        if not os.path.exists(r2_paired_file):
            r2_paired_file = os.path.join(trimmed_data_path, f"{sample_name}_R2_paired.fastq")
            if not os.path.exists(r2_paired_file):
                logging.warning(f"No matching R2 paired file for {r1_paired_file}")
                continue
        # Look for unpaired files
        unpaired_files = []
        r1_unpaired_file = os.path.join(trimmed_data_path, f"{sample_name}_R1_unpaired.fastq.gz")
        if not os.path.exists(r1_unpaired_file):
            r1_unpaired_file = os.path.join(trimmed_data_path, f"{sample_name}_R1_unpaired.fastq")
        if os.path.exists(r1_unpaired_file):
            unpaired_files.append(r1_unpaired_file)
        r2_unpaired_file = os.path.join(trimmed_data_path, f"{sample_name}_R2_unpaired.fastq.gz")
        if not os.path.exists(r2_unpaired_file):
            r2_unpaired_file = os.path.join(trimmed_data_path, f"{sample_name}_R2_unpaired.fastq")
        if os.path.exists(r2_unpaired_file):
            unpaired_files.append(r2_unpaired_file)
        sample_output_dir = os.path.join(assembly_dir, "contigs", "samples", sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        cmd = [
            spades_path,
            '--pe1-1', r1_paired_file,
            '--pe1-2', r2_paired_file,
            '-o', sample_output_dir
        ]
        # Include unpaired reads if available
        if unpaired_files:
            cmd.extend(['--s1'] + unpaired_files)
        run_subprocess(cmd, os.path.join(sample_output_dir, 'spades_output.log'))
        # Copy the resulting contigs to the all_contigs_dir
        contigs_src = os.path.join(sample_output_dir, "contigs.fasta")
        contigs_dest = os.path.join(assembly_dir, "contigs", "all_contigs", f"{sample_name}.fasta")
        if os.path.exists(contigs_src):
            shutil.copy(contigs_src, contigs_dest)
        else:
            logging.warning(f"Contigs file not found for sample {sample_name}")

def reformat_contigs(all_contigs_dir, filtered_contigs_dir, minlength=1000, reformat_path='reformat.sh'):
    contig_files = glob.glob(os.path.join(all_contigs_dir, '*.fasta'))
    for file in contig_files:
        sample_name = os.path.basename(file).split('.')[0]
        output_file = os.path.join(filtered_contigs_dir, f"{sample_name}.fasta")
        cmd = [
            reformat_path,
            f"in={file}",
            f"out={output_file}",
            f"minlength={minlength}"
        ]
        run_subprocess(cmd, 'reformat_output.log')

def run_quast(contigs_dir, output_dir, quast_path='quast'):
    contig_files = glob.glob(os.path.join(contigs_dir, '*.fasta'))
    if not contig_files:
        logging.warning(f"No contig files found in {contigs_dir}")
        return
    cmd = [quast_path] + contig_files + ['-o', output_dir]
    run_subprocess(cmd, os.path.join(output_dir, 'quast_output.log'))

def run_bbmap(trimmed_data_path, filtered_contigs_dir, coverage_dir, bbmap_path='bbmap.sh'):
    # Find all R1 paired files in the trimmed_data_path
    r1_paired_files = glob.glob(os.path.join(trimmed_data_path, '*_R1_paired.fastq*'))
    for r1_paired_file in r1_paired_files:
        base_name = os.path.basename(r1_paired_file)
        sample_name = base_name.replace('_R1_paired.fastq.gz', '').replace('_R1_paired.fastq', '')
        r2_paired_file = os.path.join(trimmed_data_path, f"{sample_name}_R2_paired.fastq.gz")
        if not os.path.exists(r2_paired_file):
            r2_paired_file = os.path.join(trimmed_data_path, f"{sample_name}_R2_paired.fastq")
            if not os.path.exists(r2_paired_file):
                logging.warning(f"No matching R2 paired file for {r1_paired_file}")
                continue
        ref_file = os.path.join(filtered_contigs_dir, f"{sample_name}.fasta")
        if not os.path.exists(ref_file):
            logging.warning(f"Reference contigs not found for sample {sample_name}")
            continue
        sample_coverage_dir = os.path.join(coverage_dir, sample_name)
        os.makedirs(sample_coverage_dir, exist_ok=True)
        cmd = [
            bbmap_path,
            f"in1={r1_paired_file}",
            f"in2={r2_paired_file}",
            f"ref={ref_file}",
            f"covstats={os.path.join(sample_coverage_dir, sample_name + '_covstats.txt')}",
            f"covhist={os.path.join(sample_coverage_dir, sample_name + '_covhist.tsv')}",
            f"basecov={os.path.join(sample_coverage_dir, sample_name + '_basecov.txt')}",
            "usejni=t",
            f"out={os.path.join(sample_coverage_dir, sample_name + '_mapped.bam')}"
        ]
        run_subprocess(cmd, os.path.join(sample_coverage_dir, 'bbmap_output.log'))

def main(args=None):
    parser = argparse.ArgumentParser(description='Genome Assembly Pipeline')
    parser.add_argument('--output_dir', required=True, help='Path to the output directory')
    parser.add_argument('--trimmed_data_path', required=True, help='Path to the trimmed FASTQ files')
    parser.add_argument('--spades_path', default='spades.py', help='Path to SPAdes executable')
    parser.add_argument('--quast_path', default='quast', help='Path to QUAST executable')
    parser.add_argument('--bbmap_path', default='bbmap.sh', help='Path to BBMap executable')
    parser.add_argument('--minlength', default=1000, type=int, help='Minimum length of contigs to keep')

    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Set the default reformat_path to the reformat.sh in the src directory
    reformat_default_path = os.path.join(script_dir, 'reformat.sh')
    parser.add_argument('--reformat_path', default=reformat_default_path, help='Path to reformat.sh script')

    args = parser.parse_args(args)

    setup_logging()

    base_dir = args.output_dir
    assembly_dir = os.path.join(base_dir, "Assembly")
    trimmed_data_path = args.trimmed_data_path
    all_contigs_dir = os.path.join(assembly_dir, "contigs", "all_contigs")
    filtered_contigs_dir = os.path.join(assembly_dir, "contigs", "filtered_contigs")
    quast_dir = os.path.join(assembly_dir, "Quast")
    pre_quast_dir = os.path.join(quast_dir, "pre_Quast")
    filtered_quast_dir = os.path.join(quast_dir, "filtered_Quast")
    coverage_dir = os.path.join(assembly_dir, "Coverage")

    # Setup directories
    setup_directories([
        os.path.join(assembly_dir, "contigs", "samples"),
        all_contigs_dir,
        filtered_contigs_dir,
        pre_quast_dir,
        filtered_quast_dir,
        coverage_dir
    ])

    # Run assembly
    run_assembly(trimmed_data_path, assembly_dir, spades_path=args.spades_path)

    # Reformat contigs
    reformat_contigs(all_contigs_dir, filtered_contigs_dir, minlength=args.minlength, reformat_path=args.reformat_path)

    # Run QUAST on all contigs
    run_quast(all_contigs_dir, pre_quast_dir, quast_path=args.quast_path)

    # Run QUAST on filtered contigs
    run_quast(filtered_contigs_dir, filtered_quast_dir, quast_path=args.quast_path)

    # Run BBMap to calculate coverage
    run_bbmap(trimmed_data_path, filtered_contigs_dir, coverage_dir, bbmap_path=args.bbmap_path)

if __name__ == "__main__":
    main()
