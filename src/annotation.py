#!/usr/bin/env python3

import subprocess
import logging
import argparse
import os
import sys
import glob

def setup_logging():
    logging.basicConfig(filename='annotation.log', filemode='a', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def run_subprocess(cmd, log_file):
    logging.info(f'Running command: {" ".join(cmd)}')
    with open(log_file, 'a') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        logging.error(f'Command failed: {" ".join(cmd)}')
        sys.exit(f'Error: Command failed. Check {log_file} for details.')

def run_prokka(fasta_file, output_dir, prefix='annotation', prokka_options=None):
    """
    Run Prokka on a specified FASTA file.
    
    Parameters:
    - fasta_file (str): Path to the FASTA file.
    - output_dir (str): Directory where Prokka output should be stored.
    - prefix (str): Prefix for output files.
    - prokka_options (list): Additional command-line options for Prokka.
    """
    prokka_options = prokka_options or []
    # Add the --force option to allow overwriting existing directories
    prokka_options.append('--force')

    cmd = [
        'prokka',
        '--outdir', output_dir,
        '--prefix', prefix,
        fasta_file
    ] + prokka_options
    run_subprocess(cmd, 'prokka_output.log')

def annotate_contigs(filtered_contigs_dir, annotation_output_dir):
    """
    Annotate all contig files in the filtered contigs directory using Prokka.
    """
    os.makedirs(annotation_output_dir, exist_ok=True)
    fasta_files = glob.glob(os.path.join(filtered_contigs_dir, '*.fasta'))
    if not fasta_files:
        logging.warning(f"No FASTA files found in {filtered_contigs_dir} for annotation.")
        return
    for fasta_file in fasta_files:
        sample_name = os.path.basename(fasta_file).split('.')[0]
        sample_output_dir = os.path.join(annotation_output_dir, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        run_prokka(fasta_file, sample_output_dir, prefix=sample_name)

def main(args=None):
    parser = argparse.ArgumentParser(description='Genome Annotation Pipeline')
    parser.add_argument('--filtered_contigs_dir', required=True, help='Path to the filtered contigs directory')
    parser.add_argument('--output_dir', required=True, help='Path to the output directory')

    args = parser.parse_args(args)

    setup_logging()

    # The annotation directory is inferred from output_dir
    annotation_output_dir = os.path.join(args.output_dir, "Annotation")

    logging.info(f'Filtered contigs directory: {args.filtered_contigs_dir}')
    logging.info(f'Annotation output directory: {annotation_output_dir}')

    # Check that filtered_contigs_dir exists and is not empty
    if not os.path.isdir(args.filtered_contigs_dir):
        logging.error(f'Filtered contigs directory does not exist: {args.filtered_contigs_dir}')
        sys.exit(f'Error: Filtered contigs directory does not exist: {args.filtered_contigs_dir}')

    fasta_files = glob.glob(os.path.join(args.filtered_contigs_dir, '*.fasta'))
    if not fasta_files:
        logging.warning(f"No FASTA files found in {args.filtered_contigs_dir}. Nothing to annotate.")
        return

    # Annotate contigs
    annotate_contigs(args.filtered_contigs_dir, annotation_output_dir)

if __name__ == "__main__":
    main()
