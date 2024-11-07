# src/annotation.py

import os
import subprocess
import logging
import argparse
import sys
import glob

def setup_logging():
    logging.basicConfig(filename='annotation.log', filemode='a', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def run_subprocess(cmd, log_file):
    logging.info(f'Running command: {" ".join(cmd)}')
    with open(log_file, 'a') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        logging.error(f'Command failed: {" ".join(cmd)}')
        sys.exit(f'Error: Command failed. Check {log_file} for details.')

def run_prokka(fasta_file, output_dir, prefix='annotation', prokka_path='prokka', prokka_options=None):
    """
    Run Prokka on a specified FASTA file.

    Parameters:
    - fasta_file (str): Path to the FASTA file.
    - output_dir (str): Directory where Prokka output should be stored.
    - prefix (str): Prefix for output files.
    - prokka_path (str): Path to the Prokka executable.
    - prokka_options (list): Additional options for Prokka.
    """
    prokka_options = prokka_options or []
    try:
        # Command to run Prokka
        command = [
            prokka_path,
            '--outdir', output_dir,
            '--prefix', prefix,
            fasta_file
        ] + prokka_options

        # Execute the command
        run_subprocess(command, 'prokka_output.log')

        logging.info(f"Prokka finished successfully for {fasta_file}")

    except Exception as e:
        logging.error(f"An error occurred while running Prokka: {e}")
        sys.exit(1)

def annotate_contigs(contigs_dir, annotation_output_dir, prokka_path='prokka'):
    """
    Run Prokka on all FASTA files in the contigs directory.

    Parameters:
    - contigs_dir (str): Directory containing contigs FASTA files.
    - annotation_output_dir (str): Directory to store Prokka outputs.
    - prokka_path (str): Path to the Prokka executable.
    """
    fasta_files = glob.glob(os.path.join(contigs_dir, '*.fasta'))
    if not fasta_files:
        logging.error(f"No FASTA files found in {contigs_dir}")
        sys.exit(1)
    for fasta_file in fasta_files:
        sample_name = os.path.splitext(os.path.basename(fasta_file))[0]
        sample_output_dir = os.path.join(annotation_output_dir, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        run_prokka(fasta_file, sample_output_dir, prefix=sample_name, prokka_path=prokka_path)

def main(args=None):
    parser = argparse.ArgumentParser(description='Annotation Pipeline using Prokka')
    parser.add_argument('--contigs_dir', required=True, help='Directory containing contigs FASTA files')
    parser.add_argument('--annotation_output_dir', required=True, help='Directory to store Prokka outputs')
    parser.add_argument('--prokka_path', default='prokka', help='Path to the Prokka executable')
    args = parser.parse_args(args)

    setup_logging()

    annotate_contigs(args.contigs_dir, args.annotation_output_dir, prokka_path=args.prokka_path)

if __name__ == "__main__":
    main()
