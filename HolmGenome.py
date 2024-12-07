#!/usr/bin/env python3
"""
HolmGenome.py

This script orchestrates the genome analysis pipeline, including quality control,
assembly, and annotation steps. It uses command-line arguments instead of a YAML file.
Now it calls qc_main twice to run QC on trimmed reads as well.
Usage:
    python HolmGenome.py [options]

Options:
    -i, --input              Input directory path
    -o, --output             Output directory path
    --trimmomatic_path       Path to the Trimmomatic executable or JAR
    --adapters_path          Path to the adapters file for Trimmomatic
    --prokka_db_path         Path to the Prokka database
    --min_contig_length      Minimum contig length (default: 100)
    --check                  Check required tools and exit
    -h, --help               Show this help message and exit

"""

import subprocess
import sys
import os
import logging
import argparse
import shutil

# Adjust the Python path to include the 'src' directory
script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(script_dir, 'src')
sys.path.append(src_dir)

from qc import main as qc_main
from assembly import main as assembly_main
from annotation import main as annotation_main

def setup_logging(log_level=logging.INFO, log_file='HolmGenome.log'):
    """
    Configure logging for the script.
    """
    logging.basicConfig(
        filename=log_file,
        filemode='a',
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    # Also log to console
    console = logging.StreamHandler()
    console.setLevel(log_level)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

def check_tool(tool):
    tool_name = os.path.basename(tool)
    if os.path.isabs(tool):
        tool_path = tool
    else:
        tool_path = shutil.which(tool)
        if tool_path is None:
            print(f"Error: {tool_name} not found in PATH.")
            sys.exit(1)

    if not os.access(tool_path, os.X_OK):
        print(f"Error: {tool_name} is not executable.")
        sys.exit(1)

    commands_to_try = [['--version'], ['-v'], ['-h'], ['--help']]
    for cmd_args in commands_to_try:
        try:
            with open(os.devnull, 'w') as devnull:
                result = subprocess.run([tool_path] + cmd_args, stdout=devnull, stderr=devnull)
            if result.returncode in [0, 1]:
                logging.info(f"{tool_name} is installed and accessible.")
                return
        except Exception as e:
            logging.error(f"Error checking {tool_name}: {e}")
            continue
    logging.error(f"{tool_name} not functioning as expected.")
    sys.exit(f"Error: {tool_name} not functioning as expected.")

def check_required_tools(tools):
    logging.info('Checking required tools...')
    for tool in tools:
        check_tool(tool)
    logging.info('All required tools are installed and accessible.')

def main():
    parser = argparse.ArgumentParser(description='HolmGenome Pipeline')
    parser.add_argument('-i', '--input', help='Path to the input directory')
    parser.add_argument('-o', '--output', help='Path to the output directory')
    parser.add_argument('--trimmomatic_path', help='Path to the Trimmomatic executable or JAR')
    parser.add_argument('--adapters_path', help='Path to the adapters file')
    parser.add_argument('--prokka_db_path', help='Path to the Prokka database')
    parser.add_argument('--fastqc_path', help='Path to the FastQC executable')
    parser.add_argument('--min_contig_length', default='100', help='Minimum contig length (default: 100)')
    parser.add_argument('--check', action='store_true', help='Check required tools and exit')
    parser.add_argument('--log_level', default='INFO', help='Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)')

    args = parser.parse_args()

    # Prompt for missing required arguments if not provided
    if not args.input:
        args.input = input("Enter the input directory: ").strip()
    if not args.output:
        args.output = input("Enter the output directory: ").strip()
    if not args.trimmomatic_path:
        args.trimmomatic_path = input("Enter the Trimmomatic path: ").strip()
    if not args.adapters_path:
        args.adapters_path = input("Enter the adapters file path: ").strip()
    if not args.prokka_db_path:
        args.prokka_db_path = input("Enter the Prokka database path: ").strip()
    if not args.fastqc_path:
        args.fastqc_path = input("Enter the FastQC path: ").strip()

    # Ensure directories exist
    if not os.path.isdir(args.input):
        print(f"Error: Input directory does not exist: {args.input}")
        sys.exit(1)
    if not os.path.isdir(args.output):
        # Attempt to create output directory if it doesn't exist
        try:
            os.makedirs(args.output, exist_ok=True)
        except Exception as e:
            print(f"Error: Could not create output directory {args.output}: {e}")
            sys.exit(1)

    # Set up logging in the output directory
    log_file = os.path.join(args.output, 'HolmGenome.log')
    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        print(f'Invalid log level: {args.log_level}')
        sys.exit('Error: Invalid log level.')
    setup_logging(log_level=numeric_level, log_file=log_file)

    logging.info('Starting HolmGenome pipeline with user-specified arguments.')

    # If --check is used, just check required tools and exit
    if args.check:
        tools = [
            'fastqc',
            'java',  # For Trimmomatic
            'spades.py',
            'prokka',
            'checkm',
            'reformat.sh',  # from BBMap suite
            'bbmap.sh',
            'quast'
        ]
        check_required_tools(tools)
        logging.info("Tool check completed successfully. Exiting as requested by --check.")
        sys.exit(0)

    # Proceed with the pipeline
    # Prepare arguments for qc.py (first run on raw data)
    qc_args = [
        '--input_dir', args.input,
        '--output_dir', args.output,
        '--trimmomatic_path', args.trimmomatic_path,
        '--adapters_path', args.adapters_path,
        '--fastqc_path', args.fastqc_path,
        '--suffix1', '_R1_001',
        '--suffix2', '_R2_001'
    ]

    logging.info('Starting Quality Control step on raw data.')
    qc_main(qc_args)
    logging.info('Quality Control on raw data completed successfully.')

    # Set the trimmed data path based on the QC output
    trimmed_data_path = os.path.join(args.output, 'Trim_data')

    # Run qc_main again for trimmed reads
    # This second call runs fastqc on trimmed data. We don't need trimming again, so no trimmomatic args needed.
    # Just ensure qc.py knows to run fastqc on existing trimmed data.
    # If qc.py requires special arguments to skip trimming, add them, otherwise assume it runs fastqc only.
    qc_args_trimmed = [
        '--input_dir', trimmed_data_path,
        '--output_dir', args.output,
        '--fastqc_path', args.fastqc_path
    ]

    logging.info('Starting Quality Control step on trimmed reads.')
    qc_main(qc_args_trimmed)
    logging.info('Quality Control on trimmed reads completed successfully.')

    # Prepare arguments for assembly.py
    assembly_args = [
        '--output_dir', args.output,
        '--spades_path', 'spades.py',  # or args.spades_path if needed
        '--quast_path', 'quast',
        '--reformat_path', 'reformat.sh',
        '--minlength', str(args.min_contig_length)
    ]

    logging.info('Starting Assembly step.')
    assembly_main(assembly_args)
    logging.info('Assembly step completed successfully.')

    filtered_contigs_dir = os.path.join(args.output, 'Assembly', 'contigs', 'filtered_contigs')

    # Prepare arguments for annotation.py
    annotation_args = [
        '--filtered_contigs_dir', filtered_contigs_dir,
        '--output_dir', args.output
    ]

    logging.info('Starting Annotation step.')
    annotation_main(annotation_args)
    logging.info('Annotation step completed successfully.')

    logging.info('HolmGenome pipeline completed successfully.')

if __name__ == "__main__":
    main()
