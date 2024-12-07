#!/usr/bin/env python3
"""
HolmGenome.py

This script orchestrates the genome analysis pipeline, including quality control,
assembly, and annotation steps.

Usage:
    python HolmGenome.py [options]

Options:
    -i, --input              Input directory path
    -o, --output             Output directory path
    --trimmomatic_path       Path to the Trimmomatic executable or JAR
    --adapters_path          Path to the adapters file for Trimmomatic
    --prokka_db_path         Path to the Prokka database
    --min_contig_length      Minimum contig length (default: 1000)
    --check                  Check if all dependencies are installed
    -v, --version            Show pipeline version number and exit
    -h, --help               Show this help message and exit
"""

import subprocess
import sys
import os
import logging
import argparse
import shutil

script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(script_dir, 'src')
sys.path.append(src_dir)

from qc import main as qc_main
from assembly import main as assembly_main
from annotation import main as annotation_main

VERSION = "1.0.0"  # Set your pipeline version here

def show_version():
    print(f"HolmGenome pipeline version {VERSION}")

def setup_logging(log_level=logging.INFO, log_file='HolmGenome.log'):
    logging.basicConfig(
        filename=log_file,
        filemode='a',
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    console = logging.StreamHandler()
    console.setLevel(log_level)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

def check_tool(tool):
    tool_name = os.path.basename(tool)
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
            if result.returncode in [0,1]:
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
    parser.add_argument('--min_contig_length', default='1000', help='Minimum contig length (default: 1000)')
    parser.add_argument('--check', action='store_true', help='Check if all dependencies are installed')
    parser.add_argument('-v', '--version', action='store_true', help='Show the pipeline version number and exit')

    args = parser.parse_args()

    # If --version is specified, show version and exit immediately
    if args.version:
        show_version()
        sys.exit(0)

    # If --check is used, run check_required_tools and exit immediately
    if args.check:
        tools = [
            'fastqc',
            'java',  # For Trimmomatic
            'spades.py',
            'prokka',
            'checkm',
            'reformat.sh',
            'bbmap.sh',
            'quast'
        ]
        # Setup a minimal logging before running check
        setup_logging()
        check_required_tools(tools)
        logging.info("Dependencies check completed successfully.")
        sys.exit(0)

    # Only prompt for directories and paths if we're not in --check mode
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

    if not os.path.isdir(args.input):
        print(f"Error: Input directory does not exist: {args.input}")
        sys.exit(1)
    if not os.path.isdir(args.output):
        try:
            os.makedirs(args.output, exist_ok=True)
        except Exception as e:
            print(f"Error: Could not create output directory {args.output}: {e}")
            sys.exit(1)

    log_file = os.path.join(args.output, 'HolmGenome.log')
    setup_logging(log_file=log_file)

    logging.info('Starting HolmGenome pipeline with user-specified arguments.')

    # Run QC on raw data
    qc_args_raw = [
        '--input_dir', args.input,
        '--output_dir', args.output,
        '--trimmomatic_path', args.trimmomatic_path,
        '--adapters_path', args.adapters_path,
        '--suffix1', '_R1_001',
        '--suffix2', '_R2_001'
    ]
    logging.info('Starting Quality Control on raw data.')
    qc_main(qc_args_raw)
    logging.info('Quality Control on raw data completed successfully.')

    trimmed_data_path = os.path.join(args.output, 'Trim_data')
    qc_args_trimmed = [
        '--input_dir', trimmed_data_path,
        '--output_dir', args.output,
        '--trimmomatic_path', args.trimmomatic_path,
        '--adapters_path', args.adapters_path,
        '--suffix1', '_R1_paired',
        '--suffix2', '_R2_paired',
        '--skip_trim'
    ]
    logging.info('Starting Quality Control on trimmed reads.')
    qc_main(qc_args_trimmed)
    logging.info('Quality Control on trimmed reads completed successfully.')

    assembly_args = [
        '--output_dir', args.output,
        '--spades_path', 'spades.py',
        '--quast_path', 'quast',
        '--reformat_path', 'reformat.sh',
        '--minlength', str(args.min_contig_length)
    ]
    logging.info('Starting Assembly step.')
    assembly_main(assembly_args)
    logging.info('Assembly step completed successfully.')

    filtered_contigs_dir = os.path.join(args.output, 'Assembly', 'contigs', 'filtered_contigs')

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
