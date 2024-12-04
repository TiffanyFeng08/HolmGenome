#!/usr/bin/env python3
"""
HolmGenome.py

This script orchestrates the genome analysis pipeline, including quality control,
assembly, and annotation steps. It reads configuration parameters from a YAML file,
checks for the required tools, and executes each step while handling exceptions
and logging the process.

Usage:
    python HolmGenome.py --config config.yaml
"""

import subprocess
import sys
import os
import logging
import yaml
import argparse
import shutil

# Adjust the Python path to include the 'src' directory
script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(script_dir, 'src')
sys.path.append(src_dir)

# Import the main functions from qc.py, assembly.py, and annotation.py
from qc import main as qc_main
from assembly import main as assembly_main
from annotation import main as annotation_main

def setup_logging(log_level=logging.INFO):
    """
    Configure logging for the script.
    """
    logging.basicConfig(
        filename='HolmGenome.log',
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
    """
    Check if a tool is installed and accessible in the system PATH or at a specified path.

    Parameters:
    - tool (str): Name of the tool executable or full path to it.

    Returns:
    - None
    """
    tool_name = os.path.basename(tool)
    # Determine if the tool is in PATH or if a full path is provided
    if os.path.isabs(tool):
        tool_path = tool
    else:
        tool_path = shutil.which(tool)
        if tool_path is None:
            logging.error(f"{tool_name} is not installed or not found in the system PATH.")
            sys.exit(f"Error: {tool_name} is not installed or not found in the system PATH.")

    if not os.access(tool_path, os.X_OK):
        logging.error(f"{tool_name} is not executable.")
        sys.exit(f"Error: {tool_name} is not executable.")

    commands_to_try = [['--version'], ['-v'], ['-h'], ['--help']]
    for cmd_args in commands_to_try:
        try:
            with open(os.devnull, 'w') as devnull:
                result = subprocess.run(
                    [tool_path] + cmd_args,
                    stdout=devnull,
                    stderr=devnull
                )
            if result.returncode in [0, 1]:  # Some tools return 1 for help/version
                logging.info(f"{tool_name} is installed and accessible.")
                return
        except Exception as e:
            logging.error(f"Error checking {tool_name}: {e}")
            continue
    logging.error(f"{tool_name} is not functioning as expected.")
    sys.exit(f"Error: {tool_name} is not functioning as expected.")

def check_required_tools(tools):
    """
    Check all required tools for the pipeline.

    Parameters:
    - tools (list): List of tool names or paths to check.

    Returns:
    - None
    """
    logging.info('Checking required tools...')
    for tool in tools:
        check_tool(tool)
    logging.info('All required tools are installed and accessible.')

def main():
    parser = argparse.ArgumentParser(description='HolmGenome Pipeline')
    parser.add_argument('--config', help='Path to the configuration file')
    parser.add_argument('--input_dir', help='Path to the input directory')
    parser.add_argument('--output_dir', help='Path to the output directory')
    parser.add_argument('--log_level', default='INFO', help='Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)')
    args = parser.parse_args()

    # Setup logging
    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        print(f'Invalid log level: {args.log_level}')
        sys.exit('Error: Invalid log level.')
    setup_logging(log_level=numeric_level)

    logging.info('Starting HolmGenome pipeline.')

    # List of tools to check
    tools = [
        'fastqc',
        'java',  # For Trimmomatic
        'spades.py',
        'prokka',
        'checkm',
        'reformat.sh',  # Use system-installed version from BBMap module
        'bbmap.sh',
        'quast'
    ]

    # Check for required tools before doing anything else
    logging.info('Checking for required tools...')
    check_required_tools(tools)

    # Load configuration
    config = load_config(args.config)

    # Merge command-line arguments into the configuration
    if args.input_dir:
        config['input_dir'] = args.input_dir
    if args.output_dir:
        config['output_dir'] = args.output_dir

    # If no arguments are provided, prompt for input and output directories
    if not config.get('input_dir') or not config.get('output_dir'):
        print("Required directories not specified. Please provide them now.")
        if not config.get('input_dir'):
            config['input_dir'] = input("Enter the input directory: ").strip()
        if not config.get('output_dir'):
            config['output_dir'] = input("Enter the output directory: ").strip()

    # Verify required parameters are present
    required_params = ['input_dir', 'output_dir']
    missing_params = [p for p in required_params if p not in config or not config[p]]
    if missing_params:
        logging.error(f"Missing required parameters: {', '.join(missing_params)}")
        sys.exit(f"Error: Missing required parameters: {', '.join(missing_params)}")

    try:
        # Prepare arguments for qc.py
        
        qc_args = [
            '--input_dir', config['input_dir'],
            '--output_dir', config['output_dir'],
            '--trimmomatic_path', config['trimmomatic_path'],
            '--adapters_path', config['adapters_path'],
            '--suffix1', '_R1_001',
            '--suffix2', '_R2_001'
        ]

        logging.info('Starting Quality Control step.')
        # Run Quality Control step
        qc_main(qc_args)
        logging.info('Quality Control step completed successfully.')

        # **Set the trimmed data path based on the QC output**
        trimmed_data_path = os.path.join(config['output_dir'], 'Trim_data')

        # Prepare arguments for assembly.py
        # Replace '--base_dir' with '--output_dir'
    
        assembly_args = [
            '--output_dir', config['output_dir'],
            '--trimmed_data_path', trimmed_data_path,
            '--spades_path', config.get('spades_path', 'spades.py'),
            '--minlength', str(config.get('min_contig_length', 1000))
        ]


        logging.info('Starting Assembly step.')
        # Run Assembly step
        assembly_main(assembly_args)
        logging.info('Assembly step completed successfully.')

        # Prepare arguments for annotation.py
        annotation_args = [
            '--filtered_contigs_dir', config.get('filtered_contigs_dir', 'path/to/filtered_contigs'),
            '--annotation_output_dir', config.get('annotation_output_dir', 'path/to/annotation_output')
        ]

        logging.info('Starting Annotation step.')
        # Run Annotation step
        annotation_main(annotation_args)
        logging.info('Annotation step completed successfully.')

        logging.info('HolmGenome pipeline completed successfully.')

    except Exception as e:
        logging.exception('An exception occurred during pipeline execution.')
        sys.exit(f"Error: {e}")

def load_config(config_file=None):
    """
    Load configuration parameters from a YAML file.

    Parameters:
    - config_file (str): Path to the configuration YAML file.

    Returns:
    - config (dict): Dictionary containing configuration parameters.
    """
    config = {}
    if config_file:
        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
            logging.info('Configuration loaded successfully.')
        except FileNotFoundError:
            logging.error(f'Configuration file {config_file} not found.')
            sys.exit(f"Error: Configuration file {config_file} not found.")
        except yaml.YAMLError as e:
            logging.error(f"Error parsing configuration file: {e}")
            sys.exit("Error: Invalid YAML syntax in configuration file.")
    return config

if __name__ == "__main__":
    main()
