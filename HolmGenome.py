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

def load_config(config_file='config.yaml'):
    """
    Load configuration parameters from a YAML file.

    Parameters:
    - config_file (str): Path to the configuration YAML file.

    Returns:
    - config (dict): Dictionary containing configuration parameters.
    """
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        logging.info('Configuration loaded successfully.')
        return config
    except FileNotFoundError:
        logging.error(f'Configuration file {config_file} not found.')
        sys.exit(f'Error: Configuration file {config_file} not found.')
    except yaml.YAMLError as e:
        logging.error(f'Error parsing configuration file: {e}')
        sys.exit('Error: Invalid YAML syntax in configuration file.')

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
    logging.info('Tool check completed.')

def main():
    parser = argparse.ArgumentParser(description='HolmGenome Pipeline')
    parser.add_argument('--config', default='config.yaml', help='Path to the configuration file')
    parser.add_argument('--log_level', default='INFO', help='Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)')
    args = parser.parse_args()

    # Setup logging
    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        print(f'Invalid log level: {args.log_level}')
        sys.exit('Error: Invalid log level.')
    setup_logging(log_level=numeric_level)

    logging.info('Starting HolmGenome pipeline.')

    # Load configuration
    config = load_config(args.config)

    # Full path to reformat.sh (if using local copy)
    # reformat_path = os.path.join(src_dir, 'reformat.sh')

    # List of tools to check
    tools = [
        config.get('fastqc_path', 'fastqc'),
        'java',  # For Trimmomatic
        'spades.py',
        'prokka',
        'checkm',
        'reformat.sh',  # Use system-installed version from BBMap module
        'bbmap.sh',
        'quast'
    ]

    # Check for required tools
    check_required_tools(tools)

    try:
        # Prepare arguments for qc.py
        qc_args = [
            '--input_dir', config['input_dir'],
            '--output_dir', config['output_dir'],
            '--trimmomatic_path', config['trimmomatic_path'],
            '--adapters_path', config['adapters_path'],
            '--fastqc_path', config.get('fastqc_path', 'fastqc'),
            '--prokka_path', config.get('prokka_path', 'prokka')
            # Include additional arguments from config if needed
        ]

        logging.info('Starting Quality Control step.')
        # Run Quality Control step
        qc_main(qc_args)
        logging.info('Quality Control step completed successfully.')

        # Prepare arguments for assembly.py
        assembly_args = [
            '--base_dir', config['base_dir'],
            '--spades_path', config.get('spades_path', 'spades.py'),
            '--quast_path', config.get('quast_path', 'quast'),
            '--bbmap_path', config.get('bbmap_path', 'bbmap.sh'),
            '--reformat_path', config.get('reformat_path', 'reformat.sh'),
            '--minlength', str(config.get('min_contig_length', 1000))
            # Include additional arguments from config if needed
        ]

        logging.info('Starting Assembly step.')
        # Run Assembly step
        assembly_main(assembly_args)
        logging.info('Assembly step completed successfully.')

        # Prepare arguments for annotation.py
        annotation_args = [
            '--filtered_contigs_dir', config['filtered_contigs_dir'],
            '--annotation_output_dir', config['annotation_output_dir']
            # Include additional arguments from config if needed
        ]

        logging.info('Starting Annotation step.')
        # Run Annotation step
        annotation_main(annotation_args)
        logging.info('Annotation step completed successfully.')

        logging.info('HolmGenome pipeline completed successfully.')

    except Exception as e:
        logging.exception('An exception occurred during pipeline execution.')
        sys.exit(f'Error: {e}')

if __name__ == "__main__":
    main()
