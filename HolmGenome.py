# HolmGenome.py

import subprocess
import sys
import os

# Adjust sys.path to include 'src' directory
current_dir = os.path.dirname(__file__)
src_dir = os.path.join(current_dir, 'src')
sys.path.append(src_dir)

# Import the main functions from QC.py and Assembly.py
from QC import main as qc_main
from Assembly import main as assembly_main

# List of tools to check
tools = ["fastqc", "multiqc", "trimmomatic", "spades.py", "prokka", "checkm"]

def check_tool(tool):
    try:
        result = subprocess.run([tool, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print(f"{tool} is installed: {result.stdout.decode().strip()}")
        else:
            print(f"{tool} is not installed or not accessible.")
    except FileNotFoundError:
        print(f"{tool} is not installed or not found in the system path.")

def main():
    print("Checking tool installation status:\n")
    for tool in tools:
        check_tool(tool)
    
    # Define base directories
    input_dir = '/path/to/input'         # Replace with actual input directory
    output_dir = '/path/to/output'       # Replace with actual output directory

    # Define arguments for QC.py
    qc_args = [
        '--input_dir', input_dir,
        '--output_dir', output_dir,
        '--trimmomatic_path', '/path/to/trimmomatic.jar',
        '--adapters_path', '/path/to/adapters.fa',
        '--fastqc_path', '/usr/local/bin/fastqc',
        '--multiqc_path', '/usr/local/bin/multiqc',
        # Include additional arguments as needed
    ]

    # Call the QC main function with the arguments
    qc_main(qc_args)

    # Define arguments for Assembly.py
    assembly_args = [
        '--base_dir', output_dir,  # Use the same output_dir from QC.py
        '--spades_path', '/usr/local/bin/spades.py',
        '--reformat_path', '/path/to/reformat.sh',
        '--quast_path', '/usr/local/bin/quast',
        '--multiqc_path', '/usr/local/bin/multiqc',
        '--bbmap_path', '/path/to/bbmap.sh',
        '--minlength', '1000',
        # Include additional arguments as needed
    ]

    # Call the Assembly main function with the arguments
    assembly_main(assembly_args)

if __name__ == "__main__":
    main()
