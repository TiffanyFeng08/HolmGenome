# HolmGenome.py

import subprocess
import sys
import os

# Import the main function from QC.py
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))
from QC import main as qc_main

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
    
    # Define arguments for QC.py
    qc_args = [
        '--input_dir', 'path/to/input',          # Replace with actual input directory
        '--output_dir', 'path/to/output',        # Replace with actual output directory
        '--trimmomatic_path', '/path/to/trimmomatic.jar',
        '--adapters_path', '/path/to/adapters.fa',
        '--fastqc_path', '/usr/local/bin/fastqc',  # Replace if necessary
        '--multiqc_path', '/usr/local/bin/multiqc' # Replace if necessary
        # Include additional arguments as needed
    ]

    # Call the QC main function with the arguments
    qc_main(qc_args)

if __name__ == "__main__":
    main()
