import subprocess

def run_prokka(fasta_file, output_dir):
    """
    Run Prokka on a specified FASTA file.
    
    Parameters:
    fasta_file (str): Path to the FASTA file.
    output_dir (str): Directory where Prokka output should be stored.
    """
    try:
        # Command to run Prokka
        command = [
            'prokka',             # Prokka command
            '--outdir', output_dir,  # Output directory
            '--prefix', 'annotation',  # Prefix for output files
            fasta_file            # FASTA file for annotation
        ]

        # Execute the command
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Check if Prokka ran successfully
        if result.returncode == 0:
            print("Prokka finished successfully.")
            print(result.stdout)
        else:
            print("Prokka encountered an error:")
            print(result.stderr)

    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
if __name__ == "__main__":
    fasta_path = "/path/to/your/fasta_file.fasta"
    output_path = "/path/to/output/directory"
    run_prokka(fasta_path, output_path)
