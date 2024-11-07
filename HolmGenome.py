# This script checks the installation status of the tools used in the pipeline.
import subprocess
# List of tools to check
tools = ["fastqc", "multiqc", "trimmomatic", "spades.py", "prokka", "checkm"]

def check_tool(tool):
    try:
        # Run the tool with --version or a version-check command
        result = subprocess.run([tool, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print(f"{tool} is installed: {result.stdout.decode().strip()}")
        else:
            print(f"{tool} is not installed or not accessible.")
    except FileNotFoundError:
        print(f"{tool} is not installed or not found in the system path.")

if __name__ == "__main__":
    print("Checking tool installation status:\n")
    for tool in tools:
        check_tool(tool)
