# HolmGenome Pipeline

**HolmGenome** is a one-step pipeline that performs **assembly**, **quality control**, and **annotation** of bacteria genomic data. It takes raw sequencing reads and produces a fully annotated genome in a single run.

## Features
- **Assembly:** Uses SPAdes to assemble raw reads.
- **QC:** Employs Trimmomatic & FastQC for trimming and quality checks.
- **Annotation:** Runs Prokka to annotate assembled contigs.

## Installation
### Requirements
- **Python Version:** Python 3.6 or higher
- **External Dependencies:**
- fastqc
- trimmomatic
- spades.py
- prokka
- checkm
- quast
- bbmap

### Development Version
```bash
git clone https://github.com/TiffanyFeng08/HolmGenome.git

``` 
## Usage 
HolmGenome.py will detect all *.fastq or *.fastq.gz files in a directory and run the assembly pipeline on each sample it can pair

```
Usage: HolmGenome.py [OPTIONS]
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input directory
  -o OUTPUT, --output OUTPUT
                        Path to the output directory
  --trimmomatic_path TRIMMOMATIC_PATH
                        Path to the Trimmomatic executable or JAR
  --adapters_path ADAPTERS_PATH
                        Path to the adapters file
  --prokka_db_path PROKKA_DB_PATH
                        Path to the Prokka database
  --min_contig_length MIN_CONTIG_LENGTH
                        Minimum contig length (default: 1000)
  --check               Check if all dependencies are installed
  -v, --version         Show the pipeline version number and exit

``` 

```
python HolmGenome.py -i INPUT -o OUTPUT --trimmomatic_path TRIMMOMATIC_PATH --adapters_path ADAPTERS_PATH --prokka_db_path PROKKA_DB_PATH
```
### Example bash script
```
#!/bin/bash
#SBATCH --account=your_account
#SBATCH --time=5:00:00
#SBATCH --job-name=HolmGenome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=your.email@somemail.com
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

module load StdEnv/2020 gcc/9.3.0
module load python/3.10.2
module load scipy-stack/2022a
module load fastqc/0.12.0
module load trimmomatic/0.39
module load spades/3.15.4
module load prokka/1.14.5
module load bbmap/38.86
module load quast/5.2.0
source ~/HolmGenome/bin/activate

python HolmGenome.py -i /path/to/data -o /path/to/output \
--trimmomatic_path /path/to/trimmomatic-0.39.jar \
--adapters_path /path/to/adapters.fa \
--prokka_db_path /path/to/prokka_db 
```