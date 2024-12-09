# HolmGenome Pipeline

**HolmGenome** is a one-step pipeline that performs **assembly**, **quality control**, and **annotation** of genomic data. It takes raw sequencing reads and produces a fully annotated genome in a single run.

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
```
python HolmGenome.py -i INPUT -o OUTPUT --trimmomatic_path TRIMMOMATIC_PATH --adapters_path ADAPTERS_PATH --prokka_db_path PROKKA_DB_PATH
```
