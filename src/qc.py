# src/QC.py

import argparse
import os
import subprocess
import logging
import glob
import sys

def setup_logging():
    logging.basicConfig(filename='trim.log', filemode='a', level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def setup_directories(base_dir):
    dirs = [
        os.path.join(base_dir, 'Trim_data'),
        os.path.join(base_dir, 'QC', 'raw_data'),
        os.path.join(base_dir, 'QC', 'Trim')
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)

def run_subprocess(cmd, log_file):
    logging.info(f'Running command: {cmd}')
    with open(log_file, 'a') as f:
        result = subprocess.run(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        logging.error(f'Command failed: {cmd}')
        sys.exit(f'Error: Command failed. Check {log_file} for details.')

def run_fastqc(fastqc_path, data_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    # Search recursively for fastq files
    fastq_files = glob.glob(os.path.join(data_dir, '**', '*.fastq*'), recursive=True)
    if not fastq_files:
        logging.warning(f"No FASTQ files found in {data_dir} for FastQC.")
    for fastq_file in fastq_files:
        cmd = f'{fastqc_path} -o {output_dir} {fastq_file}'
        run_subprocess(cmd, 'fastqc_output.log')

def run_multiqc(multiqc_path, input_dir):
    os.makedirs(input_dir, exist_ok=True)
    cmd = f'{multiqc_path} -o {input_dir} {input_dir}'
    run_subprocess(cmd, 'multiqc_output.log')

def run_trimmomatic(params):
    cmd = (
        f"java -jar {params['trimmomatic_path']} PE -phred33 "
        f"{params['input_file1']} {params['input_file2']} "
        f"{params['output_file1_paired']} {params['output_file1_unpaired']} "
        f"{params['output_file2_paired']} {params['output_file2_unpaired']} "
        f"ILLUMINACLIP:{params['adapters_path']}:{params['illuminaclip']} "
        f"LEADING:{params['leading']} TRAILING:{params['trailing']} "
        f"SLIDINGWINDOW:{params['slidingwindow']} MINLEN:{params['minlen']}"
    )
    if params['crop']:
        cmd += f" CROP:{params['crop']}"
    if params['headcrop']:
        cmd += f" HEADCROP:{params['headcrop']}"
    run_subprocess(cmd, 'trimmomatic_output.log')

def run_bbduk(params):
    cmd = (
        f"{params['bbduk_path']} in1={params['input_file1']} in2={params['input_file2']} "
        f"out1={params['output_file1']} out2={params['output_file2']} "
        f"ref={params['adapters_path']} ktrim={params['ktrim']} k={params['k']} "
        f"mink={params['mink']} hdist={params['hdist']} tpe={params['tpe']} "
        f"tbo={params['tbo']} qtrim={params['qtrim']} trimq={params['trimq']} "
        f"minlen={params['minlen']}"
    )
    run_subprocess(cmd, 'bbduk_output.log')

def pair_fastq_files(input_dir, suffix1, suffix2):
    # Recursively find files that match the suffixes
    # Using '**' and recursive=True to handle subdirectories
    paired_files = {}
    # Consider .fastq and .fastq.gz
    patterns = [f"**/*{suffix1}*.fastq*", f"**/*{suffix2}*.fastq*"]
    found_files = []
    for pattern in patterns:
        found_files.extend(glob.glob(os.path.join(input_dir, pattern), recursive=True))
    found_files = list(set(found_files))
    logging.info(f'Found files with suffixes {suffix1}, {suffix2}: {found_files}')

    # Group by sample
    for file in found_files:
        base = os.path.basename(file)
        if suffix1 in base:
            sample_name = base.replace(suffix1, '').replace('.fastq.gz', '').replace('.fastq','')
            paired_files.setdefault(sample_name, {})['R1'] = file
        elif suffix2 in base:
            sample_name = base.replace(suffix2, '').replace('.fastq.gz', '').replace('.fastq','')
            paired_files.setdefault(sample_name, {})['R2'] = file

    logging.info(f'Paired files: {paired_files}')
    return paired_files

def process_samples(args, paired_files):
    for sample, files in paired_files.items():
        if 'R1' in files and 'R2' in files:
            input_file1 = files['R1']
            input_file2 = files['R2']
            base_output_path = os.path.join(args.output_dir, 'Trim_data', sample)
            os.makedirs(base_output_path, exist_ok=True)

            if args.bbduk:
                params = {
                    'input_file1': input_file1,
                    'input_file2': input_file2,
                    'output_file1': os.path.join(base_output_path, f'{sample}_R1_trimmed.fastq.gz'),
                    'output_file2': os.path.join(base_output_path, f'{sample}_R2_trimmed.fastq.gz'),
                    'bbduk_path': args.bbduk_path if args.bbduk_path else 'bbduk.sh',
                    'adapters_path': args.adapters_path,
                    'ktrim': args.ktrim,
                    'k': args.k,
                    'mink': args.mink,
                    'hdist': args.hdist,
                    'tpe': args.tpe,
                    'tbo': args.tbo,
                    'qtrim': args.qtrim,
                    'trimq': args.trimq,
                    'minlen': args.minlen
                }
                run_bbduk(params)
            else:
                params = {
                    'input_file1': input_file1,
                    'input_file2': input_file2,
                    'output_file1_paired': os.path.join(base_output_path, f"{sample}_R1_paired.fastq.gz"),
                    'output_file1_unpaired': os.path.join(base_output_path, f"{sample}_R1_unpaired.fastq.gz"),
                    'output_file2_paired': os.path.join(base_output_path, f"{sample}_R2_paired.fastq.gz"),
                    'output_file2_unpaired': os.path.join(base_output_path, f"{sample}_R2_unpaired.fastq.gz"),
                    'trimmomatic_path': args.trimmomatic_path,
                    'adapters_path': args.adapters_path,
                    'illuminaclip': args.illuminaclip,
                    'slidingwindow': args.slidingwindow,
                    'leading': args.leading,
                    'trailing': args.trailing,
                    'crop': args.crop,
                    'headcrop': args.headcrop,
                    'minlen': args.minlen
                }
                run_trimmomatic(params)
        else:
            logging.warning(f'Missing pair for sample {sample}')

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', required=True, help='Path to the input directory')
    parser.add_argument('--output_dir', required=True, help='Path to the output directory')
    parser.add_argument('--trimmomatic_path', required=True, help='Path to the Trimmomatic jar file')
    parser.add_argument('--adapters_path', required=True, help='Path to the adapters file')
    parser.add_argument('--illuminaclip', default='2:30:10', help='ILLUMINACLIP option for Trimmomatic')
    parser.add_argument('--slidingwindow', default='4:15', help='SLIDINGWINDOW option for Trimmomatic')
    parser.add_argument('--leading', default='3', help='LEADING option for Trimmomatic')
    parser.add_argument('--trailing', default='3', help='TRAILING option for Trimmomatic')
    parser.add_argument('--crop', default=None, help='CROP option for Trimmomatic')
    parser.add_argument('--headcrop', default=None, help='HEADCROP option for Trimmomatic')
    parser.add_argument('--minlen', default='36', help='MINLEN option')
    parser.add_argument('--suffix1', default='_R1_001', help='Suffix for the first file in each pair')
    parser.add_argument('--suffix2', default='_R2_001', help='Suffix for the second file in each pair')
    parser.add_argument('--bbduk', action='store_true', help='Use BBDuk for processing')
    parser.add_argument('--bbduk_path', required=False, help='Path to the BBDuk executable')
    parser.add_argument('--ktrim', default='r', help='KTRIM option for BBDuk')
    parser.add_argument('--k', default='23', help='K option for BBDuk')
    parser.add_argument('--mink', default='11', help='MINK option for BBDuk')
    parser.add_argument('--hdist', default='1', help='HDIST option for BBDuk')
    parser.add_argument('--tpe', default='t', help='TPE option for BBDuk')
    parser.add_argument('--tbo', default='t', help='TBO option for BBDuk')
    parser.add_argument('--qtrim', default='rl', help='QTRIM option for BBDuk')
    parser.add_argument('--trimq', default='10', help='TRIMQ option for BBDuk')
    parser.add_argument('--fastqc_path', default='fastqc', help='Path to the FastQC executable')
    parser.add_argument('--skip_trim', action='store_true', help='Skip trimming and just run FastQC on input_dir')

    args = parser.parse_args(argv)

    setup_logging()
    setup_directories(args.output_dir)

    if args.skip_trim:
        # Just run FastQC on the given input_dir (assume trimmed reads here)
        trimmed_data_qc_dir = os.path.join(args.output_dir, 'QC', 'Trim')
        run_fastqc(args.fastqc_path, args.input_dir, trimmed_data_qc_dir)
    else:
        # Normal operation: pair files, trim, run fastqc on raw and trimmed
        paired_files = pair_fastq_files(args.input_dir, args.suffix1, args.suffix2)
        process_samples(args, paired_files)

        # After trimming, run fastqc on raw_data and trimmed_data
        raw_data_qc_dir = os.path.join(args.output_dir, 'QC', 'raw_data')
        trimmed_data_qc_dir = os.path.join(args.output_dir, 'QC', 'Trim')

        # FastQC on raw data
        run_fastqc(args.fastqc_path, args.input_dir, raw_data_qc_dir)

        # For trimmed reads in subdirectories with `_R1_paired` and `_R2_paired`
        # If needed, run qc.py again with --skip_trim and updated suffixes.
        # But we can also just rely on run_fastqc to find fastqs recursively:
        trimmed_data_path = os.path.join(args.output_dir, 'Trim_data')
        run_fastqc(args.fastqc_path, trimmed_data_path, trimmed_data_qc_dir)

    return

if __name__ == "__main__":
    main()
