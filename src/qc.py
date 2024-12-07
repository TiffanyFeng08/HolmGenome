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
    fastq_files = glob.glob(os.path.join(data_dir, '*.fastq*'))
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
    extensions = ['.fastq', '.fastq.gz']
    paired_files = {}
    for ext in extensions:
        full_suffix1 = suffix1 + ext
        full_suffix2 = suffix2 + ext
        logging.info(f'Looking for files ending with {full_suffix1} and {full_suffix2}')
        fastq_files = glob.glob(os.path.join(input_dir, f'*{full_suffix1}')) + \
                      glob.glob(os.path.join(input_dir, f'*{full_suffix2}'))
        logging.info(f'Found files: {fastq_files}')
        for file in fastq_files:
            base = os.path.basename(file)
            if base.endswith(full_suffix1):
                sample_name = base.replace(full_suffix1, '')
                paired_files.setdefault(sample_name, {})['R1'] = file
            elif base.endswith(full_suffix2):
                sample_name = base.replace(full_suffix2, '')
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
        # If skip_trim is True, we do NOT run process_samples().
        # We assume input_dir points to trimmed reads.
        # Just run FastQC on these reads and store in QC/Trim.

        trimmed_data_qc_dir = os.path.join(args.output_dir, 'QC', 'Trim')
        run_fastqc(args.fastqc_path, args.input_dir, trimmed_data_qc_dir)

    else:
        # Normal run: pair files, trim, then run fastqc on raw and trimmed reads
        paired_files = pair_fastq_files(args.input_dir, args.suffix1, args.suffix2)
        process_samples(args, paired_files)

        # After trimming, run fastqc on raw_data and trimmed_data
        raw_data_qc_dir = os.path.join(args.output_dir, 'QC', 'raw_data')
        trimmed_data_qc_dir = os.path.join(args.output_dir, 'QC', 'Trim')

        # FastQC on raw data (original input_dir)
        run_fastqc(args.fastqc_path, args.input_dir, raw_data_qc_dir)

        # FastQC on trimmed data: The trimmed reads are in output_dir/Trim_data
        # Suffix might differ for trimmed reads (_R1_paired, _R2_paired)
        # If needed, rerun qc.py with --skip_trim and adjusted suffixes on second run
        trimmed_data_path = os.path.join(args.output_dir, 'Trim_data')
        run_fastqc(args.fastqc_path, trimmed_data_path, trimmed_data_qc_dir)

    return

if __name__ == "__main__":
    main()
