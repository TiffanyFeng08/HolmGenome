# src/assembly.py

import os
import subprocess
import shutil
import argparse
import glob
import logging
import sys

def setup_logging():
    logging.basicConfig(filename='assembly.log', filemode='a', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def setup_directories(directories):
    for directory in directories:
        os.makedirs(directory, exist_ok=True)

def run_subprocess(cmd, log_file):
    logging.info(f'Running command: {" ".join(cmd)}')
    with open(log_file, 'a') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        logging.error(f'Command failed: {" ".join(cmd)}')
        sys.exit(f'Error: Command failed. Check {log_file} for details.')

def pair_fastq_files(input_dir, suffix1, suffix2):
    fastq_files = glob.glob(os.path.join(input_dir, f'*{suffix1}')) + \
                  glob.glob(os.path.join(input_dir, f'*{suffix2}'))
    paired_files = {}
    for file in fastq_files:
        base = os.path.basename(file)
        if base.endswith(suffix1):
            sample_name = base.replace(suffix1, '')
            paired_files.setdefault(sample_name, {})['R1'] = file
        elif base.endswith(suffix2):
            sample_name = base.replace(suffix2, '')
            paired_files.setdefault(sample_name, {})['R2'] = file
    return paired_files

def run_assembly(trimmed_data_path, assembly_dir, spades_path='spades.py'):
    paired_files = pair_fastq_files(trimmed_data_path, '_R1_001_paired.fastq.gz', '_R2_001_paired.fastq.gz')
    for sample, files in paired_files.items():
        if 'R1' in files and 'R2' in files:
            sample_output_dir = os.path.join(assembly_dir, "contigs", "samples", sample)
            os.makedirs(sample_output_dir, exist_ok=True)
            cmd = [
                spades_path,
                '--pe1-1', files['R1'],
                '--pe1-2', files['R2'],
                '-o', sample_output_dir
            ]
            run_subprocess(cmd, 'spades_output.log')
            # Copy the resulting contigs to the all_contigs_dir
            contigs_src = os.path.join(sample_output_dir, "contigs.fasta")
            contigs_dest = os.path.join(assembly_dir, "contigs", "all_contigs", f"{sample}.fasta")
            if os.path.exists(contigs_src):
                shutil.copy(contigs_src, contigs_dest)
            else:
                logging.warning(f"Contigs file not found for sample {sample}")
        else:
            logging.warning(f"Paired files not found for sample {sample}")

def reformat_contigs(all_contigs_dir, filtered_contigs_dir, minlength=1000, reformat_path='reformat.sh'):
    contig_files = glob.glob(os.path.join(all_contigs_dir, '*.fasta'))
    for file in contig_files:
        sample_name = os.path.basename(file).split('.')[0]
        output_file = os.path.join(filtered_contigs_dir, f"{sample_name}_filtered_contigs.fasta")
        cmd = [
            reformat_path,
            f"in={file}",
            f"out={output_file}",
            f"minlength={minlength}"
        ]
        run_subprocess(cmd, 'reformat_output.log')

def run_quast(contigs_dir, output_dir, quast_path='quast'):
    contig_files = glob.glob(os.path.join(contigs_dir, '*.fasta'))
    if not contig_files:
        logging.warning(f"No contig files found in {contigs_dir}")
        return
    cmd = [quast_path] + contig_files + ['-o', output_dir]
    run_subprocess(cmd, 'quast_output.log')

def run_bbmap(trimmed_data_path, filtered_contigs_dir, coverage_dir, bbmap_path='bbmap.sh'):
    paired_files = pair_fastq_files(trimmed_data_path, '_R1_001_paired.fastq.gz', '_R2_001_paired.fastq.gz')
    for sample, files in paired_files.items():
        if 'R1' in files and 'R2' in files:
            ref_file = os.path.join(filtered_contigs_dir, f"{sample}_filtered_contigs.fasta")
            if not os.path.exists(ref_file):
                logging.warning(f"Reference contigs not found for sample {sample}")
                continue
            sample_coverage_dir = os.path.join(coverage_dir, sample)
            os.makedirs(sample_coverage_dir, exist_ok=True)
            cmd = [
                bbmap_path,
                f"in1={files['R1']}",
                f"in2={files['R2']}",
                f"ref={ref_file}",
                f"covstats={os.path.join(sample_coverage_dir, sample + '_covstats.txt')}",
                f"covhist={os.path.join(sample_coverage_dir, sample + '_covhist.tsv')}",
                f"basecov={os.path.join(sample_coverage_dir, sample + '_basecov.txt')}",
                "usejni=t",
                f"out={os.path.join(sample_coverage_dir, sample + '_mapped.bam')}"
            ]
            run_subprocess(cmd, 'bbmap_output.log')
        else:
            logging.warning(f"Paired files not found for sample {sample}")

def main(args=None):
    parser = argparse.ArgumentParser(description='Genome Assembly Pipeline')
    parser.add_argument('--trimmed_data_path', required=True, help='Path to the trimmed data directory')
    parser.add_argument('--output_dir', required=True, help='Path to the output directory')
   #parser.add_argument('--base_dir', required=True, help='Base directory for the pipeline')
    parser.add_argument('--spades_path', default='spades.py', help='Path to SPAdes executable')
    parser.add_argument('--quast_path', default='quast', help='Path to QUAST executable')
    parser.add_argument('--bbmap_path', default='bbmap.sh', help='Path to BBMap executable')
    parser.add_argument('--minlength', default=1000, type=int, help='Minimum length of contigs to keep')

    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Set the default reformat_path to the reformat.sh in the src directory
    reformat_default_path = os.path.join(script_dir, 'reformat.sh')
    parser.add_argument('--reformat_path', default=reformat_default_path, help='Path to reformat.sh script')

    args = parser.parse_args(args)

    setup_logging()

    # Define directories
    base_dir = args.base_dir
    assembly_dir = os.path.join(base_dir, "Assembly")
    trimmed_data_path = os.path.join(base_dir, 'Trim_data')
    all_contigs_dir = os.path.join(assembly_dir, "contigs", "all_contigs")
    filtered_contigs_dir = os.path.join(assembly_dir, "contigs", "filtered_contigs")
    quast_dir = os.path.join(assembly_dir, "Quast")
    pre_quast_dir = os.path.join(quast_dir, "pre_Quast")
    filtered_quast_dir = os.path.join(quast_dir, "filtered_Quast")
    multi_quast_dir = os.path.join(quast_dir, "multiQuast")
    coverage_dir = os.path.join(assembly_dir, "Coverage")
    multi_coverage_dir = os.path.join(coverage_dir, "multiCoverage")

    # Setup directories
    setup_directories([
        os.path.join(assembly_dir, "contigs", "samples"),
        all_contigs_dir,
        filtered_contigs_dir,
        pre_quast_dir,
        filtered_quast_dir,
        multi_quast_dir,
        coverage_dir,
        multi_coverage_dir
    ])

    # Run assembly
    run_assembly(trimmed_data_path, assembly_dir, spades_path=args.spades_path)

    # Reformat contigs
    reformat_contigs(all_contigs_dir, filtered_contigs_dir, minlength=args.minlength, reformat_path=args.reformat_path)

    # Run QUAST on all contigs
    run_quast(all_contigs_dir, pre_quast_dir, quast_path=args.quast_path)

    # Run QUAST on filtered contigs
    run_quast(filtered_contigs_dir, filtered_quast_dir, quast_path=args.quast_path)

    # Run BBMap to calculate coverage
    run_bbmap(trimmed_data_path, filtered_contigs_dir, coverage_dir, bbmap_path=args.bbmap_path)

if __name__ == "__main__":
    main()
