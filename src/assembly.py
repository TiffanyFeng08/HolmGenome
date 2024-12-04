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

def run_assembly(trimmed_data_path, assembly_dir, spades_path='spades.py'):
    # List sample directories in trimmed_data_path
    samples = [d for d in os.listdir(trimmed_data_path) if os.path.isdir(os.path.join(trimmed_data_path, d))]
    for sample in samples:
        sample_dir = os.path.join(trimmed_data_path, sample)
        # Find R1 and R2 files within the sample directory
        r1_files = glob.glob(os.path.join(sample_dir, '*_R1_paired.fastq*'))
        r2_files = glob.glob(os.path.join(sample_dir, '*_R2_paired.fastq*'))
        for r1_file in r1_files:
            # Find the corresponding R2 file
            r2_file = r1_file.replace('_R1_paired', '_R2_paired')
            if r2_file in r2_files:
                sample_output_dir = os.path.join(assembly_dir, "contigs", "samples", sample)
                os.makedirs(sample_output_dir, exist_ok=True)
                cmd = [
                    spades_path,
                    '--pe1-1', r1_file,
                    '--pe1-2', r2_file,
                    '-o', sample_output_dir
                ]
                run_subprocess(cmd, os.path.join(sample_output_dir, 'spades_output.log'))
                # Copy the resulting contigs to the all_contigs_dir
                contigs_src = os.path.join(sample_output_dir, "contigs.fasta")
                contigs_dest = os.path.join(assembly_dir, "contigs", "all_contigs", f"{sample}.fasta")
                if os.path.exists(contigs_src):
                    shutil.copy(contigs_src, contigs_dest)
                else:
                    logging.warning(f"Contigs file not found for sample {sample}")
            else:
                logging.warning(f"No matching R2 file for {r1_file} in sample {sample}")

def reformat_contigs(all_contigs_dir, filtered_contigs_dir, minlength=1000, reformat_path='reformat.sh'):
    contig_files = glob.glob(os.path.join(all_contigs_dir, '*.fasta'))
    for file in contig_files:
        sample_name = os.path.basename(file).split('.')[0]
        output_file = os.path.join(filtered_contigs_dir, f"{sample_name}.fasta")
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
    run_subprocess(cmd, os.path.join(output_dir, 'quast_output.log'))

def run_bbmap(trimmed_data_path, filtered_contigs_dir, coverage_dir, bbmap_path='bbmap.sh'):
    samples = [d for d in os.listdir(trimmed_data_path) if os.path.isdir(os.path.join(trimmed_data_path, d))]
    for sample in samples:
        sample_dir = os.path.join(trimmed_data_path, sample)
        r1_files = glob.glob(os.path.join(sample_dir, '*_R1_paired.fastq*'))
        r2_files = glob.glob(os.path.join(sample_dir, '*_R2_paired.fastq*'))
        for r1_file in r1_files:
            r2_file = r1_file.replace('_R1_paired', '_R2_paired')
            if r2_file in r2_files:
                ref_file = os.path.join(filtered_contigs_dir, f"{sample}.fasta")
                if not os.path.exists(ref_file):
                    logging.warning(f"Reference contigs not found for sample {sample}")
                    continue
                sample_coverage_dir = os.path.join(coverage_dir, sample)
                os.makedirs(sample_coverage_dir, exist_ok=True)
                cmd = [
                    bbmap_path,
                    f"in1={r1_file}",
                    f"in2={r2_file}",
                    f"ref={ref_file}",
                    f"covstats={os.path.join(sample_coverage_dir, sample + '_covstats.txt')}",
                    f"covhist={os.path.join(sample_coverage_dir, sample + '_covhist.tsv')}",
                    f"basecov={os.path.join(sample_coverage_dir, sample + '_basecov.txt')}",
                    "usejni=t",
                    f"out={os.path.join(sample_coverage_dir, sample + '_mapped.bam')}"
                ]
                run_subprocess(cmd, os.path.join(sample_coverage_dir, 'bbmap_output.log'))
            else:
                logging.warning(f"No matching R2 file for {r1_file} in sample {sample}")

def main(args=None):
    parser = argparse.ArgumentParser(description='Genome Assembly Pipeline')
    parser.add_argument('--output_dir', required=True, help='Path to the output directory')
    parser.add_argument('--trimmed_data_path', required=True, help='Path to the trimmed FASTQ files')
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

    base_dir = args.output_dir
    assembly_dir = os.path.join(base_dir, "Assembly")
    trimmed_data_path = args.trimmed_data_path
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
