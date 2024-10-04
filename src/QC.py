import argparse
import os
import subprocess
import logging
import glob

logging.basicConfig(filename='trim.log', filemode='a', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def setup_directories(base_dir):
    # Define the directories to be created
    dirs = [
        os.path.join(base_dir, 'Trim_data'),
        os.path.join(base_dir, 'QC', 'raw_data'),
        os.path.join(base_dir, 'QC', 'Trim')
    ]

    # Create each directory if it doesn't already exist
    for dir in dirs:
        os.makedirs(dir, exist_ok=True)

def run_fastqc(fastqc_path, data_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    fastq_files = glob.glob(os.path.join(data_dir, '*.fastq'))
    for fastq_file in fastq_files:
        cmd = f'{fastqc_path} -o {output_dir} {fastq_file}'
        subprocess.run(cmd, shell=True)
    
def run_multiqc(multiqc_path, input_dir):
    os.makedirs(input_dir, exist_ok=True)
    cmd = f'{multiqc_path} -o {input_dir} {input_dir}'
    subprocess.run(cmd, shell=True)
    
def run_trimmomatic(input_file1, input_file2, output_file1_paired, output_file1_unpaired, output_file2_paired, output_file2_unpaired, trimmomatic_path, adapters_path, illuminaclip='2:30:10', slidingwindow='4:15', leading='3', trailing='3', crop=None, headcrop=None, minlen='36'):
    cmd = f"java -jar {trimmomatic_path} PE -phred33 {input_file1} {input_file2} {output_file1_paired} {output_file1_unpaired} {output_file2_paired} {output_file2_unpaired} ILLUMINACLIP:{adapters_path}:{illuminaclip} LEADING:{leading} TRAILING:{trailing} SLIDINGWINDOW:{slidingwindow} MINLEN:{minlen} "
    if crop:
        cmd += f" CROP:{crop}"
    if headcrop:
        cmd += f" HEADCROP:{headcrop}"
    logging.info(f'Running command: {cmd}')
    with open('trimmomatic_output.log', 'a') as f:
        subprocess.run(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)

def run_bbduk(input_file1, input_file2, output_file1, output_file2, bbduk_path, adapters_path, ktrim='r', k='23', mink='11', hdist='1', tpe='t', tbo='t', qtrim='rl', trimq='10', minlen='25'):
    cmd = f"{bbduk_path} in1={input_file1} in2={input_file2} out1={output_file1} out2={output_file2} ref={adapters_path} ktrim={ktrim} k={k} mink={mink} hdist={hdist} tpe={tpe} tbo={tbo} qtrim={qtrim} trimq={trimq} minlen={minlen}"
    logging.info(f'Running command: {cmd}')
    with open('bbduk_output.log', 'a') as f:
        subprocess.run(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)

def main():
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
    parser.add_argument('--suffix1', default='_R1_001.fastq.gz', help='Suffix for the first file in each pair')
    parser.add_argument('--suffix2', default='_R2_001.fastq.gz', help='Suffix for the second file in each pair')
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
    parser.add_argument('--multiqc_path', default='multiqc', help='Path to the MultiQC executable')
    args = parser.parse_args()
    
    # Setup the directories
    setup_directories(args.output_dir)
    
    for root, dirs, files in os.walk(args.input_dir):
        files.sort()
        for i in range(0, len(files), 2):
            if files[i].endswith(args.suffix1) and files[i+1].endswith(args.suffix2):
                input_file1 = os.path.join(root, files[i])
                input_file2 = os.path.join(root, files[i+1])
                
               # Split the file name at the first period to separate the .fastq.gz extension
                base_name1 = files[i].split('.', 1)[0]
                base_name2 = files[i+1].split('.', 1)[0]
                
                if args.bbduk:
                    output_file1 = os.path.join(args.output_dir, f'{base_name1}.fastq.gz')
                    output_file2 = os.path.join(args.output_dir, f'{base_name2}.fastq.gz')
                    bbduk_path = args.bbduk_path if args.bbduk_path else 'bbduk.sh'
                    run_bbduk(input_file1, input_file2, output_file1, output_file2, bbduk_path, args.adapters_path, args.ktrim, args.k, args.mink, args.hdist, args.tpe, args.tbo, args.qtrim, args.trimq, args.minlen)
                else:
                    output_file1_paired = os.path.join(args.output_dir, f"{base_name1}_paired.fastq.gz")
                    output_file1_unpaired = os.path.join(args.output_dir, f"{base_name1}_unpaired.fastq.gz")
                    output_file2_paired = os.path.join(args.output_dir, f"{base_name2}_paired.fastq.gz")
                    output_file2_unpaired = os.path.join(args.output_dir, f"{base_name2}_unpaired.fastq.gz")
                    run_trimmomatic(input_file1, input_file2, output_file1_paired, output_file1_unpaired, output_file2_paired, output_file2_unpaired, args.trimmomatic_path, args.adapters_path,args.illuminaclip, args.slidingwindow, args.leading, args.trailing, args.crop, args.headcrop, args.minlen)

# Define paths to the raw and trimmed data
    raw_data_path = args.input_dir
    trimmed_data_path = os.path.join(args.output_dir, 'Trim_data')

    # Define paths to the FastQC and MultiQC executables
    fastqc_path = args.fastqc_path
    multiqc_path = args.multiqc_path

    # Run FastQC on raw data
    run_fastqc(fastqc_path, raw_data_path, os.path.join(args.output_dir, 'QC', 'raw_data'))
    
    # Run MultiQC on raw data FastQC results
    run_multiqc(multiqc_path, os.path.join(args.output_dir, 'QC', 'raw_data'))

    # Run FastQC on trimmed data
    run_fastqc(fastqc_path, trimmed_data_path, os.path.join(args.output_dir, 'QC', 'Trim'))
    
    # Run MultiQC on trimmed data FastQC results
    run_multiqc(multiqc_path, os.path.join(args.output_dir, 'QC', 'Trim'))


if __name__ == "__main__":
    main()