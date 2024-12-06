#!/bin/bash
#SBATCH --account=def-jronho
#SBATCH --time=5:00:00
#SBATCH --job-name=PlamidSpades
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=zhangbin.cai@mail.mcgill.ca
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

python HolmGenome.py --config config.yaml