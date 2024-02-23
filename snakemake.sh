#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --output=out.snake
#SBATCH --error=err.snake
#SBATCH --account=pi-mengjiechen
#SBATCH --time=32:00:00
#SBATCH --partition=amd-hm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yliuchicago@uchicago.edu  # Where to send email
cd /home/yliuchicago/workspace/xinyuan/ubs0202
module load samtools
source activate snakemake
snakemake --core 12
