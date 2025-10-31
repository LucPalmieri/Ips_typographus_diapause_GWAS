#!/bin/bash

#SBATCH --partition=cpu-intel-01
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G
#SBATCH --account=shame
#SBATCH --job-name=ipyrad
#SBATCH --output=slurm_out-%j.out
#SBATCH --error=slurm_err-%j.err


#load conda
#module load ipyrad

## change into the directory where your params file resides
cd /data/users/lpalmierirocha/barkbeetle/RADseq_ips

## call ipyrad on your params file
ipyrad -p params-Ips_ipyrad.txt -s 1234567 -c ${SLURM_CPUS_PER_TASK}

