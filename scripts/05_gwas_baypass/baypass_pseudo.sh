#!/bin/bash

#SBATCH --partition=cpu-intel-01
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH --account=shame
#SBATCH --job-name=baypass
#SBATCH --output=slurm_out-%j.out
#SBATCH --error=slurm_err-%j.err


#load conda
module load baypass

## change working directory
cd /data/users/lpalmierirocha/barkbeetle/RADseq_ips/baypass

## simulate pseudo observations POD
g_baypass -gfile G.btapods -outprefix phenotyped_pseudo
