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

## call baypass under core model ### baypass -gfile geno.bta14 -outprefix anacore ###
g_baypass -gfile phenotyped.geno -seed 1234 -npilot 20 -burnin 2500 -nthreads 4 \
 -outprefix phenotyped_core_run1 ;

g_baypass -gfile phenotyped.geno -seed 4321 -npilot 20 -burnin 2500 -nthreads 4 \
 -outprefix phenotyped_core_run2 ;

 g_baypass -gfile phenotyped.geno -seed 5678 -npilot 20 -burnin 2500 -nthreads 4 \
 -outprefix phenotyped_core_run3 ;
