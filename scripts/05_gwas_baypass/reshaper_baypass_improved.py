#!/usr/bin/env python3

################################################################################
# VCF to BayPass Format Converter
# Source: BayPass repository (adapted)
# Purpose: Convert VCF files to BayPass input format for population genetics analysis
# 
# BayPass is a population genomics software for detecting SNPs under selection
# and/or associated with population-specific covariates.
# Reference: Gautier (2015) Genetics, doi:10.1534/genetics.115.181453
################################################################################

"""
Convert VCF file to BayPass input format

Usage:
    python reshaper_baypass.py input.vcf popmap.txt output.baypass

Input files:
    - input.vcf: VCF file with SNP genotypes
    - popmap.txt: Two-column tab-delimited file mapping samples to populations
                  Format: SAMPLE_ID    POPULATION
    
Output:
    - BayPass format file with allele counts per population per SNP
    - Each row = one SNP
    - Columns = alternating ref/alt allele counts for each population
"""

import sys
import time
from pathlib import Path


################################################################################
# FUNCTIONS
################################################################################

def update_progress(job_title, progress):
    """
    Display a simple progress bar in the terminal.
    
    Parameters
    ----------
    job_title : str
        Description of the current task
    progress : float
        Progress value between 0.0 and 1.0
    """
    length = 20  # Progress bar length
    block = int(round(length * progress))
    msg = "\r{0}: [{1}] {2}%".format(
        job_title,
        "#" * block + "-" * (length - block),
        round(progress * 100, 0)
    )
    if progress >= 1:
        msg += " DONE\r\n"
    sys.stdout.write(msg)
    sys.stdout.flush()


def parse_arguments():
    """
    Parse command-line arguments.
    
    Returns
    -------
    tuple
        (input_vcf, popmap_file, output_file)
    """
    try:
        input_vcf = sys.argv[1]
        popmap_file = sys.argv[2]
        output_file = sys.argv[3]
        
        # Validate input files exist
        if not Path(input_vcf).exists():
            raise FileNotFoundError(f"VCF file not found: {input_vcf}")
        if not Path(popmap_file).exists():
            raise FileNotFoundError(f"Population map file not found: {popmap_file}")
        
        return input_vcf, popmap_file, output_file
        
    except IndexError:
        print(__doc__)
        print("\nError: Missing required arguments")
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"\nError: {e}")
        sys.exit(1)


def load_population_map(popmap_file):
    """
    Read population map file and create a dictionary of populations.
    
    Parameters
    ----------
    popmap_file : str
        Path to population map file
    
    Returns
    -------
    dict
        Dictionary mapping population names to lists of sample IDs
    """
    pop_dict = {}
    
    with open(popmap_file) as popfile:
        for line in popfile:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
                
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"Warning: Skipping malformed line: {line}")
                continue
            
            sample_id, population = parts[:2]
            
            # Add sample to population
            if population in pop_dict:
                pop_dict[population].append(sample_id)
            else:
                pop_dict[population] = [sample_id]
    
    return pop_dict


def count_vcf_lines(vcf_file):
    """
    Count total lines in VCF file for progress tracking.
    
    Parameters
    ----------
    vcf_file : str
        Path to VCF file
    
    Returns
    -------
    int
        Number of lines in file
    """
    with open(vcf_file) as f:
        return sum(1 for _ in f)


def parse_genotype(gt_string):
    """
    Parse genotype string and return alleles.
    
    Parameters
    ----------
    gt_string : str
        Genotype string (e.g., "0/1" or "0/0:20:30")
    
    Returns
    -------
    tuple or None
        (allele1, allele2) as integers, or None if missing
    """
    gt = gt_string.split(":")[0]
    
    if gt == "./.":
        return None
    
    if "/" in gt:
        alleles = gt.split("/")
    elif "|" in gt:
        alleles = gt.split("|")
    else:
        return None
    
    try:
        return (int(alleles[0]), int(alleles[1]))
    except (ValueError, IndexError):
        return None


def process_vcf(input_vcf, pop_dict):
    """
    Process VCF file and extract allele counts per population.
    
    Parameters
    ----------
    input_vcf : str
        Path to input VCF file
    pop_dict : dict
        Dictionary of populations and their samples
    
    Returns
    -------
    dict
        Dictionary of allele counts per population per SNP
    int
        Total number of SNPs processed
    """
    # Get file statistics
    num_lines = count_vcf_lines(input_vcf)
    n_populations = len(pop_dict)
    pop_names = list(pop_dict.keys())
    samples_per_pop = [len(pop_dict[pop]) for pop in pop_names]
    
    # Storage for genotype data
    geno_dict = {}
    snp_count = 0
    i_progress = 0
    
    with open(input_vcf) as vcffile:
        for line in vcffile:
            i_progress += 1
            
            # Update progress bar
            if i_progress % 100 == 0 or i_progress == num_lines:
                update_progress("Parsing VCF", i_progress / num_lines)
            
            # Skip metadata lines
            if line.startswith("##"):
                continue
            
            # Process header line
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                n_samples_vcf = len(fields[9:])
                n_samples_popmap = sum(samples_per_pop)
                
                if n_samples_vcf != n_samples_popmap:
                    print(f"\nError: Sample count mismatch!")
                    print(f"  VCF file: {n_samples_vcf} samples")
                    print(f"  Population map: {n_samples_popmap} samples")
                    sys.exit(1)
                continue
            
            # Process SNP line
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue
            
            scaffold, position = fields[0], fields[1]
            genotypes = fields[9:]
            
            snp_count += 1
            
            # Process each population
            sample_idx = 0
            for pop_idx, pop_name in enumerate(pop_names):
                n_samples_in_pop = samples_per_pop[pop_idx]
                
                # Extract genotypes for this population
                pop_genotypes = genotypes[sample_idx:sample_idx + n_samples_in_pop]
                
                # Count alleles
                allele_counts = [0, 0]  # [ref_count, alt_count]
                
                for gt_string in pop_genotypes:
                    gt = parse_genotype(gt_string)
                    if gt is not None:
                        allele_counts[0] += (2 - gt[0] - gt[1])  # Reference alleles
                        allele_counts[1] += (gt[0] + gt[1])      # Alternate alleles
                
                # Store counts
                if pop_name in geno_dict:
                    geno_dict[pop_name].append(allele_counts)
                else:
                    geno_dict[pop_name] = [allele_counts]
                
                # Move to next population
                sample_idx += n_samples_in_pop
    
    return geno_dict, snp_count


def write_baypass_output(geno_dict, snp_count, output_file):
    """
    Write allele counts to BayPass format file.
    
    Parameters
    ----------
    geno_dict : dict
        Dictionary of allele counts per population
    snp_count : int
        Total number of SNPs
    output_file : str
        Path to output file
    """
    with open(output_file, "w") as outfile:
        for snp_idx in range(snp_count):
            # Get counts for all populations for this SNP
            counts_all_pops = []
            for pop_name in geno_dict.keys():
                counts = geno_dict[pop_name][snp_idx]
                counts_all_pops.extend(counts)
            
            # Write as tab-delimited line
            outfile.write("\t".join(map(str, counts_all_pops)) + "\n")


def print_summary(pop_dict, snp_count):
    """
    Print summary statistics.
    
    Parameters
    ----------
    pop_dict : dict
        Dictionary of populations
    snp_count : int
        Number of SNPs processed
    """
    print("\n" + "=" * 50)
    print("CONVERSION SUMMARY")
    print("=" * 50)
    print(f"Number of populations:  {len(pop_dict)}")
    print(f"Total SNPs processed:   {snp_count}")
    print("\nSamples per population:")
    for pop_name, samples in sorted(pop_dict.items()):
        print(f"  {pop_name:20s}: {len(samples)} samples")
    print("=" * 50)


################################################################################
# MAIN
################################################################################

def main():
    """Main execution function."""
    
    print("\n" + "=" * 70)
    print("VCF TO BAYPASS FORMAT CONVERTER")
    print("=" * 70 + "\n")
    
    # Parse arguments
    input_vcf, popmap_file, output_file = parse_arguments()
    
    print(f"Input VCF:        {input_vcf}")
    print(f"Population map:   {popmap_file}")
    print(f"Output file:      {output_file}\n")
    
    # Load population map
    print("Loading population map...")
    pop_dict = load_population_map(popmap_file)
    
    if not pop_dict:
        print("Error: No populations found in population map file")
        sys.exit(1)
    
    print(f"Found {len(pop_dict)} populations\n")
    
    # Process VCF file
    print("Processing VCF file...")
    geno_dict, snp_count = process_vcf(input_vcf, pop_dict)
    
    # Write output
    print("\nWriting BayPass output file...")
    write_baypass_output(geno_dict, snp_count, output_file)
    
    # Print summary
    print_summary(pop_dict, snp_count)
    
    print(f"\nOutput saved to: {output_file}")
    print("\nConversion complete!\n")


if __name__ == "__main__":
    main()
