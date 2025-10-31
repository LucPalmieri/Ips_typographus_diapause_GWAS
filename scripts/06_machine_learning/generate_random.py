#!/usr/bin/env python3

################################################################################
# Random Genotype Matrix Generator
# Purpose: Generate synthetic genotype data for testing machine learning models
# Author: Refactored version
# Description: Creates random SNP matrices with configurable dimensions and
#              genotype values for validation and testing purposes
################################################################################

import numpy as np
import pandas as pd
import argparse
import sys
from pathlib import Path


def generate_random_genotype_matrix(n_samples, 
                                    n_snps, 
                                    genotype_values=None,
                                    seed=None,
                                    snp_prefix="random",
                                    missing_rate=0.0):
    """
    Generate a random genotype matrix for testing purposes.
    
    Parameters
    ----------
    n_samples : int
        Number of samples (rows) to generate
    n_snps : int
        Number of SNP markers (columns) to generate
    genotype_values : list or None
        Possible genotype values. Default is ["0", "1", "2"] for diploid data
    seed : int or None
        Random seed for reproducibility
    snp_prefix : str
        Prefix for SNP column names (default: "random")
    missing_rate : float
        Proportion of missing data to introduce (0.0 to 1.0)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with random genotype data
    """
    
    # Set default genotype values if not provided
    if genotype_values is None:
        genotype_values = ["0", "1", "2"]
    
    # Validate inputs
    if n_samples <= 0 or n_snps <= 0:
        raise ValueError("Number of samples and SNPs must be positive integers")
    
    if not 0.0 <= missing_rate < 1.0:
        raise ValueError("Missing rate must be between 0.0 and 1.0 (exclusive)")
    
    # Set random seed if provided
    if seed is not None:
        np.random.seed(seed)
    
    print(f"Generating random genotype matrix...")
    print(f"  Samples: {n_samples}")
    print(f"  SNPs: {n_snps}")
    print(f"  Genotype values: {genotype_values}")
    print(f"  Missing rate: {missing_rate * 100:.1f}%")
    
    # Generate random matrix
    matrix = np.random.choice(genotype_values, size=(n_samples, n_snps))
    
    # Introduce missing data if specified
    if missing_rate > 0:
        n_missing = int(n_samples * n_snps * missing_rate)
        missing_indices = np.random.choice(
            n_samples * n_snps, 
            size=n_missing, 
            replace=False
        )
        flat_matrix = matrix.flatten()
        flat_matrix[missing_indices] = "NA"
        matrix = flat_matrix.reshape(n_samples, n_snps)
        print(f"  Introduced {n_missing} missing values")
    
    # Create column headers
    headers = [f"{snp_prefix}_{i:06d}" for i in range(1, n_snps + 1)]
    
    # Build DataFrame
    df = pd.DataFrame(matrix, columns=headers)
    
    print(f"Generated matrix shape: {df.shape}")
    
    return df


def add_sample_column(df, sample_prefix="Sample"):
    """
    Add a sample ID column to the genotype matrix.
    
    Parameters
    ----------
    df : pd.DataFrame
        Genotype matrix
    sample_prefix : str
        Prefix for sample names
    
    Returns
    -------
    pd.DataFrame
        DataFrame with sample column added
    """
    sample_ids = [f"{sample_prefix}_{i:04d}" for i in range(1, len(df) + 1)]
    df.insert(0, "sample", sample_ids)
    return df


def add_mock_phenotype(df, phenotype_proportion=0.5, phenotype_col="phenotype"):
    """
    Add a mock binary phenotype column to the matrix.
    
    Parameters
    ----------
    df : pd.DataFrame
        Genotype matrix
    phenotype_proportion : float
        Proportion of samples with phenotype value 1 (0.0 to 1.0)
    phenotype_col : str
        Name of the phenotype column
    
    Returns
    -------
    pd.DataFrame
        DataFrame with phenotype column added
    """
    n_samples = len(df)
    n_positive = int(n_samples * phenotype_proportion)
    
    # Create phenotype vector
    phenotypes = np.array([1] * n_positive + [0] * (n_samples - n_positive))
    np.random.shuffle(phenotypes)
    
    # Insert after sample column if it exists, otherwise at the beginning
    insert_pos = 1 if "sample" in df.columns else 0
    df.insert(insert_pos, phenotype_col, phenotypes)
    
    print(f"\nAdded mock phenotype:")
    print(f"  Class 1: {n_positive} samples ({phenotype_proportion * 100:.1f}%)")
    print(f"  Class 0: {n_samples - n_positive} samples ({(1 - phenotype_proportion) * 100:.1f}%)")
    
    return df


def save_matrix(df, output_file, show_preview=True):
    """
    Save the genotype matrix to a CSV file.
    
    Parameters
    ----------
    df : pd.DataFrame
        Genotype matrix to save
    output_file : str or Path
        Output file path
    show_preview : bool
        Whether to display a preview of the data
    """
    # Save to CSV
    df.to_csv(output_file, index=False)
    
    print(f"\nMatrix saved to: {output_file}")
    print(f"File size: {Path(output_file).stat().st_size / 1024:.1f} KB")
    
    if show_preview:
        print("\nPreview (first 5 rows, first 5 columns):")
        preview_cols = min(5, len(df.columns))
        print(df.iloc[:5, :preview_cols].to_string(index=False))
        if len(df.columns) > 5:
            print(f"... ({len(df.columns) - 5} more columns)")


def main():
    """Main execution function for command-line usage."""
    
    parser = argparse.ArgumentParser(
        description="Generate random genotype matrices for testing ML pipelines",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage: 100 samples, 500 SNPs
  %(prog)s -s 100 -n 500 -o random_matrix.csv
  
  # With sample IDs and phenotype
  %(prog)s -s 218 -n 1000 -o test_data.csv --add-samples --add-phenotype
  
  # Custom phenotype proportion (75%% positive)
  %(prog)s -s 200 -n 500 -o mock75.csv --add-samples --add-phenotype --pheno-prop 0.75
  
  # With missing data (10%% missing) and random seed for reproducibility
  %(prog)s -s 150 -n 800 -o data.csv --missing-rate 0.1 --seed 42
        """
    )
    
    # Required arguments
    parser.add_argument(
        "-s", "--samples",
        type=int,
        required=True,
        help="Number of samples (rows)"
    )
    
    parser.add_argument(
        "-n", "--snps",
        type=int,
        required=True,
        help="Number of SNP markers (columns)"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Output CSV file path"
    )
    
    # Optional arguments
    parser.add_argument(
        "--genotype-values",
        type=str,
        nargs="+",
        default=["0", "1", "2"],
        help="Possible genotype values (default: 0 1 2)"
    )
    
    parser.add_argument(
        "--snp-prefix",
        type=str,
        default="random",
        help="Prefix for SNP column names (default: random)"
    )
    
    parser.add_argument(
        "--missing-rate",
        type=float,
        default=0.0,
        help="Proportion of missing data (0.0 to 1.0, default: 0.0)"
    )
    
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility"
    )
    
    parser.add_argument(
        "--add-samples",
        action="store_true",
        help="Add sample ID column"
    )
    
    parser.add_argument(
        "--sample-prefix",
        type=str,
        default="Sample",
        help="Prefix for sample IDs (default: Sample)"
    )
    
    parser.add_argument(
        "--add-phenotype",
        action="store_true",
        help="Add mock binary phenotype column"
    )
    
    parser.add_argument(
        "--pheno-prop",
        type=float,
        default=0.5,
        help="Proportion of positive phenotypes (default: 0.5)"
    )
    
    parser.add_argument(
        "--no-preview",
        action="store_true",
        help="Don't show data preview"
    )
    
    args = parser.parse_args()
    
    # Print header
    print("\n" + "=" * 70)
    print("RANDOM GENOTYPE MATRIX GENERATOR")
    print("=" * 70 + "\n")
    
    try:
        # Generate matrix
        df = generate_random_genotype_matrix(
            n_samples=args.samples,
            n_snps=args.snps,
            genotype_values=args.genotype_values,
            seed=args.seed,
            snp_prefix=args.snp_prefix,
            missing_rate=args.missing_rate
        )
        
        # Add sample column if requested
        if args.add_samples:
            df = add_sample_column(df, sample_prefix=args.sample_prefix)
        
        # Add phenotype column if requested
        if args.add_phenotype:
            if not args.add_samples:
                print("\nWarning: --add-phenotype typically requires --add-samples")
            df = add_mock_phenotype(
                df, 
                phenotype_proportion=args.pheno_prop
            )
        
        # Save matrix
        save_matrix(df, args.output, show_preview=not args.no_preview)
        
        print("\n" + "=" * 70)
        print("Generation complete!")
        print("=" * 70 + "\n")
        
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
