# GWAS Analysis with BayPass

Bayesian genome-wide association study to identify SNPs associated with diapause phenotype.

## Overview

This stage performs a Bayesian inference approach using BayPass v2.41 to:
1. Estimate population covariance (Omega) matrix under neutrality
2. Test for associations between SNP allele frequencies and diapause phenotype
3. Account for population structure using VCA scores as a continuous covariate
4. Generate Bayes Factors for all SNPs

## Files

- `reshaper_baypass_improved.py` - Convert VCF to BayPass allele count format
- `baypass_core_model_script.R` - Extract and process Bayes Factor results
- `baypass_core.sh` - Run BayPass core model (Omega matrix estimation)
- `baypass_is.sh` - Run importance sampling with VCA covariate
- `baypass_pseudo.sh` - Generate pseudo-observations for validation

## Quick Start

### Step 1: Format Data

```bash
python reshaper_baypass_improved.py phenotyped.vcf phenotyped.geno
```

Input: `phenotyped.vcf` (VCF with phenotyped individuals only)  
Output: `phenotyped.geno` (BayPass format)

### Step 2: Run Core Model

```bash
sbatch baypass_core.sh
```

Runs 3 independent MCMC chains to estimate Omega matrix:
- Seed 1234: `phenotyped_core_run1_*`
- Seed 4321: `phenotyped_core_run2_*`  
- Seed 5678: `phenotyped_core_run3_*`

Key output: Omega matrix (population covariance)

### Step 3: Importance Sampling with Covariate

```bash
sbatch baypass_is.sh
```

Tests association with VCA scores (continuous covariate):
- Input: `contrast_PCA.txt` (VCA scores)
- Output: `phenotyped_importance_sampling_cov_run1_summary_betai_fst.txt`

**Key metric**: Bayes Factors for each SNP

### Step 4: Extract Results

```bash
Rscript baypass_core_model_script.R
```

Extracts SNPs with BF > 6 (moderate association):
- Output: `GWAS_moderateBF.txt` (for GO enrichment)
- Output: `gwas_summary.txt` (statistics)

## Parameters

### Core Model (baypass_core.sh)
- `npilot`: 20 (pilot phases for adaptation)
- `burnin`: 2500 (burn-in iterations)
- `nthreads`: 4 (CPU threads)

### Importance Sampling (baypass_is.sh)
- `npilot`: 15 (reduced for speed)
- `pilotlength`: 500
- `burnin`: 2500

## Bayes Factor Interpretation

Classification based on Jeffreys (1961):
- **BF > 2**: Weak evidence
- **BF > 6**: Moderate evidence
- **BF > 10**: Strong evidence
- **BF > 100**: Very strong evidence

## Outputs

| File | Description |
|------|-------------|
| `phenotyped_core_run*_summary_pi_xtx.txt` | Omega matrix |
| `phenotyped_core_run*_summary_fst.txt` | Fst values |
| `phenotyped_importance_sampling_cov_run1_summary_betai_fst.txt` | Bayes Factors |
| `GWAS_moderateBF.txt` | SNPs with BF > 6 |
| `gwas_summary.txt` | Results summary |

## Notes

- BayPass operates on **population-level allele counts** (not individual genotypes)
- Covariate model uses **mean VCA score per phenotype group**
- Multiple runs allow assessment of parameter stability
- Recommend using Omega from core_run1 as fixed parameter

## References

Gautier, M. (2015). Genome-wide scan for adaptive divergence and association with population-specific covariates. Genetics, 201(4), 1555-1575.

Jeffreys, H. (1961). Theory of Probability. Oxford University Press.

## Troubleshooting

**Error: Memory limit exceeded**
- Reduce `nthreads` or submit to larger memory queue

**Bayes Factors all near zero**
- Check VCA scores are normalized correctly
- Verify phenotype assignment is correct

**Long runtime**
- Increase `nthreads` if available
- Reduce `pilotlength` (though this may reduce accuracy)
