# Execution Workflow Guide

Complete step-by-step instructions for running the diapause GWAS pipeline from raw reads through final phenotype predictions.

## Prerequisites

1. **Install software** (see `requirements/environment.yml` and `docs/SOFTWARE_VERSIONS.md`)
2. **Prepare input data**:
   - Raw fastq reads (demultiplexed or ready for demultiplexing)
   - Reference genome (FASTA format)
   - Population assignment file
   - Genome annotation (GFF3 format for Liftoff)

3. **Activate environment**:
   ```bash
   conda activate ips_typographus
   ```

---

## Stage 0: Raw Read Demultiplexing & Filtering

**Purpose**: Process raw sequencing reads and apply quality filters

**Input**: Raw fastq files, VCF from initial calling

**Output**: Filtered VCF file

### Step 0.1: Run Demultiplexing & VCF Filtering Script

```bash
cd scripts/00_demultiplex_filter/

# Usage: perl filter_script.pl input.vcf
perl mec14342-sup-0004-appendixs4.pl your_input.vcf

# This creates: your_input.filtered.vcf
# Applies filters:
# - Minimum quality (QUAL ≥ 21)
# - Minimum depth (DP ≥ 50)
# - Bi-allelic loci only
# - Minimum alternate allele frequency (≥5%)
# - Binomial filter for heterozygotes (p-value threshold: 0.05)
```

**Output**: `your_input.filtered.vcf`

---

## Stage 1: Assembly & Reference Mapping

**Purpose**: Assemble ddRAD sequences against reference genome

**Input**: Demultiplexed fastq files, reference genome FASTA

**Output**: VCF with called genotypes

### Step 1.1: Update ipyrad Parameters File

Edit `params-Ips_ipyrad.txt` with your file paths:

```
[0] assembly_name:              Ips_ipyrad
[1] project_dir:                /path/to/project
[2] raw_fastq_path:             (leave blank if using demultiplexed)
[3] barcodes_path:              (leave blank if using demultiplexed)
[4] sorted_fastq_path:          /path/to/demultiplexed/*.fastq
[5] assembly_method:            reference
[6] reference_sequence:         /path/to/ref_genome.fasta
[7] datatype:                   ddrad
[28] pop_assign_file:           /path/to/pop_assign.txt
```

### Step 1.2: Run ipyrad Assembly

```bash
cd scripts/01_assembly/

# Local run (if you have sufficient CPU/memory):
ipyrad -p params-Ips_ipyrad.txt -s 1234567 -c 4

# OR submit to HPC cluster:
sbatch ipyrad.sh

# Steps performed:
# 1. Demultiplexing (if needed)
# 2. Quality filtering & adapter trimming
# 3. Clustering & consensus calling
# 4. Joint genotyping
# 5. Output in multiple formats (including VCF)
```

**Output**: `Ips_ipyrad_outfiles/Ips_ipyrad.vcf`

---

## Stage 2: Variant Filtering & Linkage Pruning

**Purpose**: Remove low-quality variants and reduce linkage disequilibrium

**Input**: VCF from ipyrad

**Output**: Pruned SNP set in multiple formats

### Step 2.1: vcftools Filtering (Remove >25% Missing Data)

```bash
cd scripts/02_variant_filtering/

# Filter for missing data
vcftools --vcf Ips_ipyrad.vcf \
  --max-missing 0.75 \
  --recode \
  --out filtered75_Ips_ipyrad

# Output: filtered75_Ips_ipyrad.recode.vcf
```

### Step 2.2: PLINK LD Pruning

```bash
# Convert VCF to PLINK format and apply LD pruning
plink --vcf filtered75_Ips_ipyrad.recode.vcf \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.2 \
  --out pruned

# Prune variants
plink --vcf filtered75_Ips_ipyrad.recode.vcf \
  --allow-extra-chr \
  --extract pruned.prune.in \
  --make-bed \
  --out pruned_filtered75_Ips_ipyrad

# Final SNP set: 11,867 variants (from example log)
```

**Outputs**:
- `pruned.prune.in` - Variant list to keep
- `pruned_filtered75_Ips_ipyrad.bed/bim/fam` - PLINK binary format
- `pruned_filtered75_Ips_ipyrad.recode.vcf` - Filtered VCF

**Note**: Save PLINK VCF output for downstream analyses.

---

## Stage 3: Population Structure Analysis

**Purpose**: Characterize genetic structure across populations

**Input**: Pruned VCF file, population assignments

**Output**: VCA scores, STRUCTURE plots, population summaries

### Step 3.1: Variance Component Analysis (VCA)

```bash
cd scripts/03_population_genetics/

Rscript VCA_script_CLEANED.R \
  --vcf pruned_filtered75_Ips_ipyrad.recode.vcf \
  --pop pop_assign.txt \
  --output VCA_results

# Script performs:
# 1. Read VCF → genlight object
# 2. Impute missing genotypes (mean allele frequency)
# 3. Calculate genomic relationship matrix
# 4. Eigenvalue decomposition for VCA scores
# 5. Output VCA scores for use in GWAS

# Outputs: VCA_results.scores, VCA_results.kinship
```

**Key Output**: VCA scores (continuous covariate for GWAS)

### Step 3.2: STRUCTURE Analysis

```bash
# Requires pre-run STRUCTURE analysis (usually done separately)
# Visualize STRUCTURE output with pophelper

Rscript STRUCTURE_plot.R \
  --structure STRUCTURE_output.q \
  --k 3 \
  --output STRUCTURE_plot.pdf

# Produces admixture proportion bar plots for K=3 populations
```

**Outputs**:
- VCA scores file (for GWAS)
- STRUCTURE visualization PDF

---

## Stage 4: Genome Annotation & SNP-to-Gene Mapping

**Purpose**: Transfer annotations and map SNPs to genomic features

**Input**: Old genome annotation (GFF3), new reference genome (FASTA)

**Output**: Feature annotations, SNP-to-gene mapping

### Step 4.1: Liftoff Annotation Transfer

```bash
cd scripts/04_annotation/

# Transfer annotations from old to new genome
liftoff -g old_annotation.gff3 \
  -o new_annotation.gff3 \
  -a agp_file.agp \
  new_genome.fasta \
  old_genome.fasta

# Output: new_annotation.gff3 (transferred features)
```

**Key Output**: `new_annotation.gff3`

### Step 4.2: SNP-to-Gene Feature Mapping

```bash
# Run R script to find SNP-feature overlaps
Rscript Gene_feature_sorting.rmd \
  --snps GWAS_moderateBF.txt \      # SNPs to map (from GWAS step)
  --features liftoff_features.txt \ # Liftoff output
  --output SNP_feature_overlaps.tsv

# Performs GenomicRanges overlap detection
# Output: SNP_feature_overlaps_GWAS_moderateBF.tsv
```

**Outputs**:
- `SNP_feature_overlaps_GWAS_moderateBF.tsv` - SNP-feature mapping

---

## Stage 5: Genome-Wide Association Study (GWAS)

**Purpose**: Identify SNPs associated with diapause phenotype via Bayesian inference

**Input**: Pruned VCF (phenotyped individuals only), VCA scores, phenotype data

**Output**: Bayes Factors, SNP associations, ranked candidate loci

### Step 5.1: Format Data for BayPass

```bash
cd scripts/05_gwas_baypass/

# Convert pruned VCF to BayPass allele count format
python reshaper_baypass_improved.py \
  pruned_filtered75_Ips_ipyrad.recode.vcf \
  phenotyped.geno

# Output: phenotyped.geno (BayPass format)
```

**Key Output**: `phenotyped.geno`

### Step 5.2: Run BayPass Core Model

```bash
# Initialize Omega matrix (population covariance) under neutrality
# Submit to cluster:
sbatch baypass_core.sh

# Runs three independent MCMC chains:
# g_baypass -gfile phenotyped.geno \
#   -seed [SEED] -npilot 20 -burnin 2500 \
#   -outprefix phenotyped_core_run[1-5]

# Sets parameters:
# - 20 pilot phases for adaptation
# - 2,500 burn-in iterations
# - Different random seeds for replication
```

**Outputs** (per run):
- `phenotyped_core_run[1-5]_summary_pi_xtx.txt` - XtX matrix
- `phenotyped_core_run[1-5]_summary_fst.txt` - Fst values
- `phenotyped_core_run[1-5]_pi_xtx_NOGT_ANC.txt` - Omega matrix

### Step 5.3: Estimate Population Covariance Matrix

```bash
# Process Omega matrix output from core model
# Use run1 Omega as fixed parameter for downstream analyses

# From output: Extract mean Omega across runs
# Save as: contrast_VCA.txt (or similar for covariate model)
```

### Step 5.4: BayPass Importance Sampling with VCA Covariate

```bash
# Submit to cluster:
sbatch baypass_is.sh

# Runs importance sampling with VCA scores as covariate:
# g_baypass -gfile phenotyped.geno \
#   -efile contrast_PCA.txt \          # VCA scores as covariate
#   -seed 1234 -npilot 15 -pilotlength 500 \
#   -burnin 2500 -outprefix phenotyped_importance_sampling_cov_run1

# Tests for association between SNP allele frequency and VCA score
```

**Key Outputs**:
- `phenotyped_importance_sampling_cov_run1_summary_betai_fst.txt` - Bayes Factors
- `phenotyped_importance_sampling_cov_run2_summary_betai_fst.txt` - Bayes Factors
- `phenotyped_importance_sampling_cov_run3_summary_betai_fst.txt` - Bayes Factors
- `phenotyped_importance_sampling_cov_run4_summary_betai_fst.txt` - Bayes Factors
- `phenotyped_importance_sampling_cov_run5_summary_betai_fst.txt` - Bayes Factors

### Step 5.5: Extract & Rank SNP Associations

```bash
# Process BayPass output with R script
Rscript baypass_core_model_script.R \
  --bayes_factor phenotyped_importance_sampling_cov_run1_summary_betai_fst.txt \
  --min_bf 6 \
  --output gwas_results

# Extracts SNPs with Bayes Factor > 6 (moderate association)
# Classification: BF > 2 (weak), > 6 (moderate), > 10 (strong)

# Outputs: 
# - GWAS_moderateBF.txt (SNPs with BF > 6)
# - gwas_summary.txt (statistics)
```

**Key Outputs**:
- `GWAS_moderateBF.txt` - Moderate-to-strong associations (for GO enrichment)
- `gwas_summary.txt` - GWAS statistics

### Step 5.6: BayPass Pseudo-Observations for Model Validation

```bash
# Run pseudo-observation test (optional, for model validation)
sbatch baypass_pseudo.sh

# g_baypass -gfile G.btapods -outprefix phenotyped_pseudo
```

---

## Stage 6: Machine Learning Classification

**Purpose**: Independently validate SNP predictive power using supervised learning

**Input**: Pruned genotypes (phenotyped individuals), phenotype labels, VCA scores

**Output**: Ranked SNP importance, trained classifier, cross-validation metrics

### Step 6.1: Encode Genotypes as Integers

```bash
cd scripts/06_machine_learning/

python geno_converter_integer.py \
  --input pruned_filtered75_Ips_ipyrad.recode.vcf \
  --output GENOTYPE_PHENOTYPED_INTEGER.csv

# Encoding:
# 0 = homozygous reference (0/0)
# 2 = heterozygous (0/1 or 1/0)
# 3 = homozygous alternate (1/1)
# NA = missing (./.)

# LightGBM handles missing values natively via histogram-based approach
```

**Output**: `GENOTYPE_PHENOTYPED_INTEGER.csv`

### Step 6.2: Hyperparameter Optimization with Optuna

```bash
python lightgbm_OPTUNA_improved.py \
  --genotypes GENOTYPE_PHENOTYPED_INTEGER.csv \
  --phenotypes phenotypes.txt \           # Obligate (0) vs Facultative (1)
  --vca_scores VCA_results.scores \
  --n_trials 2000 \
  --cv_folds 50 \
  --output_dir optuna_results

# Performs:
# 1. Bayesian optimization over parameter space:
#    - learning_rate
#    - num_leaves
#    - max_depth
#    - min_child_weight
#    - lambda_l1, lambda_l2 (regularization)
# 2. 50-fold cross-validation per trial
# 3. AUROC as optimization metric
# 4. Saves best hyperparameters
```

**Outputs**:
- `optuna_results/best_params.json` - Optimal hyperparameters
- `optuna_results/optuna_study.pkl` - Full optimization history
- `optuna_results/optimization_plots.html` - Visualization

### Step 6.3: Model Training & Feature Importance Ranking

```bash
python lightgbm_AUROC_kfold_improved.py \
  --genotypes GENOTYPE_PHENOTYPED_INTEGER.csv \
  --phenotypes phenotypes.txt \
  --vca_scores VCA_results.scores \
  --best_params optuna_results/best_params.json \
  --output_dir lightgbm_results

# Performs:
# 1. Train LightGBM on complete phenotyped dataset
# 2. Include top 3 VCA eigenvectors as covariates (geographic correction)
# 3. Evaluate via k-fold cross-validation (k=50)
# 4. Calculate feature importance (gain-based ranking)
# 5. Generate performance metrics (AUROC, accuracy, F1)

# Outputs:
# - SNP_importance_ranking.csv (all SNPs ranked)
# - cross_validation_metrics.txt (AUROC, F1, accuracy, precision, recall)
# - model_predictions.csv (predicted probabilities for all samples)
```

**Key Outputs**:
- `SNP_importance_ranking.csv` - Ranked SNP predictive power
- Cross-validation metrics file

### Step 6.4: Model Validation with Mock Populations

```bash
python generate_random.py \
  --genotypes GENOTYPE_PHENOTYPED_INTEGER.csv \
  --phenotypes phenotypes.txt \
  --output_dir mock_populations

# Generates synthetic populations with known phenotype ratios:
# - 0% facultative (100% obligate)
# - 25% facultative / 75% obligate
# - 50% facultative / 50% obligate
# - 75% facultative / 25% obligate
# - 100% facultative (0% obligate)

# For each mock population:
# 1. Predict phenotypes using trained model
# 2. Compare predicted vs true ratio
# 3. Calculate F1 score, accuracy, AUROC
# 4. Validate model reliability

# Outputs: Confusion matrices, performance metrics per mock population
```

**Key Output**: Model validation metrics on synthetic data

---

## Stage 7: Functional Annotation & GO Enrichment

**Purpose**: Interpret biology of candidate SNPs via gene function analysis

**Input**: SNP associations (GWAS & ML), gene annotations, GO term mappings

**Output**: Enriched biological pathways, candidate genes, functional categories

### Step 7.1: Prepare Input Files for GO Enrichment

```bash
cd scripts/07_functional_annotation/

# Create three input files:

# 1. Gene lengths (all genes in RAD-seq dataset)
# GeneLengths_RAD_universe.txt:
# gene_id    length
# gene_001   1500
# gene_002   2300
# ...

# 2. Gene-to-GO term mappings
# Gene2GO.txt:
# gene_id    GO_term
# gene_001   GO:0006629
# gene_001   GO:0008150
# gene_002   GO:0003677
# ...

# 3. SNP-to-gene mapping (from Stage 4)
# SNPs_with_geneID_all.tsv (or use foreground_mod_BF6.txt for GWAS SNPs)
# Contains: CHROM, POS, snpID, gene_id
```

### Step 7.2: Run GOseq Enrichment Analysis

```bash
# Run enrichment on GWAS moderate-BF SNPs (BF > 6)
Rscript enrichment_analysis_good_script.rmd

# Script performs:
# 1. Read foreground genes (hits from GWAS/ML analysis)
# 2. Build background universe (all genes in dataset)
# 3. Account for gene length bias using Wallenius model
# 4. Test for GO term over-representation
# 5. Apply FDR correction

# Configuration within script:
# - Background: all genes in GeneLengths file
# - Foreground: genes harboring SNPs with BF > 6
# - Method: Wallenius non-central hypergeometric
# - Ontologies: BP (Biological Process), MF (Molecular Function), CC (Cellular Component)

# Outputs:
# - GOseq_all_terms_modBF6.tsv (all tested terms)
# - GOseq_significant_terms_modBF6.tsv (FDR < 0.05)
```

**Key Outputs**:
- `GOseq_significant_terms_modBF6.tsv` - Enriched GO categories

---

## Stage 8: Wild Population Classification & Phenotype Prediction

**Purpose**: Apply trained classifier to wild individuals and estimate population-level phenotype frequencies

**Input**: Trained LightGBM model, genotypes from wild individuals, VCA scores

**Output**: Predicted phenotypes, population-level frequency estimates, outbreak risk assessment

### Step 8.1: Prepare Wild Individual Genotypes

```bash
# From pruned VCF, extract wild individuals only
# Use same encoding as phenotyped samples:

python geno_converter_integer.py \
  --input pruned_filtered75_Ips_ipyrad.recode.vcf \
  --samples wild_individuals.txt \    # List of wild sample IDs
  --output GENOTYPE_WILD_INTEGER.csv
```

### Step 8.2: Predict Phenotypes in Wild Populations

```bash
python phenotype_classifier_improved.py \
  --trained_model lightgbm_results/trained_model.pkl \
  --wild_genotypes GENOTYPE_WILD_INTEGER.csv \
  --wild_vca_scores VCA_wild.scores \
  --threshold 0.5 \                    # Default classification threshold
  --output_dir wild_predictions

# Performs:
# 1. Generate probability predictions (0-1 scale)
# 2. Convert to binary phenotypes (threshold = 0.5)
# 3. Calculate per-population frequencies
# 4. Generate confusion matrices (vs field phenotypes if available)

# Outputs:
# - wild_predictions.csv (individual predictions + probabilities)
# - population_frequencies.txt (% facultative per population)
# - outbreak_risk_assessment.txt
```

**Key Outputs**:
- `wild_predictions.csv` - Predicted phenotype for each wild individual
- `population_frequencies.txt` - Estimated % obligate/facultative per population

### Step 8.3: Generate Statistical Summary

```bash
Rscript scripts/utils/binomial_tests.R \
  --observed population_frequencies.txt \
  --populations LOW HIGH NORTH

# Performs binomial tests:
# H0: Expected ratio (e.g., 50% obligate)
# Compares observed vs expected proportions
# Two-sided tests with 95% confidence intervals
```

**Output**: Statistical validation of phenotype predictions

---

## Stage 9: Results Visualization

**Purpose**: Generate publication-quality figures

**Input**: GWAS results, VCA scores, ML predictions, annotations

**Output**: Publication-ready plots

### Step 9.1: Manhattan Plot (BF Values)

```bash
cd scripts/08_visualization/

Rscript BF_manhattan_plot.R \
  --gwas_results gwas_summary.txt \
  --output manhattan_plot.pdf

# Shows Bayes Factor significance across genome
# Highlights SNPs with BF > 6 (moderate) and BF > 10 (strong)
```

### Step 9.2: Covariance Structure Heatmap

```bash
Rscript Bhru.R \
  --omega_matrix phenotyped_core_run1_summary_pi_xtx.txt \
  --output covariance_heatmap.pdf

# Shows population differentiation patterns
# Hierarchical clustering of populations
```

### Step 9.3: STRUCTURE Bar Plot

```bash
Rscript STRUCTURE_plot.R \
  --structure_output STRUCTURE_run_1.q \
  --k 3 \
  --output structure_barplot.pdf

# Individual admixture proportions for K=3
```

---

## Data & File Management

### Essential Intermediate Files to Retain

```
# Variant Calling & Filtering
✓ pruned_filtered75_Ips_ipyrad.recode.vcf      # Pruned SNP set

# Population Structure
✓ VCA_results.scores                            # VCA covariates
✓ STRUCTURE_run_*.q                             # STRUCTURE output

# GWAS
✓ phenotyped.geno                               # BayPass input
✓ phenotyped_importance_sampling_cov_run1_*    # BayPass results

# Machine Learning
✓ GENOTYPE_PHENOTYPED_INTEGER.csv               # Encoded genotypes
✓ SNP_importance_ranking.csv                    # ML feature ranking
✓ optuna_results/best_params.json               # Optimal hyperparameters

# Annotations
✓ new_annotation.gff3                           # Liftoff output
✓ SNP_feature_overlaps_GWAS_moderateBF.tsv     # SNP-gene mapping

# Results
✓ GOseq_significant_terms_modBF6.tsv            # Enriched GO terms
✓ wild_predictions.csv                          # Wild phenotypes
```

## Next Steps After Analysis

1. **Validate findings**: Compare GWAS and ML SNP rankings
2. **Literature review**: Search for candidate genes in diapause literature
3. **Functional studies**: Design targeted experiments on top candidates
4. **Population monitoring**: Apply predictor to manage outbreak risk
5. **Manuscript preparation**: Generate figures from Stage 9 outputs

---
