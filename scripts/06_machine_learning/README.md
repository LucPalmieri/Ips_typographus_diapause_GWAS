# Machine Learning Classification

Supervised machine learning pipeline to validate SNP predictive power and classify diapause phenotypes.

## Overview

This stage uses LightGBM with Bayesian hyperparameter optimization to:
1. Encode genotypes as integers
2. Optimize hyperparameters with Optuna
3. Train classifier on phenotyped individuals
4. Rank SNPs by feature importance
5. Validate model with mock populations

## Files

- `geno_converter_integer.py` - Encode VCF genotypes as integers
- `lightgbm_OPTUNA_improved.py` - Bayesian hyperparameter optimization
- `lightgbm_AUROC_kfold_improved.py` - Model training & cross-validation
- `phenotype_classifier_improved.py` - Apply to wild populations
- `generate_random.py` - Generate mock populations for validation

## Quick Start

### Step 1: Encode Genotypes

```bash
python geno_converter_integer.py \
  --input pruned_filtered75_Ips_ipyrad.recode.vcf \
  --output GENOTYPE_PHENOTYPED_INTEGER.csv
```

**Encoding**:
- `0` = Homozygous reference (0/0)
- `2` = Heterozygous (0/1)
- `3` = Homozygous alternate (1/1)
- `NA` = Missing (./.)

### Step 2: Hyperparameter Optimization

```bash
python lightgbm_OPTUNA_improved.py \
  --genotypes GENOTYPE_PHENOTYPED_INTEGER.csv \
  --phenotypes phenotypes.txt \
  --vca_scores VCA_results.scores \
  --n_trials 2000 \
  --cv_folds 50 \
  --output_dir optuna_results
```

**Parameters tuned**:
- `learning_rate`: 0.001 - 0.3
- `num_leaves`: 2 - 256
- `max_depth`: 1 - 15
- `min_child_weight`: 1 - 100
- `lambda_l1`, `lambda_l2`: Regularization

**Runtime**: 2-4 hours (2000 trials × 50-fold CV)

**Output**: `optuna_results/best_params.json`

### Step 3: Model Training

```bash
python lightgbm_AUROC_kfold_improved.py \
  --genotypes GENOTYPE_PHENOTYPED_INTEGER.csv \
  --phenotypes phenotypes.txt \
  --vca_scores VCA_results.scores \
  --best_params optuna_results/best_params.json \
  --output_dir lightgbm_results
```

**Performs**:
- Trains final model on all phenotyped samples
- Includes top 3 VCA eigenvectors as geographic covariates
- 50-fold cross-validation for performance estimation
- Calculates feature importance (gain-based)

**Outputs**:
- `SNP_importance_ranking.csv` - Ranked SNP predictive power
- `cross_validation_metrics.txt` - AUROC, F1, accuracy
- `model_predictions.csv` - Predicted probabilities

### Step 4: Model Validation

```bash
python generate_random.py \
  --genotypes GENOTYPE_PHENOTYPED_INTEGER.csv \
  --phenotypes phenotypes.txt \
  --output_dir mock_populations
```

**Validation**:
- Creates 5 synthetic populations with known phenotype ratios:
  - 0%, 25%, 50%, 75%, 100% facultative
- Predicts phenotypes for each mock population
- Compares predicted vs true frequencies
- Calculates accuracy metrics

## Cross-Validation Strategy

k-fold cross-validation with k=50:
- Dataset split into 50 equal folds
- Each fold used once as test set
- Prevents overfitting and provides robust estimates
- Critical for feature importance ranking

## Performance Metrics

| Metric | Description |
|--------|-------------|
| **AUROC** | Area under ROC curve (0-1, higher = better) |
| **Accuracy** | (TP + TN) / Total |
| **F1 Score** | Harmonic mean of precision & recall |
| **Sensitivity** | TP / (TP + FN) |
| **Specificity** | TN / (TN + FP) |
| **Precision** | TP / (TP + FP) |

## Feature Importance

LightGBM provides gain-based feature importance:
- **Gain**: Total reduction in loss from splits using that feature
- **Higher gain** = More important for prediction
- Includes regularization penalty (overfitted features penalized)

SNPs ranked by importance in `SNP_importance_ranking.csv`

## Input Files Required

| File | Format | Description |
|------|--------|-------------|
| `GENOTYPE_*.csv` | CSV | Integer-encoded genotypes (0, 2, 3, NA) |
| `phenotypes.txt` | TXT | Two-column: sample_id, phenotype (0 or 1) |
| `VCA_results.scores` | TXT | VCA scores (continuous covariates) |

## Outputs

| File | Description |
|------|-------------|
| `best_params.json` | Optimal hyperparameters |
| `SNP_importance_ranking.csv` | Ranked SNP feature importance |
| `cross_validation_metrics.txt` | Performance statistics |
| `model_predictions.csv` | Predicted probabilities for training data |
| `confusion_matrix.txt` | TP, FP, TN, FN counts |
| `optimization_plots.html` | Optuna optimization history (interactive) |

## Parameters

### LightGBM Model
- `boosting_type`: "gbdt" (gradient boosting decision tree)
- `objective`: "binary" (obligate vs facultative)
- `metric`: "auc" (Area Under Curve)
- `num_leaves`: 31-100 (optimized by Optuna)
- `max_depth`: 5-10 (optimized by Optuna)
- `learning_rate`: 0.05-0.1 (optimized by Optuna)
- `num_leaves`: Number of leaves in tree
- `lambda_l1`, `lambda_l2`: L1/L2 regularization

### Optuna Optimization
- `n_trials`: 2000 (number of parameter combinations)
- `cv_folds`: 50 (k-fold cross-validation)
- `objective_metric`: "auc" (maximize AUROC)
- `sampler`: "TPEsampler" (Tree Parzen Estimator)

## Missing Data Handling

LightGBM natively handles missing values:
- Missing values (NA) treated as separate category
- Each split decides which branch minimizes loss
- Histogram-based approach efficiently processes missingness
- **No imputation needed**

## Geographic Stratification Correction

Top 3 VCA eigenvectors included as covariates:
- Corrects for population structure
- Prevents spurious associations
- Improves SNP ranking reliability

## Notes

- **Reproducibility**: Set `random_state` parameter for consistent results
- **Class imbalance**: Use class weights if needed (script handles automatically)
- **Feature scaling**: Not needed for tree-based methods
- **Regularization**: Controlled via `lambda_l1` and `lambda_l2`

## Troubleshooting

**Memory Error During Optuna**
- Reduce `n_trials` (e.g., 1000 instead of 2000)
- Reduce `cv_folds` (e.g., 30 instead of 50)
- Run on higher memory node

**Model Accuracy Low**
- Check genotype encoding (0, 2, 3 for diploid)
- Verify phenotype assignment
- Ensure sufficient samples per class

**Optimization Not Converging**
- Increase `n_trials`
- Check data quality
- Verify VCA scores are reasonable

## References

Ke, G., et al. (2017). LightGBM: A highly efficient gradient boosting decision tree. NeurIPS.

Akiba, T., et al. (2019). Optuna: A next-generation hyperparameter optimization framework. KDD.

Varma, S., & Simon, R. (2006). Bias in error estimation when using cross-validation for model selection. BMC Bioinformatics.
