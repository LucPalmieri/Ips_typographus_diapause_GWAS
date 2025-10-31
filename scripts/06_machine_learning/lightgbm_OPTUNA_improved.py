#!/usr/bin/env python3
"""
LightGBM Hyperparameter Optimization using Optuna

This script performs hyperparameter tuning for a LightGBM regression model using Optuna.
The model is trained on residualized continuous phenotypes and evaluated using binary
classification metrics after applying a custom threshold.

Features:
- K-fold cross-validation for robust evaluation
- Optuna-based Bayesian optimization
- Early stopping to prevent overfitting
- Comprehensive logging and result visualization
- Reproducible with seed control

Author: [Your Name]
Date: 2025-10-30
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import lightgbm as lgb
import optuna
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score, accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns

# Set style for plots
sns.set_style("whitegrid")


def setup_logging(log_file: str = None, verbose: bool = True) -> logging.Logger:
    """
    Configure logging for the script.
    
    Args:
        log_file: Path to log file. If None, logs only to console.
        verbose: If True, set logging level to INFO, else WARNING.
    
    Returns:
        Configured logger object.
    """
    log_level = logging.INFO if verbose else logging.WARNING
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode='w'))
    
    logging.basicConfig(
        level=log_level,
        format=log_format,
        handlers=handlers
    )
    
    return logging.getLogger(__name__)


def load_and_prepare_data(
    data_file: str,
    sample_col: str = 'sample',
    binary_col: str = 'phenotype_binary',
    resid_col: str = 'phenotype_resid'
) -> Tuple[pd.DataFrame, pd.Series, pd.Series, List[str]]:
    """
    Load and prepare data for model training.
    
    Args:
        data_file: Path to CSV file containing SNP data and phenotypes.
        sample_col: Name of sample ID column.
        binary_col: Name of binary phenotype column.
        resid_col: Name of residualized continuous phenotype column.
    
    Returns:
        Tuple of (X features, y_resid continuous, y_binary labels, SNP column names).
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Loading data from: {data_file}")
    data = pd.read_csv(data_file)
    logger.info(f"Loaded {data.shape[0]} samples with {data.shape[1]} total columns")
    
    # Identify SNP feature columns (all columns except metadata)
    snp_cols = [c for c in data.columns if c not in [sample_col, binary_col, resid_col]]
    logger.info(f"Identified {len(snp_cols)} SNP features")
    
    # Check for missing values
    missing_counts = data[snp_cols].isnull().sum().sum()
    if missing_counts > 0:
        logger.warning(f"Found {missing_counts} missing values in SNP data")
    
    # Prepare features and labels
    X = data[snp_cols]
    y_resid = data[resid_col]
    y_binary = data[binary_col]
    
    # Log class distribution
    class_dist = y_binary.value_counts()
    logger.info(f"Binary phenotype distribution:\n{class_dist}")
    logger.info(f"Class balance: {class_dist.min() / class_dist.max():.3f}")
    
    return X, y_resid, y_binary, snp_cols


def create_optuna_objective(
    X: pd.DataFrame,
    y_resid: pd.Series,
    y_binary: pd.Series,
    snp_cols: List[str],
    n_splits: int,
    threshold: float,
    early_stopping_rounds: int,
    random_state: int
):
    """
    Create an Optuna objective function for hyperparameter optimization.
    
    Args:
        X: Feature matrix (SNP data).
        y_resid: Continuous residualized phenotype (target for training).
        y_binary: Binary phenotype labels (for evaluation).
        snp_cols: List of SNP column names for categorical features.
        n_splits: Number of K-fold cross-validation splits.
        threshold: Decision threshold to convert continuous predictions to binary.
        early_stopping_rounds: Number of rounds for early stopping.
        random_state: Random seed for reproducibility.
    
    Returns:
        Objective function for Optuna optimization.
    """
    logger = logging.getLogger(__name__)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    
    def objective(trial):
        """Optuna objective function to maximize mean AUROC across folds."""
        
        # Define hyperparameter search space
        params = {
            'boosting_type': 'gbdt',
            'objective': 'regression',
            'metric': 'rmse',
            'learning_rate': trial.suggest_float('learning_rate', 0.005, 0.2, log=True),
            'num_leaves': trial.suggest_int('num_leaves', 20, 550),
            'max_depth': trial.suggest_int('max_depth', 3, 15),
            'min_data_in_leaf': trial.suggest_int('min_data_in_leaf', 10, 100),
            'feature_fraction': trial.suggest_float('feature_fraction', 0.5, 1.0),
            'bagging_fraction': trial.suggest_float('bagging_fraction', 0.5, 1.0),
            'bagging_freq': 1,
            'lambda_l1': trial.suggest_float('lambda_l1', 1e-8, 10.0, log=True),
            'lambda_l2': trial.suggest_float('lambda_l2', 1e-8, 10.0, log=True),
            'max_bin': trial.suggest_int('max_bin', 128, 2048),
            'verbose': -1,
            'num_threads': -1,
            'seed': random_state
        }
        
        auroc_scores = []
        accuracy_scores = []
        
        # K-fold cross-validation
        for fold_idx, (train_idx, test_idx) in enumerate(kf.split(X), 1):
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train_cont = y_resid.iloc[train_idx]
            y_test_cont = y_resid.iloc[test_idx]
            y_test_bin = y_binary.iloc[test_idx]
            
            # Create LightGBM datasets
            train_data = lgb.Dataset(X_train, label=y_train_cont, categorical_feature=snp_cols)
            valid_data = lgb.Dataset(X_test, label=y_test_cont, reference=train_data)
            
            # Train model with early stopping
            model = lgb.train(
                params,
                train_data,
                num_boost_round=2000,
                valid_sets=[valid_data],
                callbacks=[
                    lgb.early_stopping(stopping_rounds=early_stopping_rounds, verbose=False),
                    lgb.log_evaluation(period=0)  # Suppress training logs
                ]
            )
            
            # Predict continuous values
            y_pred_cont = model.predict(X_test)
            
            # Convert to binary predictions using threshold
            y_pred_bin = (y_pred_cont >= threshold).astype(int)
            
            # Calculate metrics
            # ROC AUC (handle case where only one class is present)
            if len(np.unique(y_test_bin)) > 1:
                auc = roc_auc_score(y_test_bin, y_pred_cont)
                auroc_scores.append(auc)
            
            # Accuracy
            acc = accuracy_score(y_test_bin, y_pred_bin)
            accuracy_scores.append(acc)
        
        # Calculate mean metrics
        mean_auroc = np.mean(auroc_scores) if len(auroc_scores) > 0 else 0.0
        mean_accuracy = np.mean(accuracy_scores) if len(accuracy_scores) > 0 else 0.0
        
        # Log trial results
        trial.set_user_attr('mean_accuracy', mean_accuracy)
        trial.set_user_attr('std_auroc', np.std(auroc_scores) if len(auroc_scores) > 0 else 0.0)
        
        return mean_auroc
    
    return objective


def visualize_optimization(study: optuna.Study, output_dir: Path):
    """
    Create and save visualization plots for Optuna optimization results.
    
    Args:
        study: Completed Optuna study object.
        output_dir: Directory to save plots.
    """
    logger = logging.getLogger(__name__)
    logger.info("Creating optimization visualizations...")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Optimization history plot
    fig, ax = plt.subplots(figsize=(10, 6))
    trials = study.trials
    trial_numbers = [t.number for t in trials]
    values = [t.value for t in trials if t.value is not None]
    
    if values:
        ax.plot(trial_numbers[:len(values)], values, 'o-', alpha=0.6, label='Trial AUROC')
        
        # Plot running best
        best_so_far = []
        current_best = -np.inf
        for v in values:
            if v > current_best:
                current_best = v
            best_so_far.append(current_best)
        
        ax.plot(trial_numbers[:len(values)], best_so_far, 'r-', linewidth=2, label='Best AUROC')
        
        ax.set_xlabel('Trial Number', fontsize=12)
        ax.set_ylabel('Mean AUROC', fontsize=12)
        ax.set_title('Optimization History', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'optimization_history.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 2. Parameter importance plot
    try:
        fig = optuna.visualization.matplotlib.plot_param_importances(study)
        plt.tight_layout()
        plt.savefig(output_dir / 'param_importance.png', dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        logger.warning(f"Could not create parameter importance plot: {e}")
    
    # 3. Parallel coordinate plot for top trials
    try:
        fig = optuna.visualization.matplotlib.plot_parallel_coordinate(
            study, params=[
                'learning_rate', 'num_leaves', 'max_depth',
                'min_data_in_leaf', 'feature_fraction', 'bagging_fraction'
            ]
        )
        plt.tight_layout()
        plt.savefig(output_dir / 'parallel_coordinate.png', dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        logger.warning(f"Could not create parallel coordinate plot: {e}")
    
    logger.info(f"Visualizations saved to: {output_dir}")


def save_results(study: optuna.Study, output_dir: Path):
    """
    Save optimization results to files.
    
    Args:
        study: Completed Optuna study object.
        output_dir: Directory to save results.
    """
    logger = logging.getLogger(__name__)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save best parameters
    best_params_file = output_dir / 'best_parameters.txt'
    with open(best_params_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("BEST HYPERPARAMETERS\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Best Trial Number: {study.best_trial.number}\n")
        f.write(f"Best Mean AUROC: {study.best_value:.6f}\n\n")
        f.write("Best Parameters:\n")
        f.write("-" * 60 + "\n")
        for key, value in study.best_params.items():
            f.write(f"{key:20s}: {value}\n")
        
        if 'mean_accuracy' in study.best_trial.user_attrs:
            f.write("\n" + "-" * 60 + "\n")
            f.write("Additional Metrics:\n")
            f.write("-" * 60 + "\n")
            f.write(f"Mean Accuracy: {study.best_trial.user_attrs['mean_accuracy']:.6f}\n")
            if 'std_auroc' in study.best_trial.user_attrs:
                f.write(f"Std AUROC: {study.best_trial.user_attrs['std_auroc']:.6f}\n")
    
    logger.info(f"Best parameters saved to: {best_params_file}")
    
    # Save trials dataframe
    trials_df = study.trials_dataframe()
    trials_file = output_dir / 'all_trials.csv'
    trials_df.to_csv(trials_file, index=False)
    logger.info(f"All trials saved to: {trials_file}")
    
    # Save summary statistics
    summary_file = output_dir / 'optimization_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("OPTIMIZATION SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total Trials: {len(study.trials)}\n")
        f.write(f"Completed Trials: {len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])}\n")
        f.write(f"Pruned Trials: {len([t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED])}\n")
        f.write(f"Failed Trials: {len([t for t in study.trials if t.state == optuna.trial.TrialState.FAIL])}\n\n")
        
        completed_values = [t.value for t in study.trials if t.value is not None]
        if completed_values:
            f.write(f"Best AUROC: {max(completed_values):.6f}\n")
            f.write(f"Worst AUROC: {min(completed_values):.6f}\n")
            f.write(f"Mean AUROC: {np.mean(completed_values):.6f}\n")
            f.write(f"Median AUROC: {np.median(completed_values):.6f}\n")
            f.write(f"Std AUROC: {np.std(completed_values):.6f}\n")
    
    logger.info(f"Summary saved to: {summary_file}")


def main():
    """Main execution function."""
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='LightGBM hyperparameter optimization using Optuna',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input/Output arguments
    parser.add_argument('--data', type=str, required=True,
                        help='Path to input CSV file with SNP data and phenotypes')
    parser.add_argument('--output-dir', type=str, default='optuna_results',
                        help='Directory to save results')
    parser.add_argument('--log-file', type=str, default=None,
                        help='Path to log file (default: no file logging)')
    
    # Data column names
    parser.add_argument('--sample-col', type=str, default='sample',
                        help='Name of sample ID column')
    parser.add_argument('--binary-col', type=str, default='phenotype_binary',
                        help='Name of binary phenotype column')
    parser.add_argument('--resid-col', type=str, default='phenotype_resid',
                        help='Name of residualized continuous phenotype column')
    
    # Model parameters
    parser.add_argument('--threshold', type=float, default=-0.09547073,
                        help='Decision threshold to convert continuous predictions to binary')
    parser.add_argument('--n-splits', type=int, default=20,
                        help='Number of K-fold cross-validation splits')
    parser.add_argument('--early-stopping-rounds', type=int, default=50,
                        help='Early stopping rounds for LightGBM training')
    
    # Optuna parameters
    parser.add_argument('--n-trials', type=int, default=100,
                        help='Number of Optuna optimization trials')
    parser.add_argument('--study-name', type=str, default='lightgbm_optimization',
                        help='Name for Optuna study')
    parser.add_argument('--db-file', type=str, default='optuna_lgbm_study.db',
                        help='SQLite database file for Optuna study persistence')
    parser.add_argument('--load-if-exists', action='store_true',
                        help='Load existing study if it exists')
    
    # Other parameters
    parser.add_argument('--random-seed', type=int, default=42,
                        help='Random seed for reproducibility')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(log_file=args.log_file, verbose=args.verbose)
    
    logger.info("=" * 60)
    logger.info("LightGBM Hyperparameter Optimization with Optuna")
    logger.info("=" * 60)
    logger.info(f"Configuration:")
    for arg, value in vars(args).items():
        logger.info(f"  {arg}: {value}")
    logger.info("=" * 60)
    
    # Set random seeds for reproducibility
    np.random.seed(args.random_seed)
    
    # Load and prepare data
    try:
        X, y_resid, y_binary, snp_cols = load_and_prepare_data(
            args.data,
            args.sample_col,
            args.binary_col,
            args.resid_col
        )
    except Exception as e:
        logger.error(f"Failed to load data: {e}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Results will be saved to: {output_dir}")
    
    # Setup Optuna study
    storage = f"sqlite:///{args.db_file}"
    logger.info(f"Setting up Optuna study: {args.study_name}")
    logger.info(f"Study storage: {storage}")
    
    try:
        study = optuna.create_study(
            direction='maximize',
            study_name=args.study_name,
            storage=storage,
            load_if_exists=args.load_if_exists
        )
        
        if args.load_if_exists and len(study.trials) > 0:
            logger.info(f"Loaded existing study with {len(study.trials)} trials")
    except Exception as e:
        logger.error(f"Failed to create Optuna study: {e}")
        sys.exit(1)
    
    # Create objective function
    objective = create_optuna_objective(
        X, y_resid, y_binary, snp_cols,
        args.n_splits, args.threshold,
        args.early_stopping_rounds, args.random_seed
    )
    
    # Run optimization
    logger.info(f"Starting hyperparameter optimization for {args.n_trials} trials...")
    logger.info("This may take a while depending on data size and number of trials.")
    
    try:
        study.optimize(
            objective,
            n_trials=args.n_trials,
            show_progress_bar=args.verbose
        )
    except KeyboardInterrupt:
        logger.warning("Optimization interrupted by user")
    except Exception as e:
        logger.error(f"Optimization failed: {e}")
        sys.exit(1)
    
    # Log best results
    logger.info("=" * 60)
    logger.info("OPTIMIZATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Best trial: {study.best_trial.number}")
    logger.info(f"Best mean AUROC: {study.best_value:.6f}")
    logger.info("\nBest hyperparameters:")
    for key, value in study.best_params.items():
        logger.info(f"  {key}: {value}")
    
    # Save results
    logger.info("\nSaving results...")
    save_results(study, output_dir)
    
    # Create visualizations
    visualize_optimization(study, output_dir)
    
    logger.info("=" * 60)
    logger.info("All done! Check the output directory for results.")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
