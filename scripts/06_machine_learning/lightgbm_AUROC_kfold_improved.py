#!/usr/bin/env python3
"""
LightGBM K-Fold Cross-Validation with AUROC Evaluation

This script performs stratified k-fold cross-validation using LightGBM for binary
classification on SNP data. It evaluates model performance using AUROC, accuracy,
and generates comprehensive classification reports with visualizations.

Features:
- Stratified k-fold cross-validation for imbalanced datasets
- Support for variance component analysis (VCA) columns
- Automatic class imbalance handling with scale_pos_weight
- Early stopping to prevent overfitting
- Comprehensive metrics tracking per fold
- Feature importance aggregation and visualization
- Model persistence with reproducible results

Author: [Your Name]
Date: 2025-10-30
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    accuracy_score, 
    classification_report, 
    roc_auc_score,
    confusion_matrix,
    precision_recall_curve,
    roc_curve
)
from joblib import dump
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# Set plotting style
sns.set_style("whitegrid")


def setup_logging(log_file: Optional[str] = None, verbose: bool = True) -> logging.Logger:
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
    sample_col: str,
    phenotype_col: str,
    vca_cols: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, pd.Series, List[str], List[str]]:
    """
    Load and prepare data for model training.
    
    Args:
        data_file: Path to CSV file containing SNP data and phenotypes.
        sample_col: Name of sample ID column.
        phenotype_col: Name of binary phenotype column.
        vca_cols: List of variance component analysis column names (optional).
    
    Returns:
        Tuple of (X features, y labels, SNP column names, VCA column names).
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Loading data from: {data_file}")
    data = pd.read_csv(data_file)
    logger.info(f"Loaded {data.shape[0]} samples with {data.shape[1]} total columns")
    
    # Handle VCA columns
    if vca_cols is None:
        vca_cols = []
    else:
        # Verify VCA columns exist
        missing_vca = [col for col in vca_cols if col not in data.columns]
        if missing_vca:
            logger.warning(f"VCA columns not found in data: {missing_vca}")
            vca_cols = [col for col in vca_cols if col in data.columns]
    
    # Identify SNP feature columns
    excluded_cols = [sample_col, phenotype_col] + vca_cols
    snp_cols = [col for col in data.columns if col not in excluded_cols]
    logger.info(f"Identified {len(snp_cols)} SNP features")
    
    if vca_cols:
        logger.info(f"Including {len(vca_cols)} VCA columns: {vca_cols}")
    
    # Check for missing values
    missing_counts = data[snp_cols + vca_cols].isnull().sum().sum()
    if missing_counts > 0:
        logger.warning(f"Found {missing_counts} missing values in feature data")
    
    # Prepare features
    if vca_cols:
        X = pd.concat([data[vca_cols], data[snp_cols]], axis=1)
    else:
        X = data[snp_cols].copy()
    
    y = data[phenotype_col]
    
    # Check phenotype is binary
    unique_values = y.unique()
    if not set(unique_values).issubset({0, 1}):
        logger.error(f"Phenotype must be binary (0/1), found values: {unique_values}")
        raise ValueError("Phenotype must be binary (0/1)")
    
    # Log class distribution
    class_counts = y.value_counts().sort_index()
    logger.info(f"Class distribution:")
    for class_label, count in class_counts.items():
        logger.info(f"  Class {class_label}: {count} samples ({count/len(y)*100:.2f}%)")
    
    # Calculate and log class imbalance ratio
    imbalance_ratio = class_counts.min() / class_counts.max()
    logger.info(f"Class imbalance ratio: {imbalance_ratio:.3f}")
    
    if imbalance_ratio < 0.5:
        logger.warning("Significant class imbalance detected. Using scale_pos_weight for correction.")
    
    return X, y, snp_cols, vca_cols


def create_lgb_params(
    scale_pos_weight: float,
    custom_params: Optional[Dict] = None,
    random_seed: int = 42
) -> Dict:
    """
    Create LightGBM parameters dictionary.
    
    Args:
        scale_pos_weight: Weight for positive class to handle imbalance.
        custom_params: Optional dictionary of custom parameters to override defaults.
        random_seed: Random seed for reproducibility.
    
    Returns:
        Dictionary of LightGBM parameters.
    """
    # Default parameters (can be overridden by custom_params)
    params = {
        'boosting_type': 'gbdt',
        'objective': 'binary',
        'metric': 'binary_logloss',
        'learning_rate': 0.005842397819147727,
        'num_leaves': 396,
        'max_depth': 6,
        'min_data_in_leaf': 53,
        'feature_fraction': 0.7151150335956341,
        'bagging_fraction': 0.9735126336808713,
        'bagging_freq': 1,
        'lambda_l1': 0.31704366636190856,
        'lambda_l2': 0.000010004695107799964,
        'max_bin': 183,
        'verbose': -1,
        'scale_pos_weight': scale_pos_weight,
        'seed': random_seed
    }
    
    # Override with custom parameters if provided
    if custom_params:
        params.update(custom_params)
    
    return params


def train_and_evaluate_fold(
    fold_num: int,
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    y_train: pd.Series,
    y_test: pd.Series,
    params: Dict,
    snp_cols: List[str],
    num_boost_round: int,
    early_stopping_rounds: Optional[int],
    threshold: float = 0.5
) -> Dict:
    """
    Train and evaluate model for a single fold.
    
    Args:
        fold_num: Fold number for logging.
        X_train: Training features.
        X_test: Test features.
        y_train: Training labels.
        y_test: Test labels.
        params: LightGBM parameters.
        snp_cols: List of SNP column names (for categorical feature specification).
        num_boost_round: Maximum number of boosting rounds.
        early_stopping_rounds: Early stopping rounds (None to disable).
        threshold: Classification threshold for converting probabilities to binary.
    
    Returns:
        Dictionary containing fold results and trained model.
    """
    logger = logging.getLogger(__name__)
    
    # Create LightGBM datasets with categorical features
    train_data = lgb.Dataset(X_train, label=y_train, categorical_feature=snp_cols)
    valid_data = lgb.Dataset(X_test, label=y_test, reference=train_data, categorical_feature=snp_cols)
    
    # Setup callbacks
    callbacks = [lgb.log_evaluation(period=0)]  # Suppress training logs
    if early_stopping_rounds:
        callbacks.append(lgb.early_stopping(stopping_rounds=early_stopping_rounds, verbose=False))
    
    # Train model
    model = lgb.train(
        params,
        train_data,
        valid_sets=[train_data, valid_data],
        num_boost_round=num_boost_round,
        callbacks=callbacks
    )
    
    # Get predictions
    y_pred_proba = model.predict(X_test)
    y_pred_binary = (y_pred_proba >= threshold).astype(int)
    
    # Calculate metrics
    accuracy = accuracy_score(y_test, y_pred_binary)
    auroc = roc_auc_score(y_test, y_pred_proba)
    conf_matrix = confusion_matrix(y_test, y_pred_binary)
    
    # Get classification report
    report = classification_report(y_test, y_pred_binary, output_dict=True, zero_division=0)
    
    # Get feature importance (split by feature type)
    features = model.feature_name()
    importances = model.feature_importance(importance_type='gain')
    importance_df = pd.DataFrame({'Feature': features, 'Importance': importances})
    
    # Separate SNP importance from VCA importance
    snp_importance = importance_df[importance_df['Feature'].isin(snp_cols)]
    vca_importance = importance_df[~importance_df['Feature'].isin(snp_cols)]
    
    logger.info(f"Fold {fold_num}: Accuracy={accuracy:.4f}, AUROC={auroc:.4f}")
    
    return {
        'fold': fold_num,
        'model': model,
        'accuracy': accuracy,
        'auroc': auroc,
        'confusion_matrix': conf_matrix,
        'classification_report': report,
        'predictions_proba': y_pred_proba,
        'predictions_binary': y_pred_binary,
        'true_labels': y_test,
        'snp_importance': snp_importance,
        'vca_importance': vca_importance,
        'best_iteration': model.best_iteration if early_stopping_rounds else num_boost_round
    }


def aggregate_metrics(fold_results: List[Dict]) -> Dict:
    """
    Aggregate metrics across all folds.
    
    Args:
        fold_results: List of fold result dictionaries.
    
    Returns:
        Dictionary containing aggregated metrics.
    """
    logger = logging.getLogger(__name__)
    
    accuracies = [r['accuracy'] for r in fold_results]
    aurocs = [r['auroc'] for r in fold_results]
    
    # Aggregate classification reports
    all_reports = [r['classification_report'] for r in fold_results]
    average_report = {}
    
    for key in all_reports[0].keys():
        if isinstance(all_reports[0][key], dict):
            nested_avg = pd.DataFrame([report[key] for report in all_reports]).mean(axis=0).to_dict()
            average_report[key] = nested_avg
        else:
            average_report[key] = np.mean([report[key] for report in all_reports])
    
    # Aggregate confusion matrices
    total_conf_matrix = sum(r['confusion_matrix'] for r in fold_results)
    
    # Aggregate feature importance
    snp_importance_list = [r['snp_importance'] for r in fold_results]
    final_snp_importance = pd.concat(snp_importance_list).groupby('Feature').mean()
    final_snp_importance = final_snp_importance.sort_values(by='Importance', ascending=False)
    
    vca_importance_list = [r['vca_importance'] for r in fold_results if not r['vca_importance'].empty]
    if vca_importance_list:
        final_vca_importance = pd.concat(vca_importance_list).groupby('Feature').mean()
        final_vca_importance = final_vca_importance.sort_values(by='Importance', ascending=False)
    else:
        final_vca_importance = None
    
    logger.info("=" * 60)
    logger.info("CROSS-VALIDATION RESULTS")
    logger.info("=" * 60)
    logger.info(f"Mean Accuracy: {np.mean(accuracies):.4f} Â± {np.std(accuracies):.4f}")
    logger.info(f"Mean AUROC: {np.mean(aurocs):.4f} Â± {np.std(aurocs):.4f}")
    logger.info(f"Min AUROC: {np.min(aurocs):.4f}")
    logger.info(f"Max AUROC: {np.max(aurocs):.4f}")
    
    return {
        'accuracies': accuracies,
        'aurocs': aurocs,
        'mean_accuracy': np.mean(accuracies),
        'std_accuracy': np.std(accuracies),
        'mean_auroc': np.mean(aurocs),
        'std_auroc': np.std(aurocs),
        'average_report': average_report,
        'total_confusion_matrix': total_conf_matrix,
        'snp_importance': final_snp_importance,
        'vca_importance': final_vca_importance
    }


def visualize_results(
    fold_results: List[Dict],
    aggregated_metrics: Dict,
    output_dir: Path
):
    """
    Create comprehensive visualizations of cross-validation results.
    
    Args:
        fold_results: List of fold result dictionaries.
        aggregated_metrics: Dictionary of aggregated metrics.
        output_dir: Directory to save plots.
    """
    logger = logging.getLogger(__name__)
    logger.info("Creating visualizations...")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Metrics across folds
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    folds = [r['fold'] for r in fold_results]
    accuracies = [r['accuracy'] for r in fold_results]
    aurocs = [r['auroc'] for r in fold_results]
    
    # Accuracy plot
    axes[0].plot(folds, accuracies, 'o-', color='steelblue', linewidth=2, markersize=6)
    axes[0].axhline(np.mean(accuracies), color='red', linestyle='--', 
                    label=f'Mean: {np.mean(accuracies):.4f}')
    axes[0].fill_between(folds, 
                         np.mean(accuracies) - np.std(accuracies),
                         np.mean(accuracies) + np.std(accuracies),
                         alpha=0.2, color='red')
    axes[0].set_xlabel('Fold Number', fontsize=12)
    axes[0].set_ylabel('Accuracy', fontsize=12)
    axes[0].set_title('Accuracy Across Folds', fontsize=14, fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # AUROC plot
    axes[1].plot(folds, aurocs, 'o-', color='darkorange', linewidth=2, markersize=6)
    axes[1].axhline(np.mean(aurocs), color='red', linestyle='--',
                    label=f'Mean: {np.mean(aurocs):.4f}')
    axes[1].fill_between(folds,
                         np.mean(aurocs) - np.std(aurocs),
                         np.mean(aurocs) + np.std(aurocs),
                         alpha=0.2, color='red')
    axes[1].set_xlabel('Fold Number', fontsize=12)
    axes[1].set_ylabel('AUROC', fontsize=12)
    axes[1].set_title('AUROC Across Folds', fontsize=14, fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'metrics_across_folds.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Confusion matrix (aggregated)
    fig, ax = plt.subplots(figsize=(8, 6))
    conf_matrix = aggregated_metrics['total_confusion_matrix']
    sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', 
                xticklabels=['Class 0', 'Class 1'],
                yticklabels=['Class 0', 'Class 1'],
                ax=ax, cbar_kws={'label': 'Count'})
    ax.set_xlabel('Predicted Label', fontsize=12)
    ax.set_ylabel('True Label', fontsize=12)
    ax.set_title('Aggregated Confusion Matrix', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / 'confusion_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. ROC curves (sample of folds to avoid clutter)
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot ROC curve for every 5th fold to avoid overcrowding
    fold_sample = fold_results[::max(1, len(fold_results)//10)]
    
    for result in fold_sample:
        fpr, tpr, _ = roc_curve(result['true_labels'], result['predictions_proba'])
        ax.plot(fpr, tpr, alpha=0.3, linewidth=1)
    
    # Plot mean ROC curve
    all_fpr = np.linspace(0, 1, 100)
    mean_tpr = np.zeros_like(all_fpr)
    
    for result in fold_results:
        fpr, tpr, _ = roc_curve(result['true_labels'], result['predictions_proba'])
        mean_tpr += np.interp(all_fpr, fpr, tpr)
    
    mean_tpr /= len(fold_results)
    mean_auroc = aggregated_metrics['mean_auroc']
    
    ax.plot(all_fpr, mean_tpr, color='red', linewidth=3, 
            label=f'Mean ROC (AUC = {mean_auroc:.4f})')
    ax.plot([0, 1], [0, 1], 'k--', linewidth=2, label='Random Classifier')
    
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title('ROC Curves Across Folds', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'roc_curves.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Feature importance (top 20 SNPs)
    fig, ax = plt.subplots(figsize=(10, 8))
    top_features = aggregated_metrics['snp_importance'].head(20)
    
    colors = sns.color_palette("viridis", len(top_features))
    ax.barh(range(len(top_features)), top_features['Importance'], color=colors)
    ax.set_yticks(range(len(top_features)))
    ax.set_yticklabels(top_features.index, fontsize=10)
    ax.set_xlabel('Mean Importance (Gain)', fontsize=12)
    ax.set_title('Top 20 SNP Feature Importance', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'feature_importance_top20.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. VCA importance (if applicable)
    if aggregated_metrics['vca_importance'] is not None and not aggregated_metrics['vca_importance'].empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        vca_imp = aggregated_metrics['vca_importance']
        
        ax.bar(range(len(vca_imp)), vca_imp['Importance'], color='coral')
        ax.set_xticks(range(len(vca_imp)))
        ax.set_xticklabels(vca_imp.index, fontsize=11)
        ax.set_ylabel('Mean Importance (Gain)', fontsize=12)
        ax.set_title('VCA Feature Importance', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'vca_importance.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    logger.info(f"Visualizations saved to: {output_dir}")


def save_results(
    fold_results: List[Dict],
    aggregated_metrics: Dict,
    final_model: lgb.Booster,
    params: Dict,
    output_dir: Path,
    model_filename: str,
    feature_importance_filename: str
):
    """
    Save all results to files.
    
    Args:
        fold_results: List of fold result dictionaries.
        aggregated_metrics: Dictionary of aggregated metrics.
        final_model: Trained model from the last fold.
        params: LightGBM parameters used.
        output_dir: Directory to save results.
        model_filename: Filename for saved model.
        feature_importance_filename: Filename for feature importance CSV.
    """
    logger = logging.getLogger(__name__)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save final model
    model_path = output_dir / model_filename
    dump(final_model, model_path)
    logger.info(f"Final model saved to: {model_path}")
    
    # Save SNP feature importance
    importance_path = output_dir / feature_importance_filename
    aggregated_metrics['snp_importance'].to_csv(importance_path)
    logger.info(f"SNP feature importance saved to: {importance_path}")
    
    # Save VCA feature importance if applicable
    if aggregated_metrics['vca_importance'] is not None:
        vca_importance_path = output_dir / 'vca_feature_importance.csv'
        aggregated_metrics['vca_importance'].to_csv(vca_importance_path)
        logger.info(f"VCA feature importance saved to: {vca_importance_path}")
    
    # Save metrics per fold
    fold_metrics = pd.DataFrame([
        {
            'Fold': r['fold'],
            'Accuracy': r['accuracy'],
            'AUROC': r['auroc'],
            'Best_Iteration': r['best_iteration']
        }
        for r in fold_results
    ])
    fold_metrics_path = output_dir / 'fold_metrics.csv'
    fold_metrics.to_csv(fold_metrics_path, index=False)
    logger.info(f"Per-fold metrics saved to: {fold_metrics_path}")
    
    # Save aggregated classification report
    report_df = pd.DataFrame(aggregated_metrics['average_report']).transpose()
    report_path = output_dir / 'classification_report.csv'
    report_df.to_csv(report_path)
    logger.info(f"Classification report saved to: {report_path}")
    
    # Save summary
    summary_path = output_dir / 'cv_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("K-FOLD CROSS-VALIDATION SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Overall Performance:\n")
        f.write("-" * 60 + "\n")
        f.write(f"Number of Folds: {len(fold_results)}\n")
        f.write(f"Mean Accuracy: {aggregated_metrics['mean_accuracy']:.6f} Â± {aggregated_metrics['std_accuracy']:.6f}\n")
        f.write(f"Mean AUROC: {aggregated_metrics['mean_auroc']:.6f} Â± {aggregated_metrics['std_auroc']:.6f}\n")
        f.write(f"Min AUROC: {min(aggregated_metrics['aurocs']):.6f}\n")
        f.write(f"Max AUROC: {max(aggregated_metrics['aurocs']):.6f}\n\n")
        
        f.write("Model Parameters:\n")
        f.write("-" * 60 + "\n")
        for key, value in params.items():
            f.write(f"{key:25s}: {value}\n")
        f.write("\n")
        
        f.write("Aggregated Confusion Matrix:\n")
        f.write("-" * 60 + "\n")
        conf_matrix = aggregated_metrics['total_confusion_matrix']
        f.write(f"                Predicted 0    Predicted 1\n")
        f.write(f"True 0          {conf_matrix[0, 0]:11d}    {conf_matrix[0, 1]:11d}\n")
        f.write(f"True 1          {conf_matrix[1, 0]:11d}    {conf_matrix[1, 1]:11d}\n\n")
        
        f.write("Top 10 SNP Features by Importance:\n")
        f.write("-" * 60 + "\n")
        top_features = aggregated_metrics['snp_importance'].head(10)
        for idx, (feature, importance) in enumerate(top_features.items(), 1):
            f.write(f"{idx:2d}. {feature:30s}: {importance:.6f}\n")
        
        if aggregated_metrics['vca_importance'] is not None:
            f.write("\nVCA Feature Importance:\n")
            f.write("-" * 60 + "\n")
            for feature, importance in aggregated_metrics['vca_importance'].items():
                f.write(f"{feature:30s}: {importance:.6f}\n")
    
    logger.info(f"Summary saved to: {summary_path}")


def main():
    """Main execution function."""
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='LightGBM k-fold cross-validation for binary classification',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input/Output arguments
    parser.add_argument('--data', type=str, required=True,
                        help='Path to input CSV file with SNP data and phenotypes')
    parser.add_argument('--output-dir', type=str, default='cv_results',
                        help='Directory to save results')
    parser.add_argument('--model-name', type=str, default='best_model_kfold_lgbm.pkl',
                        help='Filename for saved model')
    parser.add_argument('--feature-importance-name', type=str, default='snp_feature_importance.csv',
                        help='Filename for SNP feature importance CSV')
    parser.add_argument('--log-file', type=str, default=None,
                        help='Path to log file (default: no file logging)')
    
    # Data column names
    parser.add_argument('--sample-col', type=str, default='sample',
                        help='Name of sample ID column')
    parser.add_argument('--phenotype-col', type=str, default='phenotype',
                        help='Name of binary phenotype column')
    parser.add_argument('--vca-cols', type=str, nargs='*', default=None,
                        help='Names of variance component analysis columns (space-separated)')
    
    # Model parameters
    parser.add_argument('--n-splits', type=int, default=50,
                        help='Number of cross-validation folds')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Classification threshold for converting probabilities to binary')
    parser.add_argument('--num-boost-round', type=int, default=2000,
                        help='Maximum number of boosting rounds')
    parser.add_argument('--early-stopping-rounds', type=int, default=None,
                        help='Early stopping rounds (None to disable)')
    
    # LightGBM hyperparameters (optional overrides)
    parser.add_argument('--learning-rate', type=float, default=None,
                        help='Learning rate (overrides default)')
    parser.add_argument('--num-leaves', type=int, default=None,
                        help='Number of leaves (overrides default)')
    parser.add_argument('--max-depth', type=int, default=None,
                        help='Maximum tree depth (overrides default)')
    
    # Other parameters
    parser.add_argument('--random-seed', type=int, default=42,
                        help='Random seed for reproducibility')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(log_file=args.log_file, verbose=args.verbose)
    
    logger.info("=" * 60)
    logger.info("LightGBM K-Fold Cross-Validation")
    logger.info("=" * 60)
    logger.info("Configuration:")
    for arg, value in vars(args).items():
        logger.info(f"  {arg}: {value}")
    logger.info("=" * 60)
    
    # Set random seeds for reproducibility
    np.random.seed(args.random_seed)
    
    # Load and prepare data
    try:
        X, y, snp_cols, vca_cols = load_and_prepare_data(
            args.data,
            args.sample_col,
            args.phenotype_col,
            args.vca_cols
        )
    except Exception as e:
        logger.error(f"Failed to load data: {e}")
        sys.exit(1)
    
    # Calculate scale_pos_weight for class imbalance
    num_class_0 = (y == 0).sum()
    num_class_1 = (y == 1).sum()
    scale_pos_weight = num_class_0 / num_class_1
    logger.info(f"Calculated scale_pos_weight: {scale_pos_weight:.4f}")
    
    # Create LightGBM parameters
    custom_params = {}
    if args.learning_rate is not None:
        custom_params['learning_rate'] = args.learning_rate
    if args.num_leaves is not None:
        custom_params['num_leaves'] = args.num_leaves
    if args.max_depth is not None:
        custom_params['max_depth'] = args.max_depth
    
    params = create_lgb_params(scale_pos_weight, custom_params, args.random_seed)
    
    # Log parameters being used
    logger.info("\nLightGBM Parameters:")
    for key, value in params.items():
        logger.info(f"  {key}: {value}")
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"\nResults will be saved to: {output_dir}")
    
    # Perform k-fold cross-validation
    logger.info(f"\nStarting {args.n_splits}-fold cross-validation...")
    kf = StratifiedKFold(n_splits=args.n_splits, shuffle=True, random_state=args.random_seed)
    
    fold_results = []
    
    try:
        for fold_num, (train_idx, test_idx) in enumerate(kf.split(X, y), 1):
            logger.info(f"\nProcessing Fold {fold_num}/{args.n_splits}...")
            
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
            
            fold_result = train_and_evaluate_fold(
                fold_num,
                X_train, X_test,
                y_train, y_test,
                params,
                snp_cols,
                args.num_boost_round,
                args.early_stopping_rounds,
                args.threshold
            )
            
            fold_results.append(fold_result)
    
    except KeyboardInterrupt:
        logger.warning("\nCross-validation interrupted by user")
        if len(fold_results) == 0:
            logger.error("No folds completed. Exiting.")
            sys.exit(1)
        logger.info(f"Continuing with {len(fold_results)} completed folds...")
    
    except Exception as e:
        logger.error(f"Cross-validation failed: {e}")
        sys.exit(1)
    
    # Aggregate metrics
    logger.info("\nAggregating metrics across folds...")
    aggregated_metrics = aggregate_metrics(fold_results)
    
    # Display classification report
    logger.info("\nAverage Classification Report:")
    report_df = pd.DataFrame(aggregated_metrics['average_report']).transpose()
    logger.info(f"\n{report_df.to_string()}")
    
    # Display top features
    logger.info("\nTop 10 SNP Features by Importance:")
    top_snps = aggregated_metrics['snp_importance'].head(10)
    for idx, (feature, row) in enumerate(top_snps.iterrows(), 1):
        logger.info(f"  {idx:2d}. {feature:30s}: {row['Importance']:.6f}")
    
    if aggregated_metrics['vca_importance'] is not None:
        logger.info("\nVCA Feature Importance:")
        for feature, row in aggregated_metrics['vca_importance'].iterrows():
            logger.info(f"  {feature:30s}: {row['Importance']:.6f}")
    
    # Save results
    logger.info("\nSaving results...")
    final_model = fold_results[-1]['model']  # Model from last fold
    save_results(
        fold_results,
        aggregated_metrics,
        final_model,
        params,
        output_dir,
        args.model_name,
        args.feature_importance_name
    )
    
    # Create visualizations
    visualize_results(fold_results, aggregated_metrics, output_dir)
    
    logger.info("=" * 60)
    logger.info("Cross-validation complete! Check output directory for results.")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
