#!/usr/bin/env python3
"""
Phenotype Classifier using Pre-trained LightGBM Model

This script uses a pre-trained LightGBM model to predict binary phenotypes from
SNP genotype data. It's designed to classify wild or unlabeled populations after
training a model on a labeled dataset.

Features:
- Load and validate pre-trained models
- Predict phenotypes for new samples
- Calculate prediction confidence scores
- Generate comprehensive reports and visualizations
- Handle VCA (variance component analysis) columns
- Support for batch predictions
- Flexible threshold adjustment

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
from joblib import load
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


def load_trained_model(model_path: str) -> lgb.Booster:
    """
    Load a pre-trained LightGBM model.
    
    Args:
        model_path: Path to the saved model file.
    
    Returns:
        Loaded LightGBM model.
    
    Raises:
        FileNotFoundError: If model file doesn't exist.
        Exception: If model loading fails.
    """
    logger = logging.getLogger(__name__)
    
    if not Path(model_path).exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    
    logger.info(f"Loading pre-trained model from: {model_path}")
    
    try:
        model = load(model_path)
        logger.info("Model loaded successfully")
        
        # Log model information if available
        if hasattr(model, 'num_trees'):
            logger.info(f"Model contains {model.num_trees()} trees")
        
        return model
    
    except Exception as e:
        logger.error(f"Failed to load model: {e}")
        raise


def load_and_validate_data(
    data_file: str,
    sample_col: str,
    vca_cols: Optional[List[str]],
    model: lgb.Booster
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], List[str]]:
    """
    Load and validate new data for prediction.
    
    Args:
        data_file: Path to CSV file with new samples.
        sample_col: Name of sample ID column.
        vca_cols: List of VCA column names (optional).
        model: Trained model to check feature compatibility.
    
    Returns:
        Tuple of (original data, feature matrix X, SNP columns, VCA columns).
    
    Raises:
        ValueError: If data validation fails.
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Loading data from: {data_file}")
    data = pd.read_csv(data_file)
    logger.info(f"Loaded {data.shape[0]} samples with {data.shape[1]} columns")
    
    # Handle VCA columns
    if vca_cols is None:
        vca_cols = []
    else:
        # Verify VCA columns exist
        missing_vca = [col for col in vca_cols if col not in data.columns]
        if missing_vca:
            logger.warning(f"VCA columns not found in data: {missing_vca}")
            vca_cols = [col for col in vca_cols if col in data.columns]
    
    # Identify SNP columns
    excluded_cols = [sample_col] + vca_cols
    snp_cols = [col for col in data.columns if col not in excluded_cols]
    logger.info(f"Identified {len(snp_cols)} SNP features")
    
    if vca_cols:
        logger.info(f"Including {len(vca_cols)} VCA columns: {vca_cols}")
    
    # Create feature matrix
    if vca_cols:
        X = pd.concat([data[vca_cols], data[snp_cols]], axis=1)
    else:
        X = data[snp_cols].copy()
    
    # Validate features match model
    model_features = model.feature_name()
    data_features = X.columns.tolist()
    
    # Check for missing features
    missing_features = set(model_features) - set(data_features)
    if missing_features:
        logger.error(f"Data is missing {len(missing_features)} features required by the model")
        logger.error(f"Missing features: {list(missing_features)[:10]}...")  # Show first 10
        raise ValueError(f"Data missing {len(missing_features)} required features")
    
    # Check for extra features (warn but don't fail)
    extra_features = set(data_features) - set(model_features)
    if extra_features:
        logger.warning(f"Data contains {len(extra_features)} extra features not used by model")
        # Reorder to match model features
        X = X[model_features]
    
    # Check for missing values
    missing_counts = X.isnull().sum().sum()
    if missing_counts > 0:
        logger.warning(f"Found {missing_counts} missing values in feature data")
        logger.warning("Missing values will be handled by LightGBM internally")
    
    logger.info("Data validation passed")
    
    return data, X, snp_cols, vca_cols


def make_predictions(
    model: lgb.Booster,
    X: pd.DataFrame,
    threshold: float = 0.5
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Make predictions using the trained model.
    
    Args:
        model: Trained LightGBM model.
        X: Feature matrix for prediction.
        threshold: Classification threshold for converting probabilities to binary.
    
    Returns:
        Tuple of (predicted probabilities, predicted classes).
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Making predictions for {len(X)} samples...")
    logger.info(f"Using classification threshold: {threshold}")
    
    # Get probability predictions
    y_pred_proba = model.predict(X)
    
    # Convert to binary predictions
    y_pred_class = (y_pred_proba >= threshold).astype(int)
    
    # Log prediction summary
    class_counts = pd.Series(y_pred_class).value_counts().sort_index()
    logger.info("Prediction summary:")
    for class_label, count in class_counts.items():
        logger.info(f"  Class {class_label}: {count} samples ({count/len(y_pred_class)*100:.2f}%)")
    
    return y_pred_proba, y_pred_class


def calculate_prediction_confidence(
    y_pred_proba: np.ndarray,
    threshold: float = 0.5
) -> pd.DataFrame:
    """
    Calculate various confidence metrics for predictions.
    
    Args:
        y_pred_proba: Predicted probabilities.
        threshold: Classification threshold.
    
    Returns:
        DataFrame with confidence metrics.
    """
    # Distance from threshold (how confident the prediction is)
    distance_from_threshold = np.abs(y_pred_proba - threshold)
    
    # Confidence score (0 to 1, where 1 is most confident)
    confidence = distance_from_threshold / 0.5
    
    # Prediction certainty categories
    def categorize_confidence(conf):
        if conf >= 0.8:
            return 'Very High'
        elif conf >= 0.6:
            return 'High'
        elif conf >= 0.4:
            return 'Moderate'
        elif conf >= 0.2:
            return 'Low'
        else:
            return 'Very Low'
    
    confidence_categories = [categorize_confidence(c) for c in confidence]
    
    return pd.DataFrame({
        'probability': y_pred_proba,
        'distance_from_threshold': distance_from_threshold,
        'confidence_score': confidence,
        'confidence_category': confidence_categories
    })


def create_results_dataframe(
    original_data: pd.DataFrame,
    sample_col: str,
    y_pred_proba: np.ndarray,
    y_pred_class: np.ndarray,
    confidence_metrics: pd.DataFrame
) -> pd.DataFrame:
    """
    Create comprehensive results dataframe.
    
    Args:
        original_data: Original input data.
        sample_col: Name of sample ID column.
        y_pred_proba: Predicted probabilities.
        y_pred_class: Predicted classes.
        confidence_metrics: DataFrame with confidence metrics.
    
    Returns:
        DataFrame with all results.
    """
    results = original_data.copy()
    
    # Add predictions
    results['predicted_phenotype'] = y_pred_class
    results['prediction_probability'] = y_pred_proba
    results['confidence_score'] = confidence_metrics['confidence_score']
    results['confidence_category'] = confidence_metrics['confidence_category']
    
    # Reorder columns to put predictions first (after sample ID)
    prediction_cols = [
        'predicted_phenotype',
        'prediction_probability', 
        'confidence_score',
        'confidence_category'
    ]
    
    other_cols = [col for col in results.columns if col not in prediction_cols + [sample_col]]
    
    results = results[[sample_col] + prediction_cols + other_cols]
    
    return results


def generate_summary_statistics(
    results: pd.DataFrame,
    confidence_metrics: pd.DataFrame
) -> Dict:
    """
    Generate summary statistics for predictions.
    
    Args:
        results: Results dataframe with predictions.
        confidence_metrics: Confidence metrics dataframe.
    
    Returns:
        Dictionary with summary statistics.
    """
    pred_counts = results['predicted_phenotype'].value_counts().sort_index()
    
    summary = {
        'total_samples': len(results),
        'class_0_count': pred_counts.get(0, 0),
        'class_1_count': pred_counts.get(1, 0),
        'class_0_percentage': (pred_counts.get(0, 0) / len(results)) * 100,
        'class_1_percentage': (pred_counts.get(1, 0) / len(results)) * 100,
        'mean_probability': results['prediction_probability'].mean(),
        'std_probability': results['prediction_probability'].std(),
        'mean_confidence': results['confidence_score'].mean(),
        'confidence_distribution': results['confidence_category'].value_counts().to_dict()
    }
    
    return summary


def visualize_predictions(
    results: pd.DataFrame,
    confidence_metrics: pd.DataFrame,
    summary_stats: Dict,
    output_dir: Path
):
    """
    Create visualizations for predictions.
    
    Args:
        results: Results dataframe with predictions.
        confidence_metrics: Confidence metrics dataframe.
        summary_stats: Summary statistics dictionary.
        output_dir: Directory to save plots.
    """
    logger = logging.getLogger(__name__)
    logger.info("Creating prediction visualizations...")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Prediction distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Bar plot of predicted classes
    class_counts = results['predicted_phenotype'].value_counts().sort_index()
    colors = ['#3498db', '#e74c3c']
    axes[0].bar(class_counts.index, class_counts.values, color=colors, alpha=0.8, edgecolor='black')
    axes[0].set_xlabel('Predicted Phenotype', fontsize=12)
    axes[0].set_ylabel('Number of Samples', fontsize=12)
    axes[0].set_title('Predicted Phenotype Distribution', fontsize=14, fontweight='bold')
    axes[0].set_xticks([0, 1])
    axes[0].set_xticklabels(['Class 0', 'Class 1'])
    
    # Add percentage labels on bars
    for i, (idx, count) in enumerate(class_counts.items()):
        percentage = (count / len(results)) * 100
        axes[0].text(idx, count, f'{count}\n({percentage:.1f}%)', 
                    ha='center', va='bottom', fontweight='bold')
    
    axes[0].grid(True, alpha=0.3, axis='y')
    
    # Probability distribution histogram
    axes[1].hist(results['prediction_probability'], bins=30, color='steelblue', 
                alpha=0.7, edgecolor='black')
    axes[1].axvline(0.5, color='red', linestyle='--', linewidth=2, label='Threshold (0.5)')
    axes[1].set_xlabel('Prediction Probability', fontsize=12)
    axes[1].set_ylabel('Frequency', fontsize=12)
    axes[1].set_title('Distribution of Prediction Probabilities', fontsize=14, fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'prediction_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Confidence analysis
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Confidence score distribution
    axes[0].hist(results['confidence_score'], bins=30, color='mediumseagreen', 
                alpha=0.7, edgecolor='black')
    axes[0].axvline(results['confidence_score'].mean(), color='red', linestyle='--', 
                   linewidth=2, label=f'Mean: {results["confidence_score"].mean():.3f}')
    axes[0].set_xlabel('Confidence Score', fontsize=12)
    axes[0].set_ylabel('Frequency', fontsize=12)
    axes[0].set_title('Distribution of Confidence Scores', fontsize=14, fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3, axis='y')
    
    # Confidence categories
    conf_cats = results['confidence_category'].value_counts()
    cat_order = ['Very High', 'High', 'Moderate', 'Low', 'Very Low']
    conf_cats = conf_cats.reindex([c for c in cat_order if c in conf_cats.index], fill_value=0)
    
    colors_conf = ['#27ae60', '#2ecc71', '#f39c12', '#e67e22', '#e74c3c']
    axes[1].barh(range(len(conf_cats)), conf_cats.values, 
                color=colors_conf[:len(conf_cats)], alpha=0.8, edgecolor='black')
    axes[1].set_yticks(range(len(conf_cats)))
    axes[1].set_yticklabels(conf_cats.index)
    axes[1].set_xlabel('Number of Samples', fontsize=12)
    axes[1].set_title('Prediction Confidence Categories', fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='x')
    
    # Add count labels on bars
    for i, count in enumerate(conf_cats.values):
        percentage = (count / len(results)) * 100
        axes[1].text(count, i, f' {count} ({percentage:.1f}%)', 
                    va='center', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'confidence_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Scatter plot: Probability vs Confidence
    fig, ax = plt.subplots(figsize=(10, 6))
    
    scatter = ax.scatter(results['prediction_probability'], 
                        results['confidence_score'],
                        c=results['predicted_phenotype'],
                        cmap='RdYlBu_r', alpha=0.6, s=50, edgecolor='black', linewidth=0.5)
    
    ax.axvline(0.5, color='gray', linestyle='--', linewidth=2, alpha=0.7, label='Threshold')
    ax.set_xlabel('Prediction Probability', fontsize=12)
    ax.set_ylabel('Confidence Score', fontsize=12)
    ax.set_title('Prediction Probability vs Confidence Score', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Predicted Class', fontsize=12)
    cbar.set_ticks([0.25, 0.75])
    cbar.set_ticklabels(['Class 0', 'Class 1'])
    
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'probability_vs_confidence.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Box plot by predicted class
    fig, ax = plt.subplots(figsize=(10, 6))
    
    data_to_plot = [
        results[results['predicted_phenotype'] == 0]['prediction_probability'],
        results[results['predicted_phenotype'] == 1]['prediction_probability']
    ]
    
    bp = ax.boxplot(data_to_plot, labels=['Class 0', 'Class 1'], patch_artist=True,
                   boxprops=dict(facecolor='lightblue', alpha=0.7),
                   medianprops=dict(color='red', linewidth=2),
                   whiskerprops=dict(linewidth=1.5),
                   capprops=dict(linewidth=1.5))
    
    ax.axhline(0.5, color='gray', linestyle='--', linewidth=2, alpha=0.7, label='Threshold')
    ax.set_ylabel('Prediction Probability', fontsize=12)
    ax.set_xlabel('Predicted Class', fontsize=12)
    ax.set_title('Probability Distribution by Predicted Class', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'probability_by_class.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Visualizations saved to: {output_dir}")


def save_results(
    results: pd.DataFrame,
    summary_stats: Dict,
    output_dir: Path,
    output_filename: str
):
    """
    Save prediction results and summary to files.
    
    Args:
        results: Results dataframe with predictions.
        summary_stats: Summary statistics dictionary.
        output_dir: Directory to save results.
        output_filename: Filename for results CSV.
    """
    logger = logging.getLogger(__name__)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save full results
    results_path = output_dir / output_filename
    results.to_csv(results_path, index=False)
    logger.info(f"Predictions saved to: {results_path}")
    
    # Save high-confidence predictions separately
    high_conf = results[results['confidence_category'].isin(['Very High', 'High'])]
    if len(high_conf) > 0:
        high_conf_path = output_dir / f"high_confidence_{output_filename}"
        high_conf.to_csv(high_conf_path, index=False)
        logger.info(f"High-confidence predictions saved to: {high_conf_path}")
    
    # Save low-confidence predictions for review
    low_conf = results[results['confidence_category'].isin(['Low', 'Very Low'])]
    if len(low_conf) > 0:
        low_conf_path = output_dir / f"low_confidence_{output_filename}"
        low_conf.to_csv(low_conf_path, index=False)
        logger.info(f"Low-confidence predictions (for review) saved to: {low_conf_path}")
    
    # Save summary report
    summary_path = output_dir / 'prediction_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("PHENOTYPE PREDICTION SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Overall Statistics:\n")
        f.write("-" * 60 + "\n")
        f.write(f"Total samples processed: {summary_stats['total_samples']}\n\n")
        
        f.write("Predicted Phenotype Distribution:\n")
        f.write("-" * 60 + "\n")
        f.write(f"Class 0: {summary_stats['class_0_count']} samples ({summary_stats['class_0_percentage']:.2f}%)\n")
        f.write(f"Class 1: {summary_stats['class_1_count']} samples ({summary_stats['class_1_percentage']:.2f}%)\n\n")
        
        f.write("Prediction Confidence:\n")
        f.write("-" * 60 + "\n")
        f.write(f"Mean probability: {summary_stats['mean_probability']:.4f} Â± {summary_stats['std_probability']:.4f}\n")
        f.write(f"Mean confidence score: {summary_stats['mean_confidence']:.4f}\n\n")
        
        f.write("Confidence Distribution:\n")
        f.write("-" * 60 + "\n")
        for category, count in summary_stats['confidence_distribution'].items():
            percentage = (count / summary_stats['total_samples']) * 100
            f.write(f"{category:15s}: {count:4d} samples ({percentage:5.2f}%)\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("Notes:\n")
        f.write("- High confidence predictions (Very High/High) are recommended for use\n")
        f.write("- Low confidence predictions (Low/Very Low) should be verified\n")
        f.write("- Confidence score ranges from 0 (uncertain) to 1 (certain)\n")
    
    logger.info(f"Summary report saved to: {summary_path}")
    
    # Save simplified summary for quick viewing
    quick_summary_path = output_dir / 'quick_summary.txt'
    with open(quick_summary_path, 'w') as f:
        f.write("QUICK SUMMARY\n")
        f.write("=" * 40 + "\n")
        f.write(f"Total: {summary_stats['total_samples']} samples\n")
        f.write(f"Class 0: {summary_stats['class_0_count']} ({summary_stats['class_0_percentage']:.1f}%)\n")
        f.write(f"Class 1: {summary_stats['class_1_count']} ({summary_stats['class_1_percentage']:.1f}%)\n")
        f.write(f"Mean confidence: {summary_stats['mean_confidence']:.3f}\n")
    
    logger.info(f"Quick summary saved to: {quick_summary_path}")


def main():
    """Main execution function."""
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Predict phenotypes using pre-trained LightGBM model',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--model', type=str, required=True,
                        help='Path to pre-trained model file (.pkl)')
    parser.add_argument('--data', type=str, required=True,
                        help='Path to CSV file with samples to classify')
    
    # Output arguments
    parser.add_argument('--output-dir', type=str, default='predictions',
                        help='Directory to save prediction results')
    parser.add_argument('--output-filename', type=str, default='classified_samples.csv',
                        help='Filename for prediction results')
    parser.add_argument('--log-file', type=str, default=None,
                        help='Path to log file (default: no file logging)')
    
    # Data configuration
    parser.add_argument('--sample-col', type=str, default='sample',
                        help='Name of sample ID column')
    parser.add_argument('--vca-cols', type=str, nargs='*', default=None,
                        help='Names of VCA columns (space-separated, if used in training)')
    
    # Prediction parameters
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Classification threshold for converting probabilities to binary')
    
    # Flags
    parser.add_argument('--no-visualizations', action='store_true',
                        help='Skip generating visualization plots')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(log_file=args.log_file, verbose=args.verbose)
    
    logger.info("=" * 60)
    logger.info("Phenotype Prediction using Pre-trained Model")
    logger.info("=" * 60)
    logger.info("Configuration:")
    for arg, value in vars(args).items():
        logger.info(f"  {arg}: {value}")
    logger.info("=" * 60)
    
    # Load trained model
    try:
        model = load_trained_model(args.model)
    except Exception as e:
        logger.error(f"Failed to load model: {e}")
        sys.exit(1)
    
    # Load and validate data
    try:
        original_data, X, snp_cols, vca_cols = load_and_validate_data(
            args.data,
            args.sample_col,
            args.vca_cols,
            model
        )
    except Exception as e:
        logger.error(f"Data validation failed: {e}")
        sys.exit(1)
    
    # Make predictions
    try:
        y_pred_proba, y_pred_class = make_predictions(model, X, args.threshold)
    except Exception as e:
        logger.error(f"Prediction failed: {e}")
        sys.exit(1)
    
    # Calculate confidence metrics
    logger.info("Calculating prediction confidence metrics...")
    confidence_metrics = calculate_prediction_confidence(y_pred_proba, args.threshold)
    
    # Create results dataframe
    logger.info("Creating results dataframe...")
    results = create_results_dataframe(
        original_data,
        args.sample_col,
        y_pred_proba,
        y_pred_class,
        confidence_metrics
    )
    
    # Generate summary statistics
    logger.info("Generating summary statistics...")
    summary_stats = generate_summary_statistics(results, confidence_metrics)
    
    # Display summary
    logger.info("\n" + "=" * 60)
    logger.info("PREDICTION RESULTS")
    logger.info("=" * 60)
    logger.info(f"Total samples: {summary_stats['total_samples']}")
    logger.info(f"Class 0: {summary_stats['class_0_count']} ({summary_stats['class_0_percentage']:.2f}%)")
    logger.info(f"Class 1: {summary_stats['class_1_count']} ({summary_stats['class_1_percentage']:.2f}%)")
    logger.info(f"Mean confidence: {summary_stats['mean_confidence']:.4f}")
    logger.info("\nConfidence distribution:")
    for category, count in summary_stats['confidence_distribution'].items():
        percentage = (count / summary_stats['total_samples']) * 100
        logger.info(f"  {category}: {count} ({percentage:.1f}%)")
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save results
    logger.info("\nSaving results...")
    save_results(results, summary_stats, output_dir, args.output_filename)
    
    # Create visualizations
    if not args.no_visualizations:
        visualize_predictions(results, confidence_metrics, summary_stats, output_dir)
    else:
        logger.info("Skipping visualization generation (--no-visualizations flag set)")
    
    logger.info("=" * 60)
    logger.info("Prediction complete! Check output directory for results.")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
