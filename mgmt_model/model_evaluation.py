"""
Model evaluation module for MGMT methylation detection.
Provides comprehensive evaluation metrics, visualizations, and validation.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Any
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, roc_curve, precision_recall_curve, confusion_matrix,
    classification_report, average_precision_score
)
from sklearn.model_selection import (
    cross_val_score, cross_validate, learning_curve,
    validation_curve, StratifiedKFold
)
from sklearn.calibration import calibration_curve, CalibratedClassifierCV
import joblib
from pathlib import Path
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ModelEvaluator:
    """
    Comprehensive model evaluation for MGMT methylation prediction.
    """
    
    def __init__(self, output_dir: str = "evaluation_results"):
        """
        Initialize the evaluator.
        
        Args:
            output_dir (str): Directory to save evaluation results
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set plotting style
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
    
    def evaluate_single_model(self, model: Any, X_test: np.ndarray, y_test: np.ndarray,
                             model_name: str = "model") -> Dict[str, float]:
        """
        Comprehensive evaluation of a single model.
        
        Args:
            model: Trained model
            X_test (np.ndarray): Test features
            y_test (np.ndarray): Test labels
            model_name (str): Name of the model
            
        Returns:
            Dict: Comprehensive evaluation metrics
        """
        logger.info(f"Evaluating {model_name} model")
        
        # Get predictions
        y_pred = model.predict(X_test)
        y_pred_proba = model.predict_proba(X_test)[:, 1]
        
        # Calculate basic metrics
        metrics = {
            'model_name': model_name,
            'accuracy': accuracy_score(y_test, y_pred),
            'precision': precision_score(y_test, y_pred, average='binary'),
            'recall': recall_score(y_test, y_pred, average='binary'),
            'specificity': self._calculate_specificity(y_test, y_pred),
            'f1_score': f1_score(y_test, y_pred, average='binary'),
            'roc_auc': roc_auc_score(y_test, y_pred_proba),
            'average_precision': average_precision_score(y_test, y_pred_proba),
            'balanced_accuracy': self._calculate_balanced_accuracy(y_test, y_pred)
        }
        
        # Add confusion matrix metrics
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
        metrics.update({
            'true_positives': tp,
            'true_negatives': tn,
            'false_positives': fp,
            'false_negatives': fn,
            'positive_predictive_value': tp / (tp + fp) if (tp + fp) > 0 else 0,
            'negative_predictive_value': tn / (tn + fn) if (tn + fn) > 0 else 0
        })
        
        # Log metrics
        logger.info(f"{model_name} Performance:")
        logger.info(f"  Accuracy: {metrics['accuracy']:.4f}")
        logger.info(f"  Precision: {metrics['precision']:.4f}")
        logger.info(f"  Recall: {metrics['recall']:.4f}")
        logger.info(f"  F1-Score: {metrics['f1_score']:.4f}")
        logger.info(f"  ROC-AUC: {metrics['roc_auc']:.4f}")
        
        return metrics
    
    def cross_validate_model(self, model: Any, X: np.ndarray, y: np.ndarray,
                           cv_folds: int = 5, scoring: List[str] = None) -> Dict[str, np.ndarray]:
        """
        Perform cross-validation evaluation.
        
        Args:
            model: Model to evaluate
            X (np.ndarray): Features
            y (np.ndarray): Labels
            cv_folds (int): Number of CV folds
            scoring (List[str]): Scoring metrics
            
        Returns:
            Dict: Cross-validation scores for each metric
        """
        if scoring is None:
            scoring = ['accuracy', 'precision', 'recall', 'f1', 'roc_auc']
        
        logger.info(f"Performing {cv_folds}-fold cross-validation")
        
        cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
        
        cv_results = cross_validate(
            model, X, y, cv=cv, scoring=scoring, return_train_score=True
        )
        
        # Calculate statistics
        results = {}
        for metric in scoring:
            test_scores = cv_results[f'test_{metric}']
            train_scores = cv_results[f'train_{metric}']
            
            results[metric] = {
                'test_mean': np.mean(test_scores),
                'test_std': np.std(test_scores),
                'train_mean': np.mean(train_scores),
                'train_std': np.std(train_scores),
                'test_scores': test_scores,
                'train_scores': train_scores
            }
            
            logger.info(f"  {metric}: {results[metric]['test_mean']:.4f} ± {results[metric]['test_std']:.4f}")
        
        return results
    
    def compare_models(self, models: Dict[str, Any], X_test: np.ndarray, 
                      y_test: np.ndarray) -> pd.DataFrame:
        """
        Compare multiple models on the same test set.
        
        Args:
            models (Dict[str, Any]): Dictionary of model_name -> model
            X_test (np.ndarray): Test features
            y_test (np.ndarray): Test labels
            
        Returns:
            pd.DataFrame: Comparison results
        """
        logger.info(f"Comparing {len(models)} models")
        
        results = []
        for model_name, model in models.items():
            metrics = self.evaluate_single_model(model, X_test, y_test, model_name)
            results.append(metrics)
        
        comparison_df = pd.DataFrame(results)
        
        # Sort by ROC-AUC
        comparison_df = comparison_df.sort_values('roc_auc', ascending=False)
        
        # Save comparison
        output_file = self.output_dir / "model_comparison.csv"
        comparison_df.to_csv(output_file, index=False)
        logger.info(f"Model comparison saved to {output_file}")
        
        return comparison_df
    
    def plot_roc_curves(self, models: Dict[str, Any], X_test: np.ndarray, 
                       y_test: np.ndarray, save_plot: bool = True) -> plt.Figure:
        """
        Plot ROC curves for multiple models.
        
        Args:
            models (Dict[str, Any]): Dictionary of model_name -> model
            X_test (np.ndarray): Test features
            y_test (np.ndarray): Test labels
            save_plot (bool): Whether to save the plot
            
        Returns:
            plt.Figure: ROC curve plot
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        for model_name, model in models.items():
            y_pred_proba = model.predict_proba(X_test)[:, 1]
            fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
            auc = roc_auc_score(y_test, y_pred_proba)
            
            ax.plot(fpr, tpr, label=f'{model_name} (AUC = {auc:.3f})', linewidth=2)
        
        # Plot diagonal line
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
        
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title('ROC Curves - MGMT Methylation Prediction')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_plot:
            output_file = self.output_dir / "roc_curves.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"ROC curves saved to {output_file}")
        
        return fig
    
    def plot_precision_recall_curves(self, models: Dict[str, Any], X_test: np.ndarray,
                                    y_test: np.ndarray, save_plot: bool = True) -> plt.Figure:
        """
        Plot Precision-Recall curves for multiple models.
        
        Args:
            models (Dict[str, Any]): Dictionary of model_name -> model
            X_test (np.ndarray): Test features
            y_test (np.ndarray): Test labels
            save_plot (bool): Whether to save the plot
            
        Returns:
            plt.Figure: Precision-Recall curve plot
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        for model_name, model in models.items():
            y_pred_proba = model.predict_proba(X_test)[:, 1]
            precision, recall, _ = precision_recall_curve(y_test, y_pred_proba)
            avg_precision = average_precision_score(y_test, y_pred_proba)
            
            ax.plot(recall, precision, label=f'{model_name} (AP = {avg_precision:.3f})', linewidth=2)
        
        # Plot baseline
        baseline = np.sum(y_test) / len(y_test)
        ax.axhline(y=baseline, color='k', linestyle='--', alpha=0.5, 
                  label=f'Baseline (AP = {baseline:.3f})')
        
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.set_title('Precision-Recall Curves - MGMT Methylation Prediction')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_plot:
            output_file = self.output_dir / "precision_recall_curves.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Precision-Recall curves saved to {output_file}")
        
        return fig
    
    def plot_confusion_matrices(self, models: Dict[str, Any], X_test: np.ndarray,
                               y_test: np.ndarray, label_names: List[str] = None,
                               save_plot: bool = True) -> plt.Figure:
        """
        Plot confusion matrices for multiple models.
        
        Args:
            models (Dict[str, Any]): Dictionary of model_name -> model
            X_test (np.ndarray): Test features
            y_test (np.ndarray): Test labels
            label_names (List[str]): Class label names
            save_plot (bool): Whether to save the plot
            
        Returns:
            plt.Figure: Confusion matrices plot
        """
        if label_names is None:
            label_names = ['MGMT Unmethylated', 'MGMT Methylated']
        
        n_models = len(models)
        cols = min(3, n_models)
        rows = (n_models + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
        if n_models == 1:
            axes = [axes]
        elif rows == 1:
            axes = axes.reshape(1, -1)
        
        for idx, (model_name, model) in enumerate(models.items()):
            row = idx // cols
            col = idx % cols
            ax = axes[row, col] if rows > 1 else axes[col]
            
            y_pred = model.predict(X_test)
            cm = confusion_matrix(y_test, y_pred)
            
            # Normalize confusion matrix
            cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            
            sns.heatmap(cm_normalized, annot=True, fmt='.2f', cmap='Blues',
                       xticklabels=label_names, yticklabels=label_names, ax=ax)
            ax.set_title(f'{model_name}')
            ax.set_ylabel('True Label')
            ax.set_xlabel('Predicted Label')
        
        # Hide empty subplots
        for idx in range(n_models, rows * cols):
            if rows > 1:
                axes[idx // cols, idx % cols].set_visible(False)
            else:
                axes[idx].set_visible(False)
        
        plt.tight_layout()
        
        if save_plot:
            output_file = self.output_dir / "confusion_matrices.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Confusion matrices saved to {output_file}")
        
        return fig
    
    def plot_feature_importance(self, model: Any, feature_names: List[str],
                               top_n: int = 20, save_plot: bool = True) -> plt.Figure:
        """
        Plot feature importance for a model.
        
        Args:
            model: Trained model with feature importance
            feature_names (List[str]): Names of features
            top_n (int): Number of top features to plot
            save_plot (bool): Whether to save the plot
            
        Returns:
            plt.Figure: Feature importance plot
        """
        # Get feature importance
        if hasattr(model, 'feature_importances_'):
            importance = model.feature_importances_
        elif hasattr(model, 'coef_'):
            importance = np.abs(model.coef_[0])
        else:
            logger.warning("Model does not have feature importance")
            return None
        
        # Create importance DataFrame
        importance_df = pd.DataFrame({
            'feature': feature_names,
            'importance': importance
        }).sort_values('importance', ascending=False).head(top_n)
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.barplot(data=importance_df, y='feature', x='importance', ax=ax)
        ax.set_title(f'Top {top_n} Feature Importance')
        ax.set_xlabel('Importance')
        
        plt.tight_layout()
        
        if save_plot:
            output_file = self.output_dir / "feature_importance.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Feature importance plot saved to {output_file}")
        
        return fig
    
    def plot_learning_curves(self, model: Any, X: np.ndarray, y: np.ndarray,
                           train_sizes: np.ndarray = None, save_plot: bool = True) -> plt.Figure:
        """
        Plot learning curves to assess model performance vs training size.
        
        Args:
            model: Model to evaluate
            X (np.ndarray): Features
            y (np.ndarray): Labels
            train_sizes (np.ndarray): Training sizes to evaluate
            save_plot (bool): Whether to save the plot
            
        Returns:
            plt.Figure: Learning curves plot
        """
        if train_sizes is None:
            train_sizes = np.linspace(0.1, 1.0, 10)
        
        train_sizes, train_scores, val_scores = learning_curve(
            model, X, y, train_sizes=train_sizes, cv=5, scoring='roc_auc',
            random_state=42, n_jobs=-1
        )
        
        # Calculate means and standard deviations
        train_mean = np.mean(train_scores, axis=1)
        train_std = np.std(train_scores, axis=1)
        val_mean = np.mean(val_scores, axis=1)
        val_std = np.std(val_scores, axis=1)
        
        # Plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(train_sizes, train_mean, 'o-', label='Training Score', linewidth=2)
        ax.fill_between(train_sizes, train_mean - train_std, train_mean + train_std, alpha=0.1)
        
        ax.plot(train_sizes, val_mean, 'o-', label='Validation Score', linewidth=2)
        ax.fill_between(train_sizes, val_mean - val_std, val_mean + val_std, alpha=0.1)
        
        ax.set_xlabel('Training Set Size')
        ax.set_ylabel('ROC-AUC Score')
        ax.set_title('Learning Curves')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_plot:
            output_file = self.output_dir / "learning_curves.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Learning curves saved to {output_file}")
        
        return fig
    
    def calibration_analysis(self, models: Dict[str, Any], X_test: np.ndarray,
                           y_test: np.ndarray, save_plot: bool = True) -> plt.Figure:
        """
        Analyze model calibration (reliability of predicted probabilities).
        
        Args:
            models (Dict[str, Any]): Dictionary of model_name -> model
            X_test (np.ndarray): Test features
            y_test (np.ndarray): Test labels
            save_plot (bool): Whether to save the plot
            
        Returns:
            plt.Figure: Calibration plot
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        for model_name, model in models.items():
            y_pred_proba = model.predict_proba(X_test)[:, 1]
            
            # Calculate calibration curve
            fraction_of_positives, mean_predicted_value = calibration_curve(
                y_test, y_pred_proba, n_bins=10
            )
            
            ax.plot(mean_predicted_value, fraction_of_positives, 'o-', 
                   label=f'{model_name}', linewidth=2)
        
        # Perfect calibration line
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Perfect Calibration')
        
        ax.set_xlabel('Mean Predicted Probability')
        ax.set_ylabel('Fraction of Positives')
        ax.set_title('Calibration Plot - MGMT Methylation Prediction')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_plot:
            output_file = self.output_dir / "calibration_plot.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Calibration plot saved to {output_file}")
        
        return fig
    
    def generate_evaluation_report(self, models: Dict[str, Any], X_test: np.ndarray,
                                 y_test: np.ndarray, feature_names: List[str] = None) -> str:
        """
        Generate a comprehensive evaluation report.
        
        Args:
            models (Dict[str, Any]): Dictionary of model_name -> model
            X_test (np.ndarray): Test features
            y_test (np.ndarray): Test labels
            feature_names (List[str]): Names of features
            
        Returns:
            str: Path to the generated report
        """
        logger.info("Generating comprehensive evaluation report")
        
        # Model comparison
        comparison_df = self.compare_models(models, X_test, y_test)
        
        # Generate plots
        self.plot_roc_curves(models, X_test, y_test)
        self.plot_precision_recall_curves(models, X_test, y_test)
        self.plot_confusion_matrices(models, X_test, y_test)
        self.calibration_analysis(models, X_test, y_test)
        
        # Feature importance for best model
        best_model_name = comparison_df.iloc[0]['model_name']
        best_model = models[best_model_name]
        
        if feature_names:
            self.plot_feature_importance(best_model, feature_names)
        
        # Generate text report
        report_file = self.output_dir / "evaluation_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("MGMT Methylation Prediction - Model Evaluation Report\n")
            f.write("=" * 60 + "\n\n")
            
            f.write("Model Performance Summary:\n")
            f.write("-" * 30 + "\n")
            for _, row in comparison_df.iterrows():
                f.write(f"\n{row['model_name']}:\n")
                f.write(f"  ROC-AUC: {row['roc_auc']:.4f}\n")
                f.write(f"  Accuracy: {row['accuracy']:.4f}\n")
                f.write(f"  Precision: {row['precision']:.4f}\n")
                f.write(f"  Recall: {row['recall']:.4f}\n")
                f.write(f"  F1-Score: {row['f1_score']:.4f}\n")
            
            f.write(f"\nBest Model: {best_model_name}\n")
            f.write(f"Best ROC-AUC: {comparison_df.iloc[0]['roc_auc']:.4f}\n")
            
            f.write("\nGenerated Files:\n")
            f.write("- model_comparison.csv\n")
            f.write("- roc_curves.png\n")
            f.write("- precision_recall_curves.png\n")
            f.write("- confusion_matrices.png\n")
            f.write("- calibration_plot.png\n")
            if feature_names:
                f.write("- feature_importance.png\n")
        
        logger.info(f"Evaluation report saved to {report_file}")
        
        return str(report_file)
    
    def _calculate_specificity(self, y_true: np.ndarray, y_pred: np.ndarray) -> float:
        """Calculate specificity (true negative rate)."""
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        return tn / (tn + fp) if (tn + fp) > 0 else 0
    
    def _calculate_balanced_accuracy(self, y_true: np.ndarray, y_pred: np.ndarray) -> float:
        """Calculate balanced accuracy."""
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        return (sensitivity + specificity) / 2

if __name__ == "__main__":
    # Example usage would go here
    logger.info("Model evaluation module loaded successfully")

