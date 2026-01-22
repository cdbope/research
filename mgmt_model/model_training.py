"""
Machine learning model training module for MGMT methylation detection.
Supports multiple algorithms and provides comprehensive training pipeline.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from sklearn.model_selection import (
    train_test_split, cross_val_score, GridSearchCV, 
    StratifiedKFold, RandomizedSearchCV
)
from sklearn.ensemble import (
    RandomForestClassifier, GradientBoostingClassifier, 
    ExtraTreesClassifier, VotingClassifier
)
from sklearn.linear_model import LogisticRegression, ElasticNet
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score, 
    roc_auc_score, confusion_matrix, classification_report,
    roc_curve, precision_recall_curve
)
from sklearn.feature_selection import (
    SelectKBest, f_classif, RFE, SelectFromModel
)
from sklearn.preprocessing import LabelEncoder
import joblib
import logging
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MGMTModelTrainer:
    """
    Comprehensive model training pipeline for MGMT methylation status prediction.
    """
    
    def __init__(self, random_state: int = 42):
        """
        Initialize the model trainer.
        
        Args:
            random_state (int): Random state for reproducibility
        """
        self.random_state = random_state
        self.models = {}
        self.best_model = None
        self.feature_selector = None
        self.label_encoder = LabelEncoder()
        self.feature_importance = {}
        
        # Define model configurations
        self.model_configs = {
            'random_forest': {
                'model': RandomForestClassifier(random_state=random_state),
                'params': {
                    'n_estimators': [100, 200, 300, 500],
                    'max_depth': [None, 10, 20, 30],
                    'min_samples_split': [2, 5, 10],
                    'min_samples_leaf': [1, 2, 4],
                    'max_features': ['sqrt', 'log2', None]
                }
            },
            'gradient_boosting': {
                'model': GradientBoostingClassifier(random_state=random_state),
                'params': {
                    'n_estimators': [100, 200, 300],
                    'learning_rate': [0.01, 0.1, 0.2],
                    'max_depth': [3, 5, 7, 10],
                    'min_samples_split': [2, 5, 10],
                    'min_samples_leaf': [1, 2, 4]
                }
            },
            'extra_trees': {
                'model': ExtraTreesClassifier(random_state=random_state),
                'params': {
                    'n_estimators': [100, 200, 300],
                    'max_depth': [None, 10, 20, 30],
                    'min_samples_split': [2, 5, 10],
                    'min_samples_leaf': [1, 2, 4]
                }
            },
            'logistic_regression': {
                'model': LogisticRegression(random_state=random_state, max_iter=1000),
                'params': {
                    'C': [0.001, 0.01, 0.1, 1, 10, 100],
                    'penalty': ['l1', 'l2', 'elasticnet'],
                    'solver': ['liblinear', 'saga'],
                    'l1_ratio': [0.1, 0.3, 0.5, 0.7, 0.9]
                }
            },
            'svm': {
                'model': SVC(random_state=random_state, probability=True),
                'params': {
                    'C': [0.1, 1, 10, 100],
                    'gamma': ['scale', 'auto', 0.001, 0.01, 0.1, 1],
                    'kernel': ['rbf', 'linear', 'poly']
                }
            },
            'mlp': {
                'model': MLPClassifier(random_state=random_state, max_iter=500),
                'params': {
                    'hidden_layer_sizes': [(50,), (100,), (100, 50), (200, 100)],
                    'activation': ['relu', 'tanh'],
                    'alpha': [0.0001, 0.001, 0.01],
                    'learning_rate': ['constant', 'adaptive']
                }
            },
            'naive_bayes': {
                'model': GaussianNB(),
                'params': {
                    'var_smoothing': [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]
                }
            }
        }
    
    def prepare_data(self, features_df: pd.DataFrame, 
                    metadata_df: pd.DataFrame,
                    target_column: str = 'mgmt_status',
                    test_size: float = 0.2) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[str]]:
        """
        Prepare data for training.
        
        Args:
            features_df (pd.DataFrame): Feature matrix
            metadata_df (pd.DataFrame): Sample metadata with labels
            target_column (str): Target column name in metadata
            test_size (float): Test set proportion
            
        Returns:
            Tuple: X_train, X_test, y_train, y_test, feature_names
        """
        logger.info("Preparing data for training")
        
        # Merge features with metadata
        data = features_df.merge(metadata_df, left_on='sample_id', right_on='sample_id', how='inner')
        
        if len(data) == 0:
            raise ValueError("No matching samples found between features and metadata")
        
        # Prepare features and target
        feature_columns = [col for col in features_df.columns if col != 'sample_id']
        X = data[feature_columns].values
        y = data[target_column].values
        
        # Encode labels
        y_encoded = self.label_encoder.fit_transform(y)
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y_encoded, test_size=test_size, random_state=self.random_state, 
            stratify=y_encoded
        )
        
        logger.info(f"Training set: {X_train.shape[0]} samples, {X_train.shape[1]} features")
        logger.info(f"Test set: {X_test.shape[0]} samples")
        logger.info(f"Class distribution: {dict(zip(*np.unique(y_encoded, return_counts=True)))}")
        
        return X_train, X_test, y_train, y_test, feature_columns
    
    def select_features(self, X_train: np.ndarray, y_train: np.ndarray, 
                       feature_names: List[str],
                       method: str = 'rfe',
                       n_features: int = 50) -> np.ndarray:
        """
        Perform feature selection.
        
        Args:
            X_train (np.ndarray): Training features
            y_train (np.ndarray): Training labels
            feature_names (List[str]): Feature names
            method (str): Feature selection method ('rfe', 'univariate', 'model_based')
            n_features (int): Number of features to select
            
        Returns:
            np.ndarray: Indices of selected features
        """
        logger.info(f"Selecting {n_features} features using {method} method")
        
        if method == 'rfe':
            # Recursive Feature Elimination
            estimator = RandomForestClassifier(n_estimators=100, random_state=self.random_state)
            self.feature_selector = RFE(estimator, n_features_to_select=n_features)
            
        elif method == 'univariate':
            # Univariate statistical test
            self.feature_selector = SelectKBest(score_func=f_classif, k=n_features)
            
        elif method == 'model_based':
            # Model-based feature selection
            estimator = RandomForestClassifier(n_estimators=100, random_state=self.random_state)
            estimator.fit(X_train, y_train)
            self.feature_selector = SelectFromModel(estimator, max_features=n_features)
            
        else:
            raise ValueError("Method must be 'rfe', 'univariate', or 'model_based'")
        
        # Fit selector and get selected features
        X_train_selected = self.feature_selector.fit_transform(X_train, y_train)
        selected_indices = self.feature_selector.get_support(indices=True)
        
        logger.info(f"Selected features: {[feature_names[i] for i in selected_indices[:10]]}...")
        
        return selected_indices
    
    def train_single_model(self, X_train: np.ndarray, y_train: np.ndarray,
                          model_name: str,
                          hyperparameter_tuning: bool = True,
                          cv_folds: int = 5) -> Dict[str, Any]:
        """
        Train a single model with optional hyperparameter tuning.
        
        Args:
            X_train (np.ndarray): Training features
            y_train (np.ndarray): Training labels
            model_name (str): Name of the model to train
            hyperparameter_tuning (bool): Whether to perform hyperparameter tuning
            cv_folds (int): Number of cross-validation folds
            
        Returns:
            Dict: Training results including model and scores
        """
        logger.info(f"Training {model_name} model")
        
        if model_name not in self.model_configs:
            raise ValueError(f"Unknown model: {model_name}")
        
        config = self.model_configs[model_name]
        model = config['model']
        
        if hyperparameter_tuning and len(config['params']) > 0:
            # Hyperparameter tuning with cross-validation
            cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=self.random_state)
            
            # Use RandomizedSearchCV for efficiency
            search = RandomizedSearchCV(
                model, config['params'], 
                n_iter=20,  # Reduced for efficiency
                cv=cv, 
                scoring='roc_auc',
                random_state=self.random_state,
                n_jobs=-1
            )
            
            search.fit(X_train, y_train)
            best_model = search.best_estimator_
            best_params = search.best_params_
            cv_score = search.best_score_
            
        else:
            # Train with default parameters
            model.fit(X_train, y_train)
            best_model = model
            best_params = {}
            
            # Cross-validation score
            cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=self.random_state)
            cv_scores = cross_val_score(model, X_train, y_train, cv=cv, scoring='roc_auc')
            cv_score = cv_scores.mean()
        
        # Store model
        self.models[model_name] = best_model
        
        # Get feature importance if available
        if hasattr(best_model, 'feature_importances_'):
            self.feature_importance[model_name] = best_model.feature_importances_
        elif hasattr(best_model, 'coef_'):
            self.feature_importance[model_name] = np.abs(best_model.coef_[0])
        
        results = {
            'model': best_model,
            'cv_score': cv_score,
            'best_params': best_params,
            'model_name': model_name
        }
        
        logger.info(f"{model_name} - CV AUC: {cv_score:.4f}")
        
        return results
    
    def train_all_models(self, X_train: np.ndarray, y_train: np.ndarray,
                        model_subset: Optional[List[str]] = None,
                        hyperparameter_tuning: bool = True) -> Dict[str, Dict]:
        """
        Train all available models.
        
        Args:
            X_train (np.ndarray): Training features
            y_train (np.ndarray): Training labels
            model_subset (List[str], optional): Subset of models to train
            hyperparameter_tuning (bool): Whether to perform hyperparameter tuning
            
        Returns:
            Dict: Results for all trained models
        """
        logger.info("Training all models")
        
        models_to_train = model_subset if model_subset else list(self.model_configs.keys())
        all_results = {}
        
        for model_name in models_to_train:
            try:
                results = self.train_single_model(
                    X_train, y_train, model_name, hyperparameter_tuning
                )
                all_results[model_name] = results
            except Exception as e:
                logger.error(f"Error training {model_name}: {e}")
                continue
        
        # Select best model based on CV score
        if all_results:
            best_model_name = max(all_results.keys(), key=lambda x: all_results[x]['cv_score'])
            self.best_model = all_results[best_model_name]['model']
            logger.info(f"Best model: {best_model_name} (CV AUC: {all_results[best_model_name]['cv_score']:.4f})")
        
        return all_results
    
    def create_ensemble_model(self, X_train: np.ndarray, y_train: np.ndarray,
                             model_names: Optional[List[str]] = None) -> VotingClassifier:
        """
        Create an ensemble model from trained models.
        
        Args:
            X_train (np.ndarray): Training features
            y_train (np.ndarray): Training labels
            model_names (List[str], optional): Names of models to include in ensemble
            
        Returns:
            VotingClassifier: Trained ensemble model
        """
        logger.info("Creating ensemble model")
        
        if not self.models:
            raise ValueError("No models have been trained yet")
        
        # Select models for ensemble
        if model_names is None:
            model_names = list(self.models.keys())
        
        ensemble_models = [(name, self.models[name]) for name in model_names if name in self.models]
        
        if len(ensemble_models) < 2:
            raise ValueError("Need at least 2 models for ensemble")
        
        # Create voting classifier
        ensemble = VotingClassifier(
            estimators=ensemble_models,
            voting='soft'  # Use predicted probabilities
        )
        
        # Train ensemble
        ensemble.fit(X_train, y_train)
        
        # Store ensemble as best model
        self.models['ensemble'] = ensemble
        self.best_model = ensemble
        
        logger.info(f"Ensemble created with {len(ensemble_models)} models")
        
        return ensemble
    
    def evaluate_model(self, model: Any, X_test: np.ndarray, y_test: np.ndarray,
                      model_name: str = "model") -> Dict[str, float]:
        """
        Evaluate a trained model on test data.
        
        Args:
            model: Trained model
            X_test (np.ndarray): Test features
            y_test (np.ndarray): Test labels
            model_name (str): Name of the model
            
        Returns:
            Dict: Evaluation metrics
        """
        # Predictions
        y_pred = model.predict(X_test)
        y_pred_proba = model.predict_proba(X_test)[:, 1]
        
        # Calculate metrics
        metrics = {
            'accuracy': accuracy_score(y_test, y_pred),
            'precision': precision_score(y_test, y_pred, average='weighted'),
            'recall': recall_score(y_test, y_pred, average='weighted'),
            'f1_score': f1_score(y_test, y_pred, average='weighted'),
            'roc_auc': roc_auc_score(y_test, y_pred_proba)
        }
        
        logger.info(f"{model_name} test metrics:")
        for metric, value in metrics.items():
            logger.info(f"  {metric}: {value:.4f}")
        
        return metrics
    
    def get_feature_importance_ranking(self, feature_names: List[str], 
                                     top_n: int = 20) -> pd.DataFrame:
        """
        Get feature importance ranking across all models.
        
        Args:
            feature_names (List[str]): Names of features
            top_n (int): Number of top features to return
            
        Returns:
            pd.DataFrame: Feature importance ranking
        """
        if not self.feature_importance:
            return pd.DataFrame()
        
        # Combine feature importance from all models
        importance_df = pd.DataFrame(index=feature_names)
        
        for model_name, importance in self.feature_importance.items():
            if len(importance) == len(feature_names):
                importance_df[model_name] = importance
        
        # Calculate mean importance
        importance_df['mean_importance'] = importance_df.mean(axis=1)
        importance_df = importance_df.sort_values('mean_importance', ascending=False)
        
        return importance_df.head(top_n)
    
    def save_models(self, output_dir: str) -> None:
        """
        Save all trained models to disk.
        
        Args:
            output_dir (str): Output directory for saved models
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save individual models
        for model_name, model in self.models.items():
            model_file = output_path / f"{model_name}_model.joblib"
            joblib.dump(model, model_file)
            logger.info(f"Saved {model_name} model to {model_file}")
        
        # Save feature selector
        if self.feature_selector is not None:
            selector_file = output_path / "feature_selector.joblib"
            joblib.dump(self.feature_selector, selector_file)
            logger.info(f"Saved feature selector to {selector_file}")
        
        # Save label encoder
        encoder_file = output_path / "label_encoder.joblib"
        joblib.dump(self.label_encoder, encoder_file)
        logger.info(f"Saved label encoder to {encoder_file}")
    
    def load_models(self, model_dir: str) -> None:
        """
        Load saved models from disk.
        
        Args:
            model_dir (str): Directory containing saved models
        """
        model_path = Path(model_dir)
        
        # Load individual models
        for model_file in model_path.glob("*_model.joblib"):
            model_name = model_file.stem.replace("_model", "")
            model = joblib.load(model_file)
            self.models[model_name] = model
            logger.info(f"Loaded {model_name} model")
        
        # Load feature selector
        selector_file = model_path / "feature_selector.joblib"
        if selector_file.exists():
            self.feature_selector = joblib.load(selector_file)
            logger.info("Loaded feature selector")
        
        # Load label encoder
        encoder_file = model_path / "label_encoder.joblib"
        if encoder_file.exists():
            self.label_encoder = joblib.load(encoder_file)
            logger.info("Loaded label encoder")
        
        # Set best model (assume ensemble if available, otherwise use random forest)
        if 'ensemble' in self.models:
            self.best_model = self.models['ensemble']
        elif 'random_forest' in self.models:
            self.best_model = self.models['random_forest']
        elif self.models:
            self.best_model = next(iter(self.models.values()))

def train_mgmt_model_pipeline(features_df: pd.DataFrame,
                             metadata_df: pd.DataFrame,
                             output_dir: str = "mgmt_models",
                             feature_selection: bool = True,
                             n_features: int = 50,
                             hyperparameter_tuning: bool = True,
                             create_ensemble: bool = True) -> MGMTModelTrainer:
    """
    Complete pipeline for training MGMT methylation models.
    
    Args:
        features_df (pd.DataFrame): Feature matrix
        metadata_df (pd.DataFrame): Sample metadata with MGMT labels
        output_dir (str): Output directory for models
        feature_selection (bool): Whether to perform feature selection
        n_features (int): Number of features to select
        hyperparameter_tuning (bool): Whether to tune hyperparameters
        create_ensemble (bool): Whether to create ensemble model
        
    Returns:
        MGMTModelTrainer: Trained model trainer
    """
    logger.info("Starting MGMT model training pipeline")
    
    # Initialize trainer
    trainer = MGMTModelTrainer()
    
    # Prepare data
    X_train, X_test, y_train, y_test, feature_names = trainer.prepare_data(
        features_df, metadata_df
    )
    
    # Feature selection
    if feature_selection:
        selected_indices = trainer.select_features(
            X_train, y_train, feature_names, n_features=n_features
        )
        X_train = X_train[:, selected_indices]
        X_test = X_test[:, selected_indices]
        feature_names = [feature_names[i] for i in selected_indices]
    
    # Train all models
    results = trainer.train_all_models(
        X_train, y_train, hyperparameter_tuning=hyperparameter_tuning
    )
    
    # Create ensemble
    if create_ensemble and len(results) > 1:
        ensemble = trainer.create_ensemble_model(X_train, y_train)
    
    # Evaluate best model
    if trainer.best_model is not None:
        metrics = trainer.evaluate_model(trainer.best_model, X_test, y_test, "best_model")
    
    # Get feature importance
    importance_df = trainer.get_feature_importance_ranking(feature_names)
    if not importance_df.empty:
        logger.info("Top 10 important features:")
        logger.info(importance_df.head(10)['mean_importance'].to_string())
    
    # Save models
    trainer.save_models(output_dir)
    
    logger.info("MGMT model training pipeline completed")
    
    return trainer

if __name__ == "__main__":
    # Example usage would go here
    logger.info("MGMT model training module loaded successfully")
