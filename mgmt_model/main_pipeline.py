"""
Main pipeline for MGMT methylation detection using genome-wide methylation data.
Orchestrates all modules: data loading, feature engineering, model training, and evaluation.
"""

import pandas as pd
import numpy as np
import argparse
import json
from pathlib import Path
from typing import List, Dict, Optional
import logging
import warnings
warnings.filterwarnings('ignore')

# Import custom modules
from data_loader import BedMethylLoader, load_sample_metadata, create_sample_metadata_template
from feature_engineering import MethylationFeatureExtractor
from model_training import MGMTModelTrainer, train_mgmt_model_pipeline
from model_evaluation import ModelEvaluator

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MGMTPipeline:
    """
    Complete pipeline for MGMT methylation status prediction from genome-wide methylation data.
    """
    
    def __init__(self, config: Dict):
        """
        Initialize the pipeline with configuration.
        
        Args:
            config (Dict): Configuration parameters
        """
        self.config = config
        self.output_dir = Path(config.get('output_dir', 'mgmt_results'))
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize components
        self.data_loader = BedMethylLoader()
        self.feature_extractor = MethylationFeatureExtractor()
        self.model_trainer = MGMTModelTrainer(random_state=config.get('random_state', 42))
        self.evaluator = ModelEvaluator(output_dir=str(self.output_dir / 'evaluation'))
        
        logger.info(f"Pipeline initialized with output directory: {self.output_dir}")
    
    def load_data(self) -> Dict[str, pd.DataFrame]:
        """
        Load methylation data from bedmethyl files.
        
        Returns:
            Dict: Dictionary containing loaded data for each sample
        """
        logger.info("Starting data loading phase")
        
        # Get input files
        input_files = self.config['input_files']
        if isinstance(input_files, str):
            input_files = [input_files]
        
        # Load data
        if len(input_files) == 1:
            # Single file
            sample_data = self.data_loader.load_bedmethyl_file(
                input_files[0],
                chromosomes=self.config.get('chromosomes'),
                modification_types=self.config.get('modification_types', ['m']),
                min_coverage=self.config.get('min_coverage', 3)
            )
            
            # Add sample name
            sample_name = self.config.get('sample_name', Path(input_files[0]).stem.split('.')[0])
            sample_data['sample'] = sample_name
            
            samples_data = {sample_name: sample_data}
            
        else:
            # Multiple files
            sample_names = self.config.get('sample_names')
            samples_data = self.data_loader.load_multiple_samples(
                input_files,
                sample_names=sample_names,
                chromosomes=self.config.get('chromosomes'),
                modification_types=self.config.get('modification_types', ['m']),
                min_coverage=self.config.get('min_coverage', 3)
            )
        
        # Combine all samples into one DataFrame
        combined_data = pd.concat(samples_data.values(), ignore_index=True)
        
        # Validate data
        validation_results = self.data_loader.validate_data(combined_data)
        
        # Save validation results
        validation_file = self.output_dir / 'data_validation.json'
        with open(validation_file, 'w') as f:
            # Convert numpy types to native Python types for JSON serialization
            def convert_numpy(obj):
                if isinstance(obj, np.integer):
                    return int(obj)
                elif isinstance(obj, np.floating):
                    return float(obj)
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()
                return obj
            
            json.dump(validation_results, f, indent=2, default=convert_numpy)
        
        logger.info(f"Data validation results saved to {validation_file}")
        logger.info(f"Loaded {len(combined_data)} methylation sites from {len(samples_data)} samples")
        
        return {'combined': combined_data, 'samples': samples_data}
    
    def extract_features(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Extract features from methylation data.
        
        Args:
            data (pd.DataFrame): Combined methylation data
            
        Returns:
            pd.DataFrame: Feature matrix
        """
        logger.info("Starting feature extraction phase")
        
        # Set window sizes if specified
        if 'window_sizes' in self.config:
            self.feature_extractor.window_sizes = self.config['window_sizes']
        
        # Extract all features
        features_df = self.feature_extractor.combine_all_features(data)
        
        # Scale features if requested
        if self.config.get('scale_features', True):
            scaling_method = self.config.get('scaling_method', 'robust')
            scaled_features = self.feature_extractor.scale_features(
                features_df, method=scaling_method, fit_scaler=True
            )
        else:
            scaled_features = features_df
        
        # Save features
        features_file = self.output_dir / 'features.csv'
        scaled_features.to_csv(features_file, index=False)
        logger.info(f"Features saved to {features_file}")
        
        # Save feature names
        feature_names = [col for col in scaled_features.columns if col != 'sample_id']
        feature_names_file = self.output_dir / 'feature_names.txt'
        with open(feature_names_file, 'w') as f:
            for name in feature_names:
                f.write(f"{name}\n")
        
        logger.info(f"Generated {len(feature_names)} features")
        
        return scaled_features
    
    def prepare_labels(self, features_df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare or load sample labels.
        
        Args:
            features_df (pd.DataFrame): Feature matrix
            
        Returns:
            pd.DataFrame: Sample metadata with labels
        """
        logger.info("Preparing sample labels")
        
        metadata_file = self.config.get('metadata_file')
        
        if metadata_file and Path(metadata_file).exists():
            # Load existing metadata
            metadata_df = load_sample_metadata(metadata_file)
            logger.info(f"Loaded metadata from {metadata_file}")
            
        else:
            # Create template metadata file
            sample_ids = features_df['sample_id'].tolist()
            template_file = self.output_dir / 'sample_metadata_template.csv'
            
            create_sample_metadata_template(str(template_file), sample_ids)
            
            logger.warning(f"No metadata file provided. Created template: {template_file}")
            logger.warning("Please fill in the MGMT status labels and re-run the pipeline.")
            
            # Create dummy metadata for demonstration
            metadata_df = pd.DataFrame({
                'sample_id': sample_ids,
                'mgmt_status': ['unknown'] * len(sample_ids)
            })
        
        return metadata_df
    
    def train_models(self, features_df: pd.DataFrame, metadata_df: pd.DataFrame) -> MGMTModelTrainer:
        """
        Train machine learning models.
        
        Args:
            features_df (pd.DataFrame): Feature matrix
            metadata_df (pd.DataFrame): Sample metadata with labels
            
        Returns:
            MGMTModelTrainer: Trained model trainer
        """
        logger.info("Starting model training phase")
        
        # Check if we have valid labels
        if metadata_df['mgmt_status'].str.contains('unknown').any():
            logger.warning("Some samples have unknown MGMT status. Cannot train models.")
            logger.warning("Please update the metadata file with correct labels.")
            return None
        
        # Training configuration
        training_config = {
            'feature_selection': self.config.get('feature_selection', True),
            'n_features': self.config.get('n_features', 50),
            'hyperparameter_tuning': self.config.get('hyperparameter_tuning', True),
            'create_ensemble': self.config.get('create_ensemble', True)
        }
        
        # Train models using the pipeline function
        trainer = train_mgmt_model_pipeline(
            features_df=features_df,
            metadata_df=metadata_df,
            output_dir=str(self.output_dir / 'models'),
            **training_config
        )
        
        logger.info("Model training completed")
        
        return trainer
    
    def evaluate_models(self, trainer: MGMTModelTrainer, features_df: pd.DataFrame,
                       metadata_df: pd.DataFrame) -> str:
        """
        Evaluate trained models.
        
        Args:
            trainer (MGMTModelTrainer): Trained model trainer
            features_df (pd.DataFrame): Feature matrix
            metadata_df (pd.DataFrame): Sample metadata with labels
            
        Returns:
            str: Path to evaluation report
        """
        logger.info("Starting model evaluation phase")
        
        # Prepare test data (using the same data for demonstration)
        data = features_df.merge(metadata_df, on='sample_id', how='inner')
        feature_columns = [col for col in features_df.columns if col != 'sample_id']
        
        X = data[feature_columns].values
        y = trainer.label_encoder.transform(data['mgmt_status'].values)
        
        # Apply feature selection if used
        if trainer.feature_selector is not None:
            X = trainer.feature_selector.transform(X)
            feature_columns = [feature_columns[i] for i in trainer.feature_selector.get_support(indices=True)]
        
        # Generate comprehensive evaluation report
        report_path = self.evaluator.generate_evaluation_report(
            models=trainer.models,
            X_test=X,
            y_test=y,
            feature_names=feature_columns
        )
        
        logger.info("Model evaluation completed")
        
        return report_path
    
    def predict_samples(self, trainer: MGMTModelTrainer, features_df: pd.DataFrame,
                       sample_ids: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Make predictions on new samples.
        
        Args:
            trainer (MGMTModelTrainer): Trained model trainer
            features_df (pd.DataFrame): Feature matrix
            sample_ids (List[str], optional): Specific samples to predict
            
        Returns:
            pd.DataFrame: Predictions with probabilities
        """
        logger.info("Making predictions on samples")
        
        if trainer.best_model is None:
            logger.error("No trained model available for predictions")
            return pd.DataFrame()
        
        # Filter samples if specified
        if sample_ids:
            features_df = features_df[features_df['sample_id'].isin(sample_ids)]
        
        # Prepare features
        feature_columns = [col for col in features_df.columns if col != 'sample_id']
        X = features_df[feature_columns].values
        
        # Apply feature selection if used
        if trainer.feature_selector is not None:
            X = trainer.feature_selector.transform(X)
        
        # Make predictions
        predictions = trainer.best_model.predict(X)
        probabilities = trainer.best_model.predict_proba(X)
        
        # Create results DataFrame
        results_df = pd.DataFrame({
            'sample_id': features_df['sample_id'].values,
            'predicted_mgmt_status': trainer.label_encoder.inverse_transform(predictions),
            'probability_methylated': probabilities[:, 1],
            'probability_unmethylated': probabilities[:, 0],
            'confidence': np.max(probabilities, axis=1)
        })
        
        # Save predictions
        predictions_file = self.output_dir / 'predictions.csv'
        results_df.to_csv(predictions_file, index=False)
        logger.info(f"Predictions saved to {predictions_file}")
        
        return results_df
    
    def run_complete_pipeline(self) -> Dict:
        """
        Run the complete MGMT prediction pipeline.
        
        Returns:
            Dict: Pipeline results and outputs
        """
        logger.info("Starting complete MGMT methylation prediction pipeline")
        
        results = {}
        
        try:
            # 1. Load data
            data_dict = self.load_data()
            results['data_loaded'] = True
            results['n_samples'] = len(data_dict['combined']['sample'].unique())
            results['n_methylation_sites'] = len(data_dict['combined'])
            
            # 2. Extract features
            features_df = self.extract_features(data_dict['combined'])
            results['features_extracted'] = True
            results['n_features'] = len([col for col in features_df.columns if col != 'sample_id'])
            
            # 3. Prepare labels
            metadata_df = self.prepare_labels(features_df)
            results['labels_prepared'] = True
            
            # 4. Train models (if labels are available)
            if not metadata_df['mgmt_status'].str.contains('unknown').any():
                trainer = self.train_models(features_df, metadata_df)
                if trainer:
                    results['models_trained'] = True
                    results['best_model'] = type(trainer.best_model).__name__
                    
                    # 5. Evaluate models
                    report_path = self.evaluate_models(trainer, features_df, metadata_df)
                    results['models_evaluated'] = True
                    results['evaluation_report'] = report_path
                    
                    # 6. Make predictions
                    predictions_df = self.predict_samples(trainer, features_df)
                    results['predictions_made'] = True
                    results['predictions_file'] = str(self.output_dir / 'predictions.csv')
                    
                else:
                    results['models_trained'] = False
                    logger.error("Model training failed")
            else:
                results['models_trained'] = False
                results['reason'] = "No valid labels provided"
                
                # Make dummy predictions to show pipeline structure
                logger.info("Creating dummy predictions for pipeline demonstration")
                dummy_predictions = pd.DataFrame({
                    'sample_id': features_df['sample_id'],
                    'predicted_mgmt_status': ['unknown'] * len(features_df),
                    'probability_methylated': [0.5] * len(features_df),
                    'probability_unmethylated': [0.5] * len(features_df),
                    'confidence': [0.5] * len(features_df)
                })
                
                predictions_file = self.output_dir / 'dummy_predictions.csv'
                dummy_predictions.to_csv(predictions_file, index=False)
                results['dummy_predictions'] = str(predictions_file)
            
            # Save pipeline summary
            summary_file = self.output_dir / 'pipeline_summary.json'
            with open(summary_file, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            logger.info("Pipeline completed successfully")
            logger.info(f"Results summary saved to {summary_file}")
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            results['error'] = str(e)
            results['success'] = False
            
        return results

def load_config(config_file: str) -> Dict:
    """
    Load configuration from JSON file.
    
    Args:
        config_file (str): Path to configuration file
        
    Returns:
        Dict: Configuration parameters
    """
    with open(config_file, 'r') as f:
        config = json.load(f)
    return config

def create_default_config(output_file: str, input_file: str) -> None:
    """
    Create a default configuration file.
    
    Args:
        output_file (str): Path for configuration file
        input_file (str): Input bedmethyl file path
    """
    default_config = {
        "input_files": [input_file],
        "sample_names": None,
        "metadata_file": None,
        "output_dir": "mgmt_results",
        "chromosomes": None,
        "modification_types": ["m"],
        "min_coverage": 3,
        "window_sizes": [1000, 5000, 10000, 50000],
        "scale_features": True,
        "scaling_method": "robust",
        "feature_selection": True,
        "n_features": 50,
        "hyperparameter_tuning": True,
        "create_ensemble": True,
        "random_state": 42
    }
    
    with open(output_file, 'w') as f:
        json.dump(default_config, f, indent=2)
    
    logger.info(f"Default configuration created: {output_file}")

def main():
    """
    Main function for command-line interface.
    """
    parser = argparse.ArgumentParser(
        description="MGMT Methylation Prediction Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with single file (create config first)
  python main_pipeline.py --input /path/to/sample.bedmethyl.gz --create-config

  # Run with configuration file
  python main_pipeline.py --config config.json

  # Run with multiple samples
  python main_pipeline.py --input file1.gz file2.gz --sample-names S1 S2
        """
    )
    
    parser.add_argument('--input', nargs='+', 
                       help='Input bedmethyl file(s)')
    parser.add_argument('--config', 
                       help='Configuration file (JSON)')
    parser.add_argument('--sample-names', nargs='+',
                       help='Sample names (if multiple files)')
    parser.add_argument('--metadata', 
                       help='Sample metadata file with MGMT labels')
    parser.add_argument('--output-dir', default='mgmt_results',
                       help='Output directory')
    parser.add_argument('--create-config', action='store_true',
                       help='Create default configuration file')
    parser.add_argument('--min-coverage', type=int, default=3,
                       help='Minimum coverage threshold')
    
    args = parser.parse_args()
    
    if args.create_config:
        if not args.input:
            logger.error("Input file required for creating config")
            return
        
        config_file = Path(args.output_dir) / 'config.json'
        create_default_config(str(config_file), args.input[0])
        logger.info(f"Edit {config_file} and run again with --config {config_file}")
        return
    
    # Load or create configuration
    if args.config:
        config = load_config(args.config)
    else:
        if not args.input:
            logger.error("Either --config or --input is required")
            return
        
        # Create config from command line arguments
        config = {
            'input_files': args.input,
            'sample_names': args.sample_names,
            'metadata_file': args.metadata,
            'output_dir': args.output_dir,
            'modification_types': ['m'],
            'min_coverage': args.min_coverage,
            'scale_features': True,
            'feature_selection': True,
            'n_features': 50,
            'hyperparameter_tuning': True,
            'create_ensemble': True,
            'random_state': 42
        }
    
    # Run pipeline
    pipeline = MGMTPipeline(config)
    results = pipeline.run_complete_pipeline()
    
    # Print summary
    print("\n" + "="*60)
    print("MGMT METHYLATION PREDICTION PIPELINE RESULTS")
    print("="*60)
    
    if results.get('data_loaded'):
        print(f"✓ Data loaded: {results['n_samples']} samples, {results['n_methylation_sites']} sites")
    
    if results.get('features_extracted'):
        print(f"✓ Features extracted: {results['n_features']} features")
    
    if results.get('models_trained'):
        print(f"✓ Models trained: Best model = {results.get('best_model', 'Unknown')}")
    else:
        print(f"✗ Models not trained: {results.get('reason', 'Unknown error')}")
    
    if results.get('models_evaluated'):
        print(f"✓ Models evaluated: Report = {results['evaluation_report']}")
    
    if results.get('predictions_made'):
        print(f"✓ Predictions made: {results['predictions_file']}")
    elif results.get('dummy_predictions'):
        print(f"⚠ Dummy predictions created: {results['dummy_predictions']}")
    
    print(f"\nAll results saved to: {config['output_dir']}")

if __name__ == "__main__":
    main()
