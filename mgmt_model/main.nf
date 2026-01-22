#!/usr/bin/env nextflow

/*
 * MGMT Methylation Detection Pipeline
 * 
 * A Nextflow pipeline for detecting MGMT methylation status using 
 * genome-wide methylation data from Oxford Nanopore modkit bedmethyl files.
 */

nextflow.enable.dsl = 2

// Define parameters
params.input = null
params.metadata = null
params.outdir = "results"
params.min_coverage = 3
params.modification_types = "m"
params.chromosomes = null
params.feature_selection = true
params.n_features = 50
params.hyperparameter_tuning = true
params.create_ensemble = true
params.scaling_method = "robust"
params.window_sizes = "1000,5000,10000,50000"
params.help = false

// Help message
def helpMessage() {
    log.info"""
    =========================================
    MGMT Methylation Detection Pipeline
    =========================================
    
    Usage:
      nextflow run main.nf --input samples.csv --metadata metadata.csv --outdir results
    
    Required arguments:
      --input               CSV file with sample information (sample_id,bedmethyl_path) or
                           single bedmethyl file path
      --metadata           CSV file with sample metadata including MGMT status
    
    Optional arguments:
      --outdir             Output directory (default: results)
      --min_coverage       Minimum coverage threshold (default: 3)
      --modification_types Modification types to include (default: m)
      --chromosomes        Comma-separated list of chromosomes to include (default: all)
      --feature_selection  Enable feature selection (default: true)
      --n_features         Number of features to select (default: 50)
      --hyperparameter_tuning Enable hyperparameter tuning (default: true)
      --create_ensemble    Create ensemble model (default: true)
      --scaling_method     Feature scaling method: robust, standard (default: robust)
      --window_sizes       Window sizes for regional analysis (default: 1000,5000,10000,50000)
    
    Example input CSV format:
      sample_id,bedmethyl_path
      T001,/path/to/T001.wf_mods.bedmethyl.gz
      T002,/path/to/T002.wf_mods.bedmethyl.gz
    
    Example metadata CSV format:
      sample_id,mgmt_status,batch,tissue_type
      T001,methylated,batch1,tumor
      T002,unmethylated,batch1,tumor
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.input) {
    log.error "Error: --input parameter is required"
    helpMessage()
    exit 1
}

// Define workflows
workflow {
    
    // Parse input
    if (params.input.endsWith('.csv')) {
        // Input is a CSV file with sample information
        input_ch = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> 
                tuple(row.sample_id, file(row.bedmethyl_path))
            }
    } else {
        // Input is a single bedmethyl file
        sample_id = file(params.input).getBaseName().split('\\.')[0]
        input_ch = Channel.of(tuple(sample_id, file(params.input)))
    }
    
    // Load and validate data
    LOAD_DATA(input_ch)
    
    // Collect all loaded data
    loaded_data = LOAD_DATA.out.sample_data.collect()
    
    // Extract features from all samples
    EXTRACT_FEATURES(loaded_data)
    
    // Prepare metadata
    if (params.metadata) {
        metadata_ch = Channel.fromPath(params.metadata)
        PREPARE_METADATA(metadata_ch, EXTRACT_FEATURES.out.features)
    } else {
        // Create template metadata
        CREATE_METADATA_TEMPLATE(EXTRACT_FEATURES.out.features)
        metadata_ch = CREATE_METADATA_TEMPLATE.out.template
        PREPARE_METADATA(metadata_ch, EXTRACT_FEATURES.out.features)
    }
    
    // Train models (only if valid labels are available)
    TRAIN_MODELS(
        EXTRACT_FEATURES.out.features,
        PREPARE_METADATA.out.metadata,
        EXTRACT_FEATURES.out.feature_names
    )
    
    // Evaluate models
    EVALUATE_MODELS(
        TRAIN_MODELS.out.models,
        EXTRACT_FEATURES.out.features,
        PREPARE_METADATA.out.metadata,
        EXTRACT_FEATURES.out.feature_names
    )
    
    // Make predictions
    PREDICT_SAMPLES(
        TRAIN_MODELS.out.best_model,
        EXTRACT_FEATURES.out.features,
        EXTRACT_FEATURES.out.feature_names
    )
    
    // Generate final report
    GENERATE_REPORT(
        LOAD_DATA.out.validation.collect(),
        EXTRACT_FEATURES.out.features,
        TRAIN_MODELS.out.training_results,
        EVALUATE_MODELS.out.evaluation_report,
        PREDICT_SAMPLES.out.predictions
    )
}

// Process definitions
process LOAD_DATA {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(bedmethyl_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_data.csv"), emit: sample_data
    path("${sample_id}_validation.json"), emit: validation
    
    script:
    chromosomes_arg = params.chromosomes ? "--chromosomes ${params.chromosomes}" : ""
    """
    #!/usr/bin/env python3
    
    import sys
    sys.path.append('${projectDir}')
    
    from data_loader import BedMethylLoader
    import json
    import pandas as pd
    
    # Load data
    loader = BedMethylLoader()
    
    # Parse chromosomes if provided
    chromosomes = None
    if "${params.chromosomes}" != "null" and "${params.chromosomes}" != "":
        chromosomes = "${params.chromosomes}".split(',')
    
    # Load bedmethyl file
    data = loader.load_bedmethyl_file(
        "${bedmethyl_file}",
        chromosomes=chromosomes,
        modification_types=["${params.modification_types}"],
        min_coverage=${params.min_coverage}
    )
    
    # Add sample ID
    data['sample'] = "${sample_id}"
    
    # Save data
    data.to_csv("${sample_id}_data.csv", index=False)
    
    # Validate and save validation results
    validation_results = loader.validate_data(data)
    
    # Convert numpy types for JSON serialization
    def convert_numpy(obj):
        import numpy as np
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj
    
    with open("${sample_id}_validation.json", 'w') as f:
        json.dump(validation_results, f, indent=2, default=convert_numpy)
    
    print(f"Loaded {len(data)} methylation sites for sample ${sample_id}")
    """
}

process EXTRACT_FEATURES {
    
    input:
    path(sample_files)
    
    output:
    path("features.csv"), emit: features
    path("feature_names.txt"), emit: feature_names
    path("feature_summary.json"), emit: summary
    
    script:
    window_sizes_list = params.window_sizes.split(',').collect { it as Integer }
    """
    #!/usr/bin/env python3
    
    import sys
    sys.path.append('${projectDir}')
    
    from feature_engineering import MethylationFeatureExtractor
    import pandas as pd
    import json
    
    # Read all sample data files
    sample_files = "${sample_files}".split()
    all_data = []
    
    for file_path in sample_files:
        if file_path.endswith('.csv'):
            data = pd.read_csv(file_path)
            all_data.append(data)
    
    # Combine all samples
    combined_data = pd.concat(all_data, ignore_index=True)
    
    print(f"Combined data from {len(all_data)} samples: {len(combined_data)} methylation sites")
    
    # Initialize feature extractor
    window_sizes = [${window_sizes_list.join(',')}]
    extractor = MethylationFeatureExtractor(window_sizes=window_sizes)
    
    # Extract all features
    features_df = extractor.combine_all_features(combined_data)
    
    # Scale features
    if "${params.scaling_method}" != "none":
        scaled_features = extractor.scale_features(
            features_df, 
            method="${params.scaling_method}", 
            fit_scaler=True
        )
    else:
        scaled_features = features_df
    
    # Save features
    scaled_features.to_csv("features.csv", index=False)
    
    # Save feature names
    feature_names = [col for col in scaled_features.columns if col != 'sample_id']
    with open("feature_names.txt", 'w') as f:
        for name in feature_names:
            f.write(f"{name}\\n")
    
    # Save summary
    summary = {
        'n_samples': len(scaled_features),
        'n_features': len(feature_names),
        'samples': scaled_features['sample_id'].tolist(),
        'scaling_method': "${params.scaling_method}",
        'window_sizes': window_sizes
    }
    
    with open("feature_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Generated {len(feature_names)} features for {len(scaled_features)} samples")
    """
}

process PREPARE_METADATA {
    
    input:
    path(metadata_file)
    path(features_file)
    
    output:
    path("prepared_metadata.csv"), emit: metadata
    path("metadata_summary.json"), emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    
    import sys
    sys.path.append('${projectDir}')
    
    import pandas as pd
    import json
    
    # Load metadata and features
    metadata = pd.read_csv("${metadata_file}")
    features = pd.read_csv("${features_file}")
    
    # Get sample IDs from features
    feature_samples = set(features['sample_id'])
    
    # Filter metadata to samples present in features
    valid_metadata = metadata[metadata['sample_id'].isin(feature_samples)]
    
    # Check for missing samples
    metadata_samples = set(valid_metadata['sample_id'])
    missing_in_metadata = feature_samples - metadata_samples
    missing_in_features = metadata_samples - feature_samples
    
    if missing_in_metadata:
        print(f"Warning: Samples in features but not in metadata: {missing_in_metadata}")
    
    if missing_in_features:
        print(f"Warning: Samples in metadata but not in features: {missing_in_features}")
    
    # Save prepared metadata
    valid_metadata.to_csv("prepared_metadata.csv", index=False)
    
    # Create summary
    summary = {
        'total_samples_in_metadata': len(metadata),
        'valid_samples': len(valid_metadata),
        'samples_with_features': len(feature_samples),
        'mgmt_status_distribution': valid_metadata['mgmt_status'].value_counts().to_dict() if 'mgmt_status' in valid_metadata.columns else {},
        'missing_in_metadata': list(missing_in_metadata),
        'missing_in_features': list(missing_in_features)
    }
    
    with open("metadata_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Prepared metadata for {len(valid_metadata)} samples")
    """
}

process CREATE_METADATA_TEMPLATE {
    
    input:
    path(features_file)
    
    output:
    path("metadata_template.csv"), emit: template
    
    script:
    """
    #!/usr/bin/env python3
    
    import sys
    sys.path.append('${projectDir}')
    
    from data_loader import create_sample_metadata_template
    import pandas as pd
    
    # Load features to get sample IDs
    features = pd.read_csv("${features_file}")
    sample_ids = features['sample_id'].tolist()
    
    # Create metadata template
    create_sample_metadata_template("metadata_template.csv", sample_ids)
    
    print(f"Created metadata template for {len(sample_ids)} samples")
    print("Please fill in the MGMT status labels in metadata_template.csv")
    """
}

process TRAIN_MODELS {
    
    input:
    path(features_file)
    path(metadata_file)
    path(feature_names_file)
    
    output:
    path("models/"), emit: models
    path("best_model.joblib"), emit: best_model
    path("training_results.json"), emit: training_results
    path("feature_selector.joblib"), emit: feature_selector, optional: true
    path("label_encoder.joblib"), emit: label_encoder, optional: true
    
    script:
    """
    #!/usr/bin/env python3
    
    import sys
    sys.path.append('${projectDir}')
    
    from model_training import train_mgmt_model_pipeline
    import pandas as pd
    import json
    import joblib
    from pathlib import Path
    
    # Load data
    features_df = pd.read_csv("${features_file}")
    metadata_df = pd.read_csv("${metadata_file}")
    
    # Check if we have valid labels
    if metadata_df['mgmt_status'].str.contains('unknown').any():
        print("Warning: Some samples have unknown MGMT status.")
        print("Creating dummy models for pipeline demonstration.")
        
        # Create dummy outputs
        Path("models").mkdir(exist_ok=True)
        
        # Create dummy results
        results = {
            'models_trained': False,
            'reason': 'No valid labels provided',
            'n_samples': len(features_df),
            'n_features': len([col for col in features_df.columns if col != 'sample_id'])
        }
        
        with open("training_results.json", 'w') as f:
            json.dump(results, f, indent=2)
        
        # Create dummy model file
        import pickle
        with open("best_model.joblib", 'wb') as f:
            pickle.dump(None, f)
            
    else:
        # Train models with valid labels
        trainer = train_mgmt_model_pipeline(
            features_df=features_df,
            metadata_df=metadata_df,
            output_dir="models",
            feature_selection=${params.feature_selection},
            n_features=${params.n_features},
            hyperparameter_tuning=${params.hyperparameter_tuning},
            create_ensemble=${params.create_ensemble}
        )
        
        # Save best model separately
        if trainer.best_model is not None:
            joblib.dump(trainer.best_model, "best_model.joblib")
        
        # Save feature selector and label encoder
        if trainer.feature_selector is not None:
            joblib.dump(trainer.feature_selector, "feature_selector.joblib")
        
        if hasattr(trainer, 'label_encoder'):
            joblib.dump(trainer.label_encoder, "label_encoder.joblib")
        
        # Create training results summary
        results = {
            'models_trained': True,
            'best_model_type': type(trainer.best_model).__name__ if trainer.best_model else None,
            'n_samples': len(features_df),
            'n_features': len([col for col in features_df.columns if col != 'sample_id']),
            'models_available': list(trainer.models.keys()) if trainer.models else [],
            'feature_selection_used': ${params.feature_selection},
            'n_selected_features': ${params.n_features} if ${params.feature_selection} else 'all'
        }
        
        with open("training_results.json", 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"Training completed. Best model: {results['best_model_type']}")
    """
}

process EVALUATE_MODELS {
    
    input:
    path(models_dir)
    path(features_file)
    path(metadata_file)
    path(feature_names_file)
    
    output:
    path("evaluation/"), emit: evaluation_results
    path("evaluation_report.txt"), emit: evaluation_report
    
    script:
    """
    #!/usr/bin/env python3
    
    import sys
    sys.path.append('${projectDir}')
    
    from model_evaluation import ModelEvaluator
    import pandas as pd
    import joblib
    from pathlib import Path
    import json
    
    # Load data
    features_df = pd.read_csv("${features_file}")
    metadata_df = pd.read_csv("${metadata_file}")
    
    # Check if we have valid models
    models_path = Path("${models_dir}")
    model_files = list(models_path.glob("*_model.joblib"))
    
    if not model_files:
        print("No trained models found. Creating dummy evaluation.")
        
        # Create dummy evaluation directory and report
        Path("evaluation").mkdir(exist_ok=True)
        
        with open("evaluation_report.txt", 'w') as f:
            f.write("MGMT Methylation Prediction - Evaluation Report\\n")
            f.write("=" * 60 + "\\n\\n")
            f.write("No trained models available for evaluation.\\n")
            f.write("This may be because no valid MGMT labels were provided.\\n")
            f.write("\\nPlease provide a metadata file with MGMT status labels to enable model training and evaluation.\\n")
    
    else:
        # Load models
        models = {}
        for model_file in model_files:
            model_name = model_file.stem.replace("_model", "")
            try:
                model = joblib.load(model_file)
                models[model_name] = model
            except Exception as e:
                print(f"Error loading model {model_name}: {e}")
        
        if models:
            # Initialize evaluator
            evaluator = ModelEvaluator(output_dir="evaluation")
            
            # Prepare test data
            data = features_df.merge(metadata_df, on='sample_id', how='inner')
            feature_columns = [col for col in features_df.columns if col != 'sample_id']
            
            X = data[feature_columns].values
            
            # Load label encoder if available
            label_encoder_file = models_path / "label_encoder.joblib"
            if label_encoder_file.exists():
                label_encoder = joblib.load(label_encoder_file)
                y = label_encoder.transform(data['mgmt_status'].values)
            else:
                # Create simple encoding
                from sklearn.preprocessing import LabelEncoder
                le = LabelEncoder()
                y = le.fit_transform(data['mgmt_status'].values)
            
            # Apply feature selection if used
            feature_selector_file = models_path / "feature_selector.joblib"
            if feature_selector_file.exists():
                feature_selector = joblib.load(feature_selector_file)
                X = feature_selector.transform(X)
                feature_columns = [feature_columns[i] for i in feature_selector.get_support(indices=True)]
            
            # Generate evaluation report
            report_path = evaluator.generate_evaluation_report(
                models=models,
                X_test=X,
                y_test=y,
                feature_names=feature_columns
            )
            
            # Copy main report
            import shutil
            shutil.copy(report_path, "evaluation_report.txt")
            
            print(f"Model evaluation completed. Report: {report_path}")
        
        else:
            print("No valid models could be loaded for evaluation.")
    """
}

process PREDICT_SAMPLES {
    
    input:
    path(best_model_file)
    path(features_file)
    path(feature_names_file)
    
    output:
    path("predictions.csv"), emit: predictions
    
    script:
    """
    #!/usr/bin/env python3
    
    import sys
    sys.path.append('${projectDir}')
    
    import pandas as pd
    import joblib
    import numpy as np
    from pathlib import Path
    
    # Load features
    features_df = pd.read_csv("${features_file}")
    
    # Try to load the best model
    try:
        best_model = joblib.load("${best_model_file}")
        
        if best_model is None:
            raise ValueError("No trained model available")
        
        # Prepare features
        feature_columns = [col for col in features_df.columns if col != 'sample_id']
        X = features_df[feature_columns].values
        
        # Load feature selector if available
        feature_selector_file = Path("${projectDir}") / "feature_selector.joblib"
        if feature_selector_file.exists():
            feature_selector = joblib.load(feature_selector_file)
            X = feature_selector.transform(X)
        
        # Load label encoder if available
        label_encoder_file = Path("${projectDir}") / "label_encoder.joblib"
        if label_encoder_file.exists():
            label_encoder = joblib.load(label_encoder_file)
        else:
            # Create dummy label encoder
            from sklearn.preprocessing import LabelEncoder
            label_encoder = LabelEncoder()
            label_encoder.fit(['methylated', 'unmethylated'])
        
        # Make predictions
        predictions = best_model.predict(X)
        probabilities = best_model.predict_proba(X)
        
        # Create results DataFrame
        results_df = pd.DataFrame({
            'sample_id': features_df['sample_id'].values,
            'predicted_mgmt_status': label_encoder.inverse_transform(predictions),
            'probability_methylated': probabilities[:, 1],
            'probability_unmethylated': probabilities[:, 0],
            'confidence': np.max(probabilities, axis=1)
        })
        
        print(f"Generated predictions for {len(results_df)} samples")
        
    except Exception as e:
        print(f"Error making predictions: {e}")
        print("Creating dummy predictions")
        
        # Create dummy predictions
        results_df = pd.DataFrame({
            'sample_id': features_df['sample_id'],
            'predicted_mgmt_status': ['unknown'] * len(features_df),
            'probability_methylated': [0.5] * len(features_df),
            'probability_unmethylated': [0.5] * len(features_df),
            'confidence': [0.5] * len(features_df)
        })
    
    # Save predictions
    results_df.to_csv("predictions.csv", index=False)
    """
}

process GENERATE_REPORT {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path(validation_files)
    path(features_file)
    path(training_results)
    path(evaluation_report)
    path(predictions_file)
    
    output:
    path("final_report.html"), emit: report
    path("pipeline_summary.json"), emit: summary
    path("*"), emit: all_results
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    import pandas as pd
    from pathlib import Path
    from datetime import datetime
    
    # Load data
    features_df = pd.read_csv("${features_file}")
    predictions_df = pd.read_csv("${predictions_file}")
    
    with open("${training_results}", 'r') as f:
        training_results = json.load(f)
    
    # Read evaluation report
    with open("${evaluation_report}", 'r') as f:
        evaluation_text = f.read()
    
    # Collect validation results
    validation_files = "${validation_files}".split()
    total_sites = 0
    sample_summaries = []
    
    for val_file in validation_files:
        if val_file.endswith('.json'):
            with open(val_file, 'r') as f:
                val_data = json.load(f)
                total_sites += val_data.get('total_records', 0)
                sample_summaries.append(val_data)
    
    # Create pipeline summary
    summary = {
        'pipeline_version': '1.0.0',
        'execution_time': datetime.now().isoformat(),
        'parameters': {
            'min_coverage': ${params.min_coverage},
            'modification_types': "${params.modification_types}",
            'feature_selection': ${params.feature_selection},
            'n_features': ${params.n_features},
            'scaling_method': "${params.scaling_method}",
            'hyperparameter_tuning': ${params.hyperparameter_tuning},
            'create_ensemble': ${params.create_ensemble}
        },
        'data_summary': {
            'n_samples': len(features_df),
            'total_methylation_sites': total_sites,
            'n_features_generated': len([col for col in features_df.columns if col != 'sample_id'])
        },
        'training_summary': training_results,
        'predictions': {
            'n_samples_predicted': len(predictions_df),
            'methylated_samples': (predictions_df['predicted_mgmt_status'] == 'methylated').sum(),
            'unmethylated_samples': (predictions_df['predicted_mgmt_status'] == 'unmethylated').sum(),
            'unknown_samples': (predictions_df['predicted_mgmt_status'] == 'unknown').sum()
        }
    }
    
    # Save summary
    with open("pipeline_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Generate HTML report
    html_content = f'''
    <!DOCTYPE html>
    <html>
    <head>
        <title>MGMT Methylation Detection - Pipeline Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 10px; }}
            .section {{ margin: 20px 0; padding: 15px; border-left: 4px solid #007cba; }}
            .summary-table {{ border-collapse: collapse; width: 100%; }}
            .summary-table th, .summary-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            .summary-table th {{ background-color: #f2f2f2; }}
            .success {{ color: green; }}
            .warning {{ color: orange; }}
            .error {{ color: red; }}
            pre {{ background-color: #f8f8f8; padding: 10px; border-radius: 5px; overflow-x: auto; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>MGMT Methylation Detection Pipeline Report</h1>
            <p><strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p><strong>Pipeline Version:</strong> 1.0.0</p>
        </div>
        
        <div class="section">
            <h2>Data Summary</h2>
            <table class="summary-table">
                <tr><th>Metric</th><th>Value</th></tr>
                <tr><td>Number of Samples</td><td>{len(features_df)}</td></tr>
                <tr><td>Total Methylation Sites</td><td>{total_sites:,}</td></tr>
                <tr><td>Features Generated</td><td>{len([col for col in features_df.columns if col != 'sample_id'])}</td></tr>
                <tr><td>Minimum Coverage</td><td>{params.min_coverage}</td></tr>
                <tr><td>Modification Types</td><td>{params.modification_types}</td></tr>
            </table>
        </div>
        
        <div class="section">
            <h2>Model Training Results</h2>
            <p><strong>Models Trained:</strong> <span class="{'success' if training_results.get('models_trained') else 'warning'}">{training_results.get('models_trained', False)}</span></p>
            {f'<p><strong>Best Model:</strong> {training_results.get("best_model_type", "N/A")}</p>' if training_results.get('models_trained') else ''}
            {f'<p><strong>Reason:</strong> <span class="warning">{training_results.get("reason", "")}</span></p>' if not training_results.get('models_trained') else ''}
        </div>
        
        <div class="section">
            <h2>Predictions Summary</h2>
            <table class="summary-table">
                <tr><th>MGMT Status</th><th>Number of Samples</th></tr>
                <tr><td>Methylated</td><td>{(predictions_df['predicted_mgmt_status'] == 'methylated').sum()}</td></tr>
                <tr><td>Unmethylated</td><td>{(predictions_df['predicted_mgmt_status'] == 'unmethylated').sum()}</td></tr>
                <tr><td>Unknown</td><td>{(predictions_df['predicted_mgmt_status'] == 'unknown').sum()}</td></tr>
            </table>
        </div>
        
        <div class="section">
            <h2>Model Evaluation</h2>
            <pre>{evaluation_text}</pre>
        </div>
        
        <div class="section">
            <h2>Pipeline Parameters</h2>
            <table class="summary-table">
                <tr><th>Parameter</th><th>Value</th></tr>
                <tr><td>Feature Selection</td><td>{params.feature_selection}</td></tr>
                <tr><td>Number of Features</td><td>{params.n_features}</td></tr>
                <tr><td>Scaling Method</td><td>{params.scaling_method}</td></tr>
                <tr><td>Hyperparameter Tuning</td><td>{params.hyperparameter_tuning}</td></tr>
                <tr><td>Create Ensemble</td><td>{params.create_ensemble}</td></tr>
            </table>
        </div>
        
        <div class="section">
            <h2>Output Files</h2>
            <ul>
                <li><strong>pipeline_summary.json</strong> - Complete pipeline execution summary</li>
                <li><strong>predictions.csv</strong> - MGMT methylation predictions for all samples</li>
                <li><strong>features.csv</strong> - Extracted methylation features</li>
                <li><strong>evaluation/</strong> - Model evaluation results and plots</li>
                <li><strong>models/</strong> - Trained machine learning models</li>
            </ul>
        </div>
    </body>
    </html>
    '''
    
    # Save HTML report
    with open("final_report.html", 'w') as f:
        f.write(html_content)
    
    print("Pipeline completed successfully!")
    print(f"Results summary: {len(features_df)} samples, {total_sites:,} methylation sites")
    print(f"Predictions: {(predictions_df['predicted_mgmt_status'] == 'methylated').sum()} methylated, {(predictions_df['predicted_mgmt_status'] == 'unmethylated').sum()} unmethylated")
    """
}
