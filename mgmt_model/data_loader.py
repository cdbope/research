"""
Data loading module for MGMT methylation detection.
Handles bedmethyl files from modkit package.
"""

import pandas as pd
import numpy as np
import gzip
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class BedMethylLoader:
    """
    Class to load and preprocess bedmethyl files for MGMT methylation analysis.
    """
    
    def __init__(self):
        self.bedmethyl_columns = [
            'chromosome', 'start', 'end', 'modification_type', 'score', 'strand',
            'thick_start', 'thick_end', 'rgb', 'coverage', 'percent_modified',
            'n_modified', 'n_canonical', 'n_other_mod', 'n_delete', 'n_fail', 'n_diff'
        ]
        
    def load_bedmethyl_file(self, filepath: str, 
                           chromosomes: Optional[List[str]] = None,
                           modification_types: Optional[List[str]] = None,
                           min_coverage: int = 1) -> pd.DataFrame:
        """
        Load bedmethyl file with optional filtering.
        
        Args:
            filepath (str): Path to the bedmethyl file (.gz supported)
            chromosomes (List[str], optional): List of chromosomes to include
            modification_types (List[str], optional): Modification types to include ('m', 'h')
            min_coverage (int): Minimum coverage threshold
            
        Returns:
            pd.DataFrame: Loaded and filtered methylation data
        """
        logger.info(f"Loading bedmethyl file: {filepath}")
        
        try:
            # Determine if file is gzipped
            if filepath.endswith('.gz'):
                open_func = gzip.open
                mode = 'rt'
            else:
                open_func = open
                mode = 'r'
            
            # Read the data
            data = []
            with open_func(filepath, mode) as f:
                for line in f:
                    if line.startswith('#'):  # Skip header lines
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 11:  # Ensure minimum required columns
                        data.append(fields[:17])  # Take first 17 columns
            
            # Create DataFrame
            df = pd.DataFrame(data, columns=self.bedmethyl_columns[:len(data[0]) if data else 0])
            
            # Convert numeric columns
            numeric_columns = ['start', 'end', 'score', 'thick_start', 'thick_end', 
                              'coverage', 'percent_modified', 'n_modified', 'n_canonical']
            for col in numeric_columns:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            logger.info(f"Loaded {len(df)} records")
            
            # Apply filters
            if chromosomes:
                df = df[df['chromosome'].isin(chromosomes)]
                logger.info(f"Filtered to chromosomes {chromosomes}: {len(df)} records")
            
            if modification_types:
                df = df[df['modification_type'].isin(modification_types)]
                logger.info(f"Filtered to modification types {modification_types}: {len(df)} records")
            
            if min_coverage > 1:
                df = df[df['coverage'] >= min_coverage]
                logger.info(f"Filtered to min coverage {min_coverage}: {len(df)} records")
            
            return df
            
        except Exception as e:
            logger.error(f"Error loading bedmethyl file: {e}")
            raise
    
    def load_multiple_samples(self, file_paths: List[str], 
                             sample_names: Optional[List[str]] = None,
                             **kwargs) -> Dict[str, pd.DataFrame]:
        """
        Load multiple bedmethyl files.
        
        Args:
            file_paths (List[str]): List of file paths
            sample_names (List[str], optional): Sample names for each file
            **kwargs: Additional arguments for load_bedmethyl_file
            
        Returns:
            Dict[str, pd.DataFrame]: Dictionary mapping sample names to data
        """
        if sample_names is None:
            sample_names = [Path(fp).stem.replace('.wf_mods.bedmethyl', '') for fp in file_paths]
        
        if len(sample_names) != len(file_paths):
            raise ValueError("Number of sample names must match number of file paths")
        
        samples_data = {}
        for sample_name, file_path in zip(sample_names, file_paths):
            logger.info(f"Loading sample: {sample_name}")
            df = self.load_bedmethyl_file(file_path, **kwargs)
            df['sample'] = sample_name
            samples_data[sample_name] = df
        
        return samples_data
    
    def get_chromosome_summary(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Get summary statistics per chromosome.
        
        Args:
            df (pd.DataFrame): Methylation data
            
        Returns:
            pd.DataFrame: Summary statistics per chromosome
        """
        summary = df.groupby(['chromosome', 'modification_type']).agg({
            'percent_modified': ['count', 'mean', 'std', 'median'],
            'coverage': ['mean', 'std'],
            'start': ['min', 'max']
        }).round(2)
        
        # Flatten column names
        summary.columns = ['_'.join(col).strip() for col in summary.columns.values]
        summary = summary.reset_index()
        
        return summary
    
    def get_positional_data(self, df: pd.DataFrame, 
                           chromosome: str, 
                           start_pos: int, 
                           end_pos: int) -> pd.DataFrame:
        """
        Extract methylation data for a specific genomic region.
        
        Args:
            df (pd.DataFrame): Methylation data
            chromosome (str): Chromosome name
            start_pos (int): Start position
            end_pos (int): End position
            
        Returns:
            pd.DataFrame: Filtered data for the region
        """
        region_data = df[
            (df['chromosome'] == chromosome) &
            (df['start'] >= start_pos) &
            (df['end'] <= end_pos)
        ].copy()
        
        return region_data.sort_values('start')
    
    def validate_data(self, df: pd.DataFrame) -> Dict[str, any]:
        """
        Validate the loaded data and return quality metrics.
        
        Args:
            df (pd.DataFrame): Methylation data
            
        Returns:
            Dict: Validation results and quality metrics
        """
        validation_results = {
            'total_records': len(df),
            'chromosomes': sorted(df['chromosome'].unique().tolist()),
            'modification_types': sorted(df['modification_type'].unique().tolist()),
            'coverage_stats': {
                'mean': df['coverage'].mean(),
                'median': df['coverage'].median(),
                'min': df['coverage'].min(),
                'max': df['coverage'].max()
            },
            'methylation_stats': {
                'mean': df['percent_modified'].mean(),
                'median': df['percent_modified'].median(),
                'min': df['percent_modified'].min(),
                'max': df['percent_modified'].max()
            },
            'missing_values': df.isnull().sum().to_dict(),
            'genomic_span': {
                chr_name: {
                    'start': chr_data['start'].min(),
                    'end': chr_data['end'].max(),
                    'span': chr_data['end'].max() - chr_data['start'].min()
                }
                for chr_name, chr_data in df.groupby('chromosome')
            }
        }
        
        return validation_results

def load_sample_metadata(metadata_file: str) -> pd.DataFrame:
    """
    Load sample metadata file containing MGMT status labels.
    
    Args:
        metadata_file (str): Path to metadata CSV file
        
    Returns:
        pd.DataFrame: Sample metadata with MGMT labels
    """
    try:
        metadata = pd.read_csv(metadata_file)
        required_columns = ['sample_id', 'mgmt_status']
        
        if not all(col in metadata.columns for col in required_columns):
            raise ValueError(f"Metadata file must contain columns: {required_columns}")
        
        logger.info(f"Loaded metadata for {len(metadata)} samples")
        return metadata
    
    except Exception as e:
        logger.error(f"Error loading metadata: {e}")
        raise

def create_sample_metadata_template(output_file: str, sample_names: List[str]) -> None:
    """
    Create a template metadata file for sample labels.
    
    Args:
        output_file (str): Output CSV file path
        sample_names (List[str]): List of sample names
    """
    template_df = pd.DataFrame({
        'sample_id': sample_names,
        'mgmt_status': ['unknown'] * len(sample_names),  # To be filled by user
        'batch': ['batch1'] * len(sample_names),
        'tissue_type': ['tumor'] * len(sample_names),
        'notes': [''] * len(sample_names)
    })
    
    template_df.to_csv(output_file, index=False)
    logger.info(f"Created metadata template: {output_file}")

if __name__ == "__main__":
    # Example usage
    loader = BedMethylLoader()
    
    # Load single sample
    sample_file = "/home/chbope/extension/nWGS_manuscript_data/data/results/epi2me/modkit/T001.wf_mods.bedmethyl.gz"
    
    # Load with basic filtering
    df = loader.load_bedmethyl_file(
        sample_file,
        modification_types=['m'],  # Only methylation
        min_coverage=3
    )
    
    print(f"Loaded {len(df)} methylation sites")
    
    # Get summary
    summary = loader.get_chromosome_summary(df)
    print("\nChromosome summary:")
    print(summary.head())
    
    # Validate data
    validation = loader.validate_data(df)
    print(f"\nValidation results:")
    print(f"Chromosomes: {validation['chromosomes'][:5]}...")
    print(f"Mean coverage: {validation['coverage_stats']['mean']:.2f}")
    print(f"Mean methylation: {validation['methylation_stats']['mean']:.2f}%")

