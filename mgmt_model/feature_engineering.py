"""
Feature engineering module for MGMT methylation detection.
Creates genome-wide methylation features for machine learning.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from sklearn.preprocessing import StandardScaler, RobustScaler
from scipy import stats
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MethylationFeatureExtractor:
    """
    Extract meaningful features from genome-wide methylation data for MGMT prediction.
    """
    
    def __init__(self, window_sizes: List[int] = [1000, 5000, 10000, 50000]):
        """
        Initialize the feature extractor.
        
        Args:
            window_sizes (List[int]): Different window sizes for regional analysis
        """
        self.window_sizes = window_sizes
        self.scaler = None
        
    def extract_basic_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Extract basic methylation features per sample.
        
        Args:
            df (pd.DataFrame): Methylation data with 'sample' column
            
        Returns:
            pd.DataFrame: Basic features per sample
        """
        logger.info("Extracting basic methylation features")
        
        features_list = []
        
        for sample_id, sample_data in df.groupby('sample'):
            # Filter to methylation sites only
            meth_data = sample_data[sample_data['modification_type'] == 'm'].copy()
            
            if len(meth_data) == 0:
                continue
                
            sample_features = {'sample_id': sample_id}
            
            # Global methylation statistics
            sample_features.update({
                'global_mean_methylation': meth_data['percent_modified'].mean(),
                'global_median_methylation': meth_data['percent_modified'].median(),
                'global_std_methylation': meth_data['percent_modified'].std(),
                'global_methylation_range': meth_data['percent_modified'].max() - meth_data['percent_modified'].min(),
                'global_mean_coverage': meth_data['coverage'].mean(),
                'global_median_coverage': meth_data['coverage'].median(),
                'total_cpg_sites': len(meth_data),
                'high_methylation_sites': (meth_data['percent_modified'] >= 80).sum(),
                'low_methylation_sites': (meth_data['percent_modified'] <= 20).sum(),
                'moderate_methylation_sites': ((meth_data['percent_modified'] > 20) & 
                                              (meth_data['percent_modified'] < 80)).sum()
            })
            
            # Methylation percentiles
            percentiles = [10, 25, 75, 90, 95, 99]
            for p in percentiles:
                sample_features[f'methylation_p{p}'] = np.percentile(meth_data['percent_modified'], p)
            
            # Coverage statistics
            sample_features.update({
                'high_coverage_sites': (meth_data['coverage'] >= 10).sum(),
                'low_coverage_sites': (meth_data['coverage'] < 5).sum(),
                'coverage_p90': np.percentile(meth_data['coverage'], 90),
                'coverage_p10': np.percentile(meth_data['coverage'], 10)
            })
            
            features_list.append(sample_features)
        
        return pd.DataFrame(features_list)
    
    def extract_chromosomal_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Extract methylation features per chromosome per sample.
        
        Args:
            df (pd.DataFrame): Methylation data
            
        Returns:
            pd.DataFrame: Chromosomal features per sample
        """
        logger.info("Extracting chromosomal methylation features")
        
        features_list = []
        
        for sample_id, sample_data in df.groupby('sample'):
            meth_data = sample_data[sample_data['modification_type'] == 'm'].copy()
            
            if len(meth_data) == 0:
                continue
                
            sample_features = {'sample_id': sample_id}
            
            # Per-chromosome statistics
            for chromosome, chr_data in meth_data.groupby('chromosome'):
                if len(chr_data) < 5:  # Skip chromosomes with too few sites
                    continue
                    
                prefix = f'chr_{chromosome}_'
                
                sample_features.update({
                    f'{prefix}mean_methylation': chr_data['percent_modified'].mean(),
                    f'{prefix}median_methylation': chr_data['percent_modified'].median(),
                    f'{prefix}std_methylation': chr_data['percent_modified'].std(),
                    f'{prefix}cpg_count': len(chr_data),
                    f'{prefix}high_meth_ratio': (chr_data['percent_modified'] >= 80).mean(),
                    f'{prefix}low_meth_ratio': (chr_data['percent_modified'] <= 20).mean(),
                    f'{prefix}mean_coverage': chr_data['coverage'].mean(),
                    f'{prefix}coverage_std': chr_data['coverage'].std()
                })
                
                # Genomic density (sites per Mb)
                if len(chr_data) > 1:
                    genomic_span = chr_data['end'].max() - chr_data['start'].min()
                    if genomic_span > 0:
                        sample_features[f'{prefix}density_per_mb'] = len(chr_data) / (genomic_span / 1e6)
            
            features_list.append(sample_features)
        
        return pd.DataFrame(features_list)
    
    def extract_regional_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Extract features based on genomic regions and windows.
        
        Args:
            df (pd.DataFrame): Methylation data
            
        Returns:
            pd.DataFrame: Regional features per sample
        """
        logger.info("Extracting regional methylation features")
        
        features_list = []
        
        for sample_id, sample_data in df.groupby('sample'):
            meth_data = sample_data[sample_data['modification_type'] == 'm'].copy()
            
            if len(meth_data) == 0:
                continue
                
            sample_features = {'sample_id': sample_id}
            
            # Window-based features
            for window_size in self.window_sizes:
                window_features = self._extract_window_features(meth_data, window_size)
                for key, value in window_features.items():
                    sample_features[f'window_{window_size}_{key}'] = value
            
            # CpG island-like features (high-density regions)
            island_features = self._extract_cpg_island_features(meth_data)
            for key, value in island_features.items():
                sample_features[f'island_{key}'] = value
            
            features_list.append(sample_features)
        
        return pd.DataFrame(features_list)
    
    def _extract_window_features(self, meth_data: pd.DataFrame, window_size: int) -> Dict:
        """
        Extract features within sliding windows.
        
        Args:
            meth_data (pd.DataFrame): Methylation data for one sample
            window_size (int): Window size in base pairs
            
        Returns:
            Dict: Window-based features
        """
        window_features = {}
        window_methylations = []
        window_coverages = []
        window_site_counts = []
        
        for chromosome, chr_data in meth_data.groupby('chromosome'):
            if len(chr_data) < 2:
                continue
                
            chr_data = chr_data.sort_values('start')
            start_pos = chr_data['start'].min()
            end_pos = chr_data['end'].max()
            
            # Sliding windows
            for window_start in range(int(start_pos), int(end_pos), window_size // 2):
                window_end = window_start + window_size
                
                window_data = chr_data[
                    (chr_data['start'] >= window_start) &
                    (chr_data['end'] <= window_end)
                ]
                
                if len(window_data) >= 3:  # Minimum sites per window
                    window_methylations.append(window_data['percent_modified'].mean())
                    window_coverages.append(window_data['coverage'].mean())
                    window_site_counts.append(len(window_data))
        
        if window_methylations:
            window_features.update({
                'mean_window_methylation': np.mean(window_methylations),
                'std_window_methylation': np.std(window_methylations),
                'max_window_methylation': np.max(window_methylations),
                'min_window_methylation': np.min(window_methylations),
                'window_methylation_range': np.max(window_methylations) - np.min(window_methylations),
                'mean_window_coverage': np.mean(window_coverages),
                'mean_sites_per_window': np.mean(window_site_counts),
                'high_meth_windows': np.sum(np.array(window_methylations) >= 70),
                'low_meth_windows': np.sum(np.array(window_methylations) <= 30),
                'variable_windows': np.sum(np.array(window_methylations) > 30) & 
                                  np.sum(np.array(window_methylations) < 70)
            })
        
        return window_features
    
    def _extract_cpg_island_features(self, meth_data: pd.DataFrame) -> Dict:
        """
        Extract features resembling CpG island characteristics.
        
        Args:
            meth_data (pd.DataFrame): Methylation data for one sample
            
        Returns:
            Dict: CpG island-like features
        """
        island_features = {}
        
        # Find high-density regions (CpG island-like)
        density_threshold = 10  # sites per kb
        high_density_regions = []
        
        for chromosome, chr_data in meth_data.groupby('chromosome'):
            if len(chr_data) < 5:
                continue
                
            chr_data = chr_data.sort_values('start')
            
            # Calculate local density
            for i in range(len(chr_data) - 4):
                region_data = chr_data.iloc[i:i+5]
                region_span = region_data['end'].max() - region_data['start'].min()
                
                if region_span > 0:
                    density = len(region_data) / (region_span / 1000)  # sites per kb
                    
                    if density >= density_threshold:
                        high_density_regions.append({
                            'methylation': region_data['percent_modified'].mean(),
                            'coverage': region_data['coverage'].mean(),
                            'span': region_span,
                            'site_count': len(region_data)
                        })
        
        if high_density_regions:
            island_df = pd.DataFrame(high_density_regions)
            island_features.update({
                'n_high_density_regions': len(high_density_regions),
                'mean_island_methylation': island_df['methylation'].mean(),
                'std_island_methylation': island_df['methylation'].std(),
                'mean_island_coverage': island_df['coverage'].mean(),
                'mean_island_span': island_df['span'].mean(),
                'hypermethylated_islands': (island_df['methylation'] >= 80).sum(),
                'hypomethylated_islands': (island_df['methylation'] <= 20).sum()
            })
        else:
            # Default values if no high-density regions found
            island_features = {
                'n_high_density_regions': 0,
                'mean_island_methylation': 0,
                'std_island_methylation': 0,
                'mean_island_coverage': 0,
                'mean_island_span': 0,
                'hypermethylated_islands': 0,
                'hypomethylated_islands': 0
            }
        
        return island_features
    
    def extract_mgmt_specific_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Extract features specifically relevant to MGMT methylation detection.
        
        Args:
            df (pd.DataFrame): Methylation data
            
        Returns:
            pd.DataFrame: MGMT-specific features
        """
        logger.info("Extracting MGMT-specific features")
        
        features_list = []
        
        # MGMT gene region (chr10:129,466,683-129,467,448)
        mgmt_chr = 'chr10'
        mgmt_start = 129466683
        mgmt_end = 129467448
        
        # Extended MGMT promoter region
        mgmt_extended_start = mgmt_start - 2000
        mgmt_extended_end = mgmt_end + 2000
        
        for sample_id, sample_data in df.groupby('sample'):
            meth_data = sample_data[sample_data['modification_type'] == 'm'].copy()
            
            if len(meth_data) == 0:
                continue
                
            sample_features = {'sample_id': sample_id}
            
            # MGMT promoter region features
            mgmt_data = meth_data[
                (meth_data['chromosome'] == mgmt_chr) &
                (meth_data['start'] >= mgmt_start) &
                (meth_data['end'] <= mgmt_end)
            ]
            
            if len(mgmt_data) > 0:
                sample_features.update({
                    'mgmt_promoter_mean_meth': mgmt_data['percent_modified'].mean(),
                    'mgmt_promoter_median_meth': mgmt_data['percent_modified'].median(),
                    'mgmt_promoter_max_meth': mgmt_data['percent_modified'].max(),
                    'mgmt_promoter_min_meth': mgmt_data['percent_modified'].min(),
                    'mgmt_promoter_std_meth': mgmt_data['percent_modified'].std(),
                    'mgmt_promoter_coverage': mgmt_data['coverage'].mean(),
                    'mgmt_promoter_sites': len(mgmt_data),
                    'mgmt_promoter_high_meth': (mgmt_data['percent_modified'] >= 80).sum(),
                    'mgmt_promoter_hypermeth_ratio': (mgmt_data['percent_modified'] >= 80).mean()
                })
            else:
                # Default values if no MGMT data
                sample_features.update({
                    'mgmt_promoter_mean_meth': 0,
                    'mgmt_promoter_median_meth': 0,
                    'mgmt_promoter_max_meth': 0,
                    'mgmt_promoter_min_meth': 0,
                    'mgmt_promoter_std_meth': 0,
                    'mgmt_promoter_coverage': 0,
                    'mgmt_promoter_sites': 0,
                    'mgmt_promoter_high_meth': 0,
                    'mgmt_promoter_hypermeth_ratio': 0
                })
            
            # Extended MGMT region features
            mgmt_extended_data = meth_data[
                (meth_data['chromosome'] == mgmt_chr) &
                (meth_data['start'] >= mgmt_extended_start) &
                (meth_data['end'] <= mgmt_extended_end)
            ]
            
            if len(mgmt_extended_data) > 0:
                sample_features.update({
                    'mgmt_extended_mean_meth': mgmt_extended_data['percent_modified'].mean(),
                    'mgmt_extended_sites': len(mgmt_extended_data),
                    'mgmt_extended_hypermeth_ratio': (mgmt_extended_data['percent_modified'] >= 80).mean()
                })
            else:
                sample_features.update({
                    'mgmt_extended_mean_meth': 0,
                    'mgmt_extended_sites': 0,
                    'mgmt_extended_hypermeth_ratio': 0
                })
            
            # Chr10-specific features (MGMT chromosome)
            chr10_data = meth_data[meth_data['chromosome'] == mgmt_chr]
            if len(chr10_data) > 0:
                sample_features.update({
                    'chr10_mean_meth': chr10_data['percent_modified'].mean(),
                    'chr10_std_meth': chr10_data['percent_modified'].std(),
                    'chr10_sites': len(chr10_data),
                    'chr10_hypermeth_ratio': (chr10_data['percent_modified'] >= 80).mean(),
                    'chr10_vs_global_meth_ratio': chr10_data['percent_modified'].mean() / 
                                                 meth_data['percent_modified'].mean() if len(meth_data) > 0 else 0
                })
            else:
                sample_features.update({
                    'chr10_mean_meth': 0,
                    'chr10_std_meth': 0,
                    'chr10_sites': 0,
                    'chr10_hypermeth_ratio': 0,
                    'chr10_vs_global_meth_ratio': 0
                })
            
            features_list.append(sample_features)
        
        return pd.DataFrame(features_list)
    
    def extract_statistical_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Extract advanced statistical features from methylation patterns.
        
        Args:
            df (pd.DataFrame): Methylation data
            
        Returns:
            pd.DataFrame: Statistical features
        """
        logger.info("Extracting statistical methylation features")
        
        features_list = []
        
        for sample_id, sample_data in df.groupby('sample'):
            meth_data = sample_data[sample_data['modification_type'] == 'm'].copy()
            
            if len(meth_data) == 0:
                continue
                
            sample_features = {'sample_id': sample_id}
            
            methylation_values = meth_data['percent_modified'].values
            coverage_values = meth_data['coverage'].values
            
            # Distribution features
            sample_features.update({
                'methylation_skewness': stats.skew(methylation_values),
                'methylation_kurtosis': stats.kurtosis(methylation_values),
                'methylation_entropy': stats.entropy(np.histogram(methylation_values, bins=10)[0] + 1),
                'coverage_skewness': stats.skew(coverage_values),
                'coverage_kurtosis': stats.kurtosis(coverage_values)
            })
            
            # Quantile-based features
            q25, q50, q75 = np.percentile(methylation_values, [25, 50, 75])
            sample_features.update({
                'iqr_methylation': q75 - q25,
                'q75_q25_ratio': q75 / q25 if q25 > 0 else 0,
                'median_mad': np.median(np.abs(methylation_values - q50))
            })
            
            # Bimodality features (characteristic of methylation data)
            hist, bins = np.histogram(methylation_values, bins=20, range=(0, 100))
            sample_features.update({
                'bimodality_coefficient': self._calculate_bimodality_coefficient(methylation_values),
                'mode_count': len([i for i in range(1, len(hist)-1) if hist[i] > hist[i-1] and hist[i] > hist[i+1]]),
                'valley_depth': self._calculate_valley_depth(hist)
            })
            
            # Spatial correlation features (if data is sorted)
            if len(meth_data) > 10:
                sorted_data = meth_data.sort_values(['chromosome', 'start'])
                meth_values = sorted_data['percent_modified'].values
                
                # Auto-correlation at different lags
                for lag in [1, 5, 10]:
                    if len(meth_values) > lag:
                        correlation = np.corrcoef(meth_values[:-lag], meth_values[lag:])[0, 1]
                        sample_features[f'autocorr_lag_{lag}'] = correlation if not np.isnan(correlation) else 0
            
            features_list.append(sample_features)
        
        return pd.DataFrame(features_list)
    
    def _calculate_bimodality_coefficient(self, data: np.ndarray) -> float:
        """Calculate bimodality coefficient."""
        try:
            skewness = stats.skew(data)
            kurtosis = stats.kurtosis(data)
            n = len(data)
            
            # Bimodality coefficient formula
            numerator = skewness**2 + 1
            denominator = kurtosis + 3 * (n-1)**2 / ((n-2)*(n-3))
            
            return numerator / denominator if denominator != 0 else 0
        except:
            return 0
    
    def _calculate_valley_depth(self, histogram: np.ndarray) -> float:
        """Calculate the depth of valleys in histogram."""
        try:
            if len(histogram) < 3:
                return 0
            
            # Find local minima
            minima = []
            for i in range(1, len(histogram)-1):
                if histogram[i] < histogram[i-1] and histogram[i] < histogram[i+1]:
                    minima.append(histogram[i])
            
            # Find local maxima
            maxima = []
            for i in range(1, len(histogram)-1):
                if histogram[i] > histogram[i-1] and histogram[i] > histogram[i+1]:
                    maxima.append(histogram[i])
            
            if minima and maxima:
                return np.mean(maxima) - np.mean(minima)
            else:
                return 0
        except:
            return 0
    
    def combine_all_features(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Extract all feature types and combine them.
        
        Args:
            df (pd.DataFrame): Methylation data
            
        Returns:
            pd.DataFrame: Combined feature matrix
        """
        logger.info("Extracting all methylation features")
        
        # Extract different feature types
        basic_features = self.extract_basic_features(df)
        regional_features = self.extract_regional_features(df)
        mgmt_features = self.extract_mgmt_specific_features(df)
        stat_features = self.extract_statistical_features(df)
        
        # Merge all features
        all_features = basic_features
        
        for feature_df in [regional_features, mgmt_features, stat_features]:
            all_features = all_features.merge(feature_df, on='sample_id', how='outer')
        
        # Fill missing values
        all_features = all_features.fillna(0)
        
        logger.info(f"Generated {all_features.shape[1]-1} features for {len(all_features)} samples")
        
        return all_features
    
    def scale_features(self, features_df: pd.DataFrame, 
                      method: str = 'robust', 
                      fit_scaler: bool = True) -> pd.DataFrame:
        """
        Scale features for machine learning.
        
        Args:
            features_df (pd.DataFrame): Feature matrix
            method (str): Scaling method ('standard', 'robust')
            fit_scaler (bool): Whether to fit scaler or use existing
            
        Returns:
            pd.DataFrame: Scaled features
        """
        feature_columns = [col for col in features_df.columns if col != 'sample_id']
        
        if fit_scaler:
            if method == 'standard':
                self.scaler = StandardScaler()
            elif method == 'robust':
                self.scaler = RobustScaler()
            else:
                raise ValueError("Method must be 'standard' or 'robust'")
            
            scaled_features = self.scaler.fit_transform(features_df[feature_columns])
        else:
            if self.scaler is None:
                raise ValueError("Scaler must be fitted first")
            scaled_features = self.scaler.transform(features_df[feature_columns])
        
        # Create scaled DataFrame
        scaled_df = pd.DataFrame(scaled_features, columns=feature_columns)
        scaled_df['sample_id'] = features_df['sample_id'].values
        
        return scaled_df

if __name__ == "__main__":
    # Example usage
    from data_loader import BedMethylLoader
    
    # Load data
    loader = BedMethylLoader()
    sample_file = "/home/chbope/extension/nWGS_manuscript_data/data/results/epi2me/modkit/T001.wf_mods.bedmethyl.gz"
    
    df = loader.load_bedmethyl_file(sample_file, modification_types=['m'], min_coverage=3)
    df['sample'] = 'T001'
    
    # Extract features
    extractor = MethylationFeatureExtractor()
    
    # Get all features
    all_features = extractor.combine_all_features(df)
    print(f"Generated feature matrix: {all_features.shape}")
    
    # Display some key features
    key_features = ['sample_id', 'global_mean_methylation', 'mgmt_promoter_mean_meth', 
                   'chr10_mean_meth', 'methylation_skewness']
    print("\nKey features:")
    print(all_features[key_features])

