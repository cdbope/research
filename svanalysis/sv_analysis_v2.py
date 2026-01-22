#!/usr/bin/env python3
"""
Structural Variant Analysis for GBM Samples (Version 2)
Analyzes Sniffles VCF files for SV counts, VAF, and clustering

VERSION 2 IMPROVEMENTS:
- Focus on SIGNIFICANT VAF only (high-quality, high-support variants)
- Clonal vs subclonal VAF analysis
- Weighted VAF metrics
- Top percentile VAF (top 10%, 25%)
- Chromosome-specific VAF tracking
- Quality filtering (PASS, PRECISE, min support)
"""

import os
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE PATHS FOR YOUR SYSTEM
# ============================================================================

# Path to sample.txt file (contains sample_id and purity values)
# Format: sample_id\tpurity (tab-separated)
SAMPLE_FILE = "/home/chbope/extension/script/svanalysis/sample.txt"

# Directory containing VCF files (*.vcf or *.vcf.gz)
# VCF files should be named: {sample_id}.vcf or {sample_id}.vcf.gz
VCF_DIR = "/home/chbope/extension/script/svanalysis"

# Output directory for all results (tables, plots, reports)
# Will be created if it doesn't exist
OUTPUT_DIR = "/home/chbope/extension/script/svanalysis/results_v2"

# ============================================================================
# QUALITY FILTERING PARAMETERS (OPTIONAL - Advanced users only)
# ============================================================================

# Minimum read support for a variant to be considered high-quality
MIN_SUPPORT = 10

# Minimum VAF for a variant to be included in analysis
MIN_VAF = 0.1

# Clonality threshold: variants with VAF >= this are considered "clonal"
CLONAL_THRESHOLD = 0.3

# ============================================================================

class SVAnalyzer:
    def __init__(self, sample_file, vcf_dir, output_dir):
        """
        Initialize SV Analyzer

        Args:
            sample_file: Path to sample.txt with sample_id and purity
            vcf_dir: Directory containing VCF files
            output_dir: Directory for output files
        """
        self.sample_file = sample_file
        self.vcf_dir = vcf_dir
        self.output_dir = output_dir

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Load sample information
        self.samples = self.load_samples()

    def load_samples(self):
        """Load sample information from sample.txt"""
        samples = {}
        with open(self.sample_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    if len(parts) == 2:
                        sample_id = parts[0]
                        purity = float(parts[1])
                        samples[sample_id] = {'purity': purity}
        return samples

    def parse_vcf(self, vcf_path, min_support=10, min_vaf=0.1):
        """
        Parse a VCF file and extract SV information with quality filtering

        Args:
            vcf_path: Path to VCF file (gzipped)
            min_support: Minimum read support for high-quality SVs (default: 10)
            min_vaf: Minimum VAF to consider variant (default: 0.1)

        Returns:
            Dictionary with SV counts and VAF information
        """
        sv_counts = {'DEL': 0, 'DUP': 0, 'INV': 0, 'BND': 0, 'INS': 0}
        sv_counts_hq = {'DEL': 0, 'DUP': 0, 'INV': 0, 'BND': 0, 'INS': 0}  # High quality

        all_vaf_values = []
        high_quality_vaf = []  # VAF from high-support variants
        clonal_vaf = []        # VAF >= 0.3 (likely clonal)
        subclonal_vaf = []     # 0.1 <= VAF < 0.3 (likely subclonal)

        # For chromosome-specific analysis
        chrom_vaf = {}

        open_func = gzip.open if vcf_path.endswith('.gz') else open

        with open_func(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue

                chrom = fields[0]
                info = fields[7]
                filter_field = fields[6]

                # Extract SV type, AF, and SUPPORT
                svtype = None
                af = None
                support = 0
                is_precise = 'IMPRECISE' not in info

                for item in info.split(';'):
                    if item.startswith('SVTYPE='):
                        svtype = item.split('=')[1]
                    elif item.startswith('AF='):
                        try:
                            af = float(item.split('=')[1])
                        except:
                            pass
                    elif item.startswith('SUPPORT='):
                        try:
                            support = int(item.split('=')[1])
                        except:
                            pass

                # Count SV types (all variants)
                if svtype in sv_counts:
                    sv_counts[svtype] += 1

                # Skip low-quality variants for detailed analysis
                if af is None or af < min_vaf:
                    continue

                # Collect all VAF values
                all_vaf_values.append(af)

                # High-quality variants: good support + PASS filter
                is_high_quality = (support >= min_support and
                                 filter_field == 'PASS' and
                                 is_precise)

                if is_high_quality:
                    if svtype in sv_counts_hq:
                        sv_counts_hq[svtype] += 1
                    high_quality_vaf.append(af)

                # Categorize by clonality (using threshold from configuration)
                if af >= CLONAL_THRESHOLD:
                    clonal_vaf.append(af)
                elif af >= min_vaf:
                    subclonal_vaf.append(af)

                # Track chromosome-specific VAF
                if chrom not in chrom_vaf:
                    chrom_vaf[chrom] = []
                chrom_vaf[chrom].append(af)

        # Calculate comprehensive VAF metrics
        vaf_metrics = {
            # All variants
            'mean_vaf_all': np.mean(all_vaf_values) if all_vaf_values else 0,
            'median_vaf_all': np.median(all_vaf_values) if all_vaf_values else 0,

            # High-quality variants (RECOMMENDED for clustering)
            'mean_vaf_hq': np.mean(high_quality_vaf) if high_quality_vaf else 0,
            'median_vaf_hq': np.median(high_quality_vaf) if high_quality_vaf else 0,
            'hq_variant_count': len(high_quality_vaf),

            # Clonal variants (VAF >= 0.3, likely driver events)
            'mean_vaf_clonal': np.mean(clonal_vaf) if clonal_vaf else 0,
            'median_vaf_clonal': np.median(clonal_vaf) if clonal_vaf else 0,
            'clonal_fraction': len(clonal_vaf) / len(all_vaf_values) if all_vaf_values else 0,
            'clonal_count': len(clonal_vaf),

            # Subclonal variants (0.1 <= VAF < 0.3)
            'mean_vaf_subclonal': np.mean(subclonal_vaf) if subclonal_vaf else 0,
            'subclonal_count': len(subclonal_vaf),

            # Weighted VAF (emphasizes high-VAF variants)
            'weighted_vaf': self._calculate_weighted_vaf(high_quality_vaf) if high_quality_vaf else 0,

            # Top 10% VAF (focus on most significant variants)
            'top10_vaf': self._calculate_top_percentile_vaf(high_quality_vaf, 10) if high_quality_vaf else 0,

            # Top 25% VAF
            'top25_vaf': self._calculate_top_percentile_vaf(high_quality_vaf, 25) if high_quality_vaf else 0,
        }

        return {
            'sv_counts': sv_counts,
            'sv_counts_hq': sv_counts_hq,
            'vaf_values': all_vaf_values,
            'high_quality_vaf': high_quality_vaf,
            'clonal_vaf': clonal_vaf,
            'subclonal_vaf': subclonal_vaf,
            'chrom_vaf': chrom_vaf,
            **vaf_metrics
        }

    def _calculate_weighted_vaf(self, vaf_list):
        """
        Calculate weighted mean VAF, giving more weight to higher VAF variants
        (more likely to be significant/driver events)
        """
        if not vaf_list:
            return 0
        vaf_array = np.array(vaf_list)
        # Use VAF itself as weight (higher VAF = more weight)
        weights = vaf_array
        return np.average(vaf_array, weights=weights)

    def _calculate_top_percentile_vaf(self, vaf_list, percentile):
        """
        Calculate mean VAF of top N percentile variants
        Focuses on the most significant variants
        """
        if not vaf_list:
            return 0
        threshold = np.percentile(vaf_list, 100 - percentile)
        top_vaf = [v for v in vaf_list if v >= threshold]
        return np.mean(top_vaf) if top_vaf else 0

    def analyze_all_samples(self):
        """Analyze all VCF files and create summary table with enhanced VAF metrics"""
        results = []

        for sample_id, sample_info in self.samples.items():
            vcf_file = os.path.join(self.vcf_dir, f"{sample_id}.wf_sv.vcf.gz")

            if not os.path.exists(vcf_file):
                print(f"Warning: VCF file not found for {sample_id}")
                continue

            print(f"Processing {sample_id}...")
            vcf_data = self.parse_vcf(vcf_file, min_support=MIN_SUPPORT, min_vaf=MIN_VAF)

            result = {
                'Sample': sample_id,
                'Purity': sample_info['purity'],

                # SV counts (all variants)
                'DEL': vcf_data['sv_counts']['DEL'],
                'DUP': vcf_data['sv_counts']['DUP'],
                'INV': vcf_data['sv_counts']['INV'],
                'BND': vcf_data['sv_counts']['BND'],
                'INS': vcf_data['sv_counts']['INS'],

                # High-quality SV counts
                'DEL_HQ': vcf_data['sv_counts_hq']['DEL'],
                'DUP_HQ': vcf_data['sv_counts_hq']['DUP'],
                'INV_HQ': vcf_data['sv_counts_hq']['INV'],
                'BND_HQ': vcf_data['sv_counts_hq']['BND'],
                'INS_HQ': vcf_data['sv_counts_hq']['INS'],

                # VAF metrics - all variants (for reference)
                'Mean_VAF_All': vcf_data['mean_vaf_all'],
                'Median_VAF_All': vcf_data['median_vaf_all'],

                # VAF metrics - HIGH QUALITY (RECOMMENDED)
                'Mean_VAF_HQ': vcf_data['mean_vaf_hq'],
                'Median_VAF_HQ': vcf_data['median_vaf_hq'],
                'HQ_Variant_Count': vcf_data['hq_variant_count'],

                # VAF metrics - CLONAL (most significant)
                'Mean_VAF_Clonal': vcf_data['mean_vaf_clonal'],
                'Median_VAF_Clonal': vcf_data['median_vaf_clonal'],
                'Clonal_Fraction': vcf_data['clonal_fraction'],
                'Clonal_Count': vcf_data['clonal_count'],

                # VAF metrics - SUBCLONAL
                'Mean_VAF_Subclonal': vcf_data['mean_vaf_subclonal'],
                'Subclonal_Count': vcf_data['subclonal_count'],

                # Advanced VAF metrics
                'Weighted_VAF': vcf_data['weighted_vaf'],
                'Top10_VAF': vcf_data['top10_vaf'],
                'Top25_VAF': vcf_data['top25_vaf'],
            }
            results.append(result)

        self.df = pd.DataFrame(results)
        return self.df

    def save_summary_table(self):
        """Save summary table to CSV and display-friendly formats"""
        # Save full table with VAF
        self.df.to_csv(os.path.join(self.output_dir, 'sv_summary_table.csv'), index=False)

        # Create display table (like the figure)
        display_df = self.df[['Sample', 'Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS']].copy()
        display_df['Purity'] = display_df['Purity'].round(2)

        # Save display table
        display_df.to_csv(os.path.join(self.output_dir, 'sv_summary_display.csv'), index=False)

        print("\n" + "="*80)
        print("STRUCTURAL VARIANT SUMMARY TABLE")
        print("="*80)
        print(display_df.to_string(index=False))
        print("="*80 + "\n")

        return display_df

    def perform_clustering(self, method='ward', distance_metric='euclidean', n_clusters=None,
                          use_strategy='high_quality'):
        """
        Perform hierarchical clustering optimized for GBM samples

        For GBM samples, we consider:
        - SV burden (total counts)
        - SV type distribution
        - Purity
        - VAF metrics (clonality indicators) - USING SIGNIFICANT VAF ONLY

        Args:
            method: Linkage method ('ward', 'average', 'complete')
            distance_metric: Distance metric for clustering
            n_clusters: Number of clusters (if None, will be determined by dendrogram)
            use_strategy: VAF strategy for clustering
                - 'high_quality': Use high-quality VAF (RECOMMENDED)
                - 'clonal': Use clonal VAF only (VAF >= 0.3)
                - 'weighted': Use weighted VAF (emphasizes high VAF)
                - 'top10': Use top 10% VAF
                - 'all': Use all VAF (not recommended)
        """
        # Select features based on strategy
        if use_strategy == 'high_quality':
            vaf_cols = ['Mean_VAF_HQ', 'Median_VAF_HQ', 'Clonal_Fraction']
            title_suffix = "Using High-Quality VAF"
        elif use_strategy == 'clonal':
            vaf_cols = ['Mean_VAF_Clonal', 'Median_VAF_Clonal', 'Clonal_Fraction']
            title_suffix = "Using Clonal VAF (≥0.3)"
        elif use_strategy == 'weighted':
            vaf_cols = ['Weighted_VAF', 'Clonal_Fraction']
            title_suffix = "Using Weighted VAF"
        elif use_strategy == 'top10':
            vaf_cols = ['Top10_VAF', 'Top25_VAF', 'Clonal_Fraction']
            title_suffix = "Using Top Percentile VAF"
        else:  # 'all'
            vaf_cols = ['Mean_VAF_All', 'Median_VAF_All']
            title_suffix = "Using All VAF"

        # Prepare features for clustering
        feature_cols = ['Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS'] + vaf_cols
        X = self.df[feature_cols].values
        sample_names = self.df['Sample'].values

        print(f"\nClustering with strategy: {use_strategy}")
        print(f"Features used: {', '.join(feature_cols)}")

        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Perform hierarchical clustering
        linkage_matrix = linkage(X_scaled, method=method, metric=distance_metric)

        # Create dendrogram
        plt.figure(figsize=(12, 6))
        dendrogram(linkage_matrix, labels=sample_names, leaf_rotation=45, leaf_font_size=10)
        plt.title(f'Hierarchical Clustering of GBM Samples\n{title_suffix} (Method: {method})')
        plt.xlabel('Sample')
        plt.ylabel('Distance')
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, f'clustering_dendrogram_{use_strategy}.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()

        # Determine optimal number of clusters if not specified
        if n_clusters is None:
            # Use elbow method or default to 3 clusters for GBM
            n_clusters = min(3, len(self.df) // 2)

        # Assign cluster labels
        clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
        self.df['Cluster'] = clusters

        # Save clustering results
        self.df.to_csv(os.path.join(self.output_dir, 'clustering_results.csv'), index=False)

        print(f"\nClustering Analysis (n_clusters={n_clusters}):")
        print("="*80)
        for cluster_id in range(1, n_clusters + 1):
            cluster_samples = self.df[self.df['Cluster'] == cluster_id]['Sample'].tolist()
            print(f"Cluster {cluster_id}: {', '.join(cluster_samples)}")
        print("="*80 + "\n")

        return linkage_matrix, clusters

    def perform_pca_clustering(self, use_hq_vaf=True):
        """
        Perform PCA for dimensionality reduction and visualization
        Useful for understanding sample relationships in GBM cohort

        Args:
            use_hq_vaf: If True, use high-quality VAF metrics (recommended)
        """
        if use_hq_vaf:
            feature_cols = ['Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS',
                          'Mean_VAF_HQ', 'Median_VAF_HQ', 'Clonal_Fraction']
            title_vaf = "High-Quality VAF"
        else:
            feature_cols = ['Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS',
                          'Mean_VAF_All', 'Median_VAF_All']
            title_vaf = "All VAF"

        X = self.df[feature_cols].values
        sample_names = self.df['Sample'].values

        # Standardize
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # PCA
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X_scaled)

        # Plot
        plt.figure(figsize=(10, 8))

        if 'Cluster' in self.df.columns:
            # Color by cluster if available
            scatter = plt.scatter(X_pca[:, 0], X_pca[:, 1],
                                c=self.df['Cluster'], cmap='Set1',
                                s=100, alpha=0.7, edgecolors='black')
            plt.colorbar(scatter, label='Cluster')
        else:
            plt.scatter(X_pca[:, 0], X_pca[:, 1], s=100, alpha=0.7, edgecolors='black')

        # Annotate samples
        for i, sample in enumerate(sample_names):
            plt.annotate(sample, (X_pca[i, 0], X_pca[i, 1]),
                        xytext=(5, 5), textcoords='offset points', fontsize=9)

        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        plt.title(f'PCA of GBM Samples (SV Features + Purity + {title_vaf})')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'pca_clustering.png'), dpi=300, bbox_inches='tight')
        plt.close()

        print(f"PCA Explained Variance: PC1={pca.explained_variance_ratio_[0]:.2%}, PC2={pca.explained_variance_ratio_[1]:.2%}")

    def plot_sv_stacked_histogram(self):
        """
        Create AGGREGATE stacked histogram of SV types (DEL, DUP, INV, BND, INS)
        Shows total counts across ALL samples combined
        """
        print("  Creating aggregate SV type stacked histogram...")

        fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        sv_types = ['DEL', 'DUP', 'INV', 'BND', 'INS']
        sv_types_hq = ['DEL_HQ', 'DUP_HQ', 'INV_HQ', 'BND_HQ', 'INS_HQ']
        colors = plt.cm.Set3(np.linspace(0, 1, len(sv_types)))

        # Calculate total counts across all samples
        total_counts = [self.df[sv_type].sum() for sv_type in sv_types]
        total_counts_hq = [self.df[sv_type_hq].sum() for sv_type_hq in sv_types_hq]

        # Plot 1: All SVs - aggregate stacked bar
        bottom = 0
        for i, (sv_type, count) in enumerate(zip(sv_types, total_counts)):
            axes[0].bar(0, count, bottom=bottom, label=sv_type,
                       alpha=0.8, color=colors[i], edgecolor='black', linewidth=1.5, width=0.6)
            # Add count label in the middle of each segment
            if count > 1000:  # Only label if segment is large enough
                axes[0].text(0, bottom + count/2, f'{count:,}',
                           ha='center', va='center', fontsize=11, fontweight='bold')
            bottom += count

        axes[0].set_xlim(-0.5, 0.5)
        axes[0].set_ylabel('Total SV Count (All Samples)', fontsize=12, fontweight='bold')
        axes[0].set_title('Aggregate Stacked Histogram:\nAll Structural Variants by Type',
                         fontsize=14, fontweight='bold')
        axes[0].set_xticks([])
        axes[0].legend(title='SV Type', fontsize=11, title_fontsize=12, loc='upper right')
        axes[0].grid(True, alpha=0.3, axis='y')
        axes[0].set_ylim(0, bottom * 1.02)

        # Plot 2: High-Quality SVs - aggregate stacked bar
        bottom_hq = 0
        for i, (sv_type, count_hq) in enumerate(zip(sv_types, total_counts_hq)):
            axes[1].bar(0, count_hq, bottom=bottom_hq, label=sv_type,
                       alpha=0.8, color=colors[i], edgecolor='black', linewidth=1.5, width=0.6)
            # Add count label in the middle of each segment
            if count_hq > 1000:  # Only label if segment is large enough
                axes[1].text(0, bottom_hq + count_hq/2, f'{count_hq:,}',
                           ha='center', va='center', fontsize=11, fontweight='bold')
            bottom_hq += count_hq

        axes[1].set_xlim(-0.5, 0.5)
        axes[1].set_ylabel('Total High-Quality SV Count (All Samples)', fontsize=12, fontweight='bold')
        axes[1].set_title('Aggregate Stacked Histogram:\nHigh-Quality Structural Variants\n(PASS + PRECISE + Support ≥10)',
                         fontsize=14, fontweight='bold')
        axes[1].set_xticks([])
        axes[1].legend(title='SV Type', fontsize=11, title_fontsize=12, loc='upper right')
        axes[1].grid(True, alpha=0.3, axis='y')
        axes[1].set_ylim(0, bottom_hq * 1.02)

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'sv_stacked_histogram.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()

        # Print summary to console
        print(f"\n  Aggregate SV Counts (All Samples):")
        print(f"    Total All SVs: {sum(total_counts):,}")
        for sv_type, count in zip(sv_types, total_counts):
            print(f"      {sv_type}: {count:,}")
        print(f"\n  Aggregate High-Quality SV Counts (All Samples):")
        print(f"    Total HQ SVs: {sum(total_counts_hq):,}")
        for sv_type, count_hq in zip(sv_types, total_counts_hq):
            print(f"      {sv_type}_HQ: {count_hq:,}")

    def plot_comprehensive_vaf_analysis(self):
        """
        Create comprehensive VAF analysis plots comparing different VAF metrics
        """
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))

        # Plot 1: All VAF vs High-Quality VAF
        axes[0, 0].scatter(self.df['Mean_VAF_All'], self.df['Mean_VAF_HQ'],
                          s=100, alpha=0.7, edgecolors='black')
        for i, sample in enumerate(self.df['Sample']):
            axes[0, 0].annotate(sample,
                              (self.df.iloc[i]['Mean_VAF_All'], self.df.iloc[i]['Mean_VAF_HQ']),
                              xytext=(3, 3), textcoords='offset points', fontsize=7)
        axes[0, 0].plot([0, 1], [0, 1], 'r--', alpha=0.5, label='y=x')
        axes[0, 0].set_xlabel('Mean VAF (All Variants)')
        axes[0, 0].set_ylabel('Mean VAF (High Quality)')
        axes[0, 0].set_title('All VAF vs High-Quality VAF')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)

        # Plot 2: Clonal Fraction vs Purity
        axes[0, 1].scatter(self.df['Purity'], self.df['Clonal_Fraction'],
                          s=100, alpha=0.7, c=self.df['Mean_VAF_Clonal'],
                          cmap='RdYlGn', edgecolors='black')
        for i, sample in enumerate(self.df['Sample']):
            axes[0, 1].annotate(sample,
                              (self.df.iloc[i]['Purity'], self.df.iloc[i]['Clonal_Fraction']),
                              xytext=(3, 3), textcoords='offset points', fontsize=7)
        axes[0, 1].set_xlabel('Tumor Purity')
        axes[0, 1].set_ylabel('Clonal Fraction (VAF ≥ 0.3)')
        axes[0, 1].set_title('Clonal Fraction vs Tumor Purity')
        axes[0, 1].grid(True, alpha=0.3)
        cbar = plt.colorbar(axes[0, 1].scatter(self.df['Purity'], self.df['Clonal_Fraction'],
                                              c=self.df['Mean_VAF_Clonal'], cmap='RdYlGn'),
                           ax=axes[0, 1])
        cbar.set_label('Mean Clonal VAF')

        # Plot 3: VAF Metrics Comparison
        vaf_metrics = ['Mean_VAF_All', 'Mean_VAF_HQ', 'Mean_VAF_Clonal', 'Weighted_VAF', 'Top10_VAF']
        vaf_df = self.df[vaf_metrics].copy()
        vaf_df.index = self.df['Sample']
        vaf_df.T.plot(kind='bar', ax=axes[0, 2], alpha=0.7)
        axes[0, 2].set_xlabel('VAF Metric')
        axes[0, 2].set_ylabel('VAF Value')
        axes[0, 2].set_title('Comparison of VAF Metrics')
        axes[0, 2].legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
        axes[0, 2].tick_params(axis='x', rotation=45)
        axes[0, 2].grid(True, alpha=0.3, axis='y')

        # Plot 4: Clonal vs Subclonal Counts
        x = np.arange(len(self.df))
        width = 0.35
        axes[1, 0].bar(x - width/2, self.df['Clonal_Count'], width,
                      label='Clonal (VAF ≥ 0.3)', alpha=0.7, color='darkgreen')
        axes[1, 0].bar(x + width/2, self.df['Subclonal_Count'], width,
                      label='Subclonal (0.1 ≤ VAF < 0.3)', alpha=0.7, color='orange')
        axes[1, 0].set_xlabel('Sample')
        axes[1, 0].set_ylabel('Variant Count')
        axes[1, 0].set_title('Clonal vs Subclonal Variant Counts')
        axes[1, 0].set_xticks(x)
        axes[1, 0].set_xticklabels(self.df['Sample'], rotation=45, ha='right')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3, axis='y')

        # Plot 5: High-Quality Variant Counts by Type
        hq_cols = ['DEL_HQ', 'DUP_HQ', 'INV_HQ', 'BND_HQ', 'INS_HQ']
        self.df.set_index('Sample')[hq_cols].plot(kind='bar', stacked=True,
                                                   ax=axes[1, 1], colormap='Set3')
        axes[1, 1].set_xlabel('Sample')
        axes[1, 1].set_ylabel('High-Quality SV Count')
        axes[1, 1].set_title('High-Quality SV Distribution (PASS + PRECISE + Support ≥ 10)')
        axes[1, 1].legend(title='SV Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        axes[1, 1].tick_params(axis='x', rotation=45)
        axes[1, 1].grid(True, alpha=0.3, axis='y')

        # Plot 6: Top Percentile VAF Comparison
        axes[1, 2].scatter(self.df['Top25_VAF'], self.df['Top10_VAF'],
                          s=100, alpha=0.7, c=self.df['Purity'],
                          cmap='viridis', edgecolors='black')
        for i, sample in enumerate(self.df['Sample']):
            axes[1, 2].annotate(sample,
                              (self.df.iloc[i]['Top25_VAF'], self.df.iloc[i]['Top10_VAF']),
                              xytext=(3, 3), textcoords='offset points', fontsize=7)
        axes[1, 2].set_xlabel('Top 25% VAF')
        axes[1, 2].set_ylabel('Top 10% VAF')
        axes[1, 2].set_title('Top Percentile VAF Comparison')
        axes[1, 2].grid(True, alpha=0.3)
        cbar2 = plt.colorbar(axes[1, 2].scatter(self.df['Top25_VAF'], self.df['Top10_VAF'],
                                               c=self.df['Purity'], cmap='viridis'),
                            ax=axes[1, 2])
        cbar2.set_label('Purity')

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'comprehensive_vaf_analysis.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()

    def create_heatmap(self):
        """Create heatmap of SV counts and features using high-quality VAF"""
        # Prepare data for heatmap
        feature_cols = ['Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS',
                       'Mean_VAF_HQ', 'Median_VAF_HQ', 'Clonal_Fraction']
        heatmap_data = self.df[feature_cols].copy()
        heatmap_data.index = self.df['Sample']

        # Normalize each column for better visualization
        heatmap_data_norm = (heatmap_data - heatmap_data.min()) / (heatmap_data.max() - heatmap_data.min())

        # Create heatmap
        plt.figure(figsize=(12, 8))
        sns.heatmap(heatmap_data_norm.T, annot=heatmap_data.T.values,
                   fmt='.2f', cmap='YlOrRd', cbar_kws={'label': 'Normalized Value'},
                   linewidths=0.5)
        plt.title('Heatmap of SV Features and High-Quality VAF Across GBM Samples')
        plt.xlabel('Sample')
        plt.ylabel('Feature')
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'sv_heatmap.png'), dpi=300, bbox_inches='tight')
        plt.close()

    def plot_sv_distribution(self):
        """Plot SV type distribution across samples"""
        sv_cols = ['DEL', 'DUP', 'INV', 'BND', 'INS']

        fig, axes = plt.subplots(2, 1, figsize=(12, 10))

        # Stacked bar plot
        self.df.set_index('Sample')[sv_cols].plot(kind='bar', stacked=True,
                                                   ax=axes[0], colormap='Set3')
        axes[0].set_title('Structural Variant Distribution by Sample (Stacked)')
        axes[0].set_xlabel('Sample')
        axes[0].set_ylabel('SV Count')
        axes[0].legend(title='SV Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        axes[0].tick_params(axis='x', rotation=45)

        # Grouped bar plot
        self.df.set_index('Sample')[sv_cols].plot(kind='bar', ax=axes[1], colormap='Set3')
        axes[1].set_title('Structural Variant Distribution by Sample (Grouped)')
        axes[1].set_xlabel('Sample')
        axes[1].set_ylabel('SV Count')
        axes[1].legend(title='SV Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        axes[1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'sv_distribution.png'), dpi=300, bbox_inches='tight')
        plt.close()

    def plot_vaf_analysis(self):
        """Plot high-quality VAF analysis"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # VAF vs Purity (using High-Quality VAF)
        axes[0].scatter(self.df['Purity'], self.df['Mean_VAF_HQ'], s=100, alpha=0.7, edgecolors='black')
        for i, sample in enumerate(self.df['Sample']):
            axes[0].annotate(sample,
                           (self.df.iloc[i]['Purity'], self.df.iloc[i]['Mean_VAF_HQ']),
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
        axes[0].set_xlabel('Tumor Purity')
        axes[0].set_ylabel('Mean VAF (High Quality)')
        axes[0].set_title('High-Quality Mean VAF vs Tumor Purity')
        axes[0].grid(True, alpha=0.3)

        # VAF distribution (High-Quality)
        axes[1].bar(self.df['Sample'], self.df['Mean_VAF_HQ'], alpha=0.7,
                   label='Mean VAF (HQ)', color='steelblue')
        axes[1].bar(self.df['Sample'], self.df['Median_VAF_HQ'], alpha=0.5,
                   label='Median VAF (HQ)', color='coral')
        axes[1].set_xlabel('Sample')
        axes[1].set_ylabel('VAF (High Quality)')
        axes[1].set_title('High-Quality VAF Distribution Across Samples')
        axes[1].tick_params(axis='x', rotation=45)
        axes[1].legend()
        axes[1].grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'vaf_analysis.png'), dpi=300, bbox_inches='tight')
        plt.close()

    def generate_report(self):
        """Generate comprehensive analysis report"""
        report_path = os.path.join(self.output_dir, 'analysis_report.txt')

        with open(report_path, 'w') as f:
            f.write("="*80 + "\n")
            f.write("STRUCTURAL VARIANT ANALYSIS REPORT - GBM SAMPLES\n")
            f.write("="*80 + "\n\n")

            f.write("1. SAMPLE SUMMARY\n")
            f.write("-"*80 + "\n")
            f.write(f"Total samples analyzed: {len(self.df)}\n")
            f.write(f"Sample IDs: {', '.join(self.df['Sample'].tolist())}\n\n")

            f.write("2. SV COUNTS SUMMARY\n")
            f.write("-"*80 + "\n")
            for sv_type in ['DEL', 'DUP', 'INV', 'BND', 'INS']:
                total = self.df[sv_type].sum()
                mean = self.df[sv_type].mean()
                f.write(f"{sv_type}: Total={total}, Mean={mean:.1f}\n")
            f.write("\n")

            f.write("3. PURITY AND VAF SUMMARY\n")
            f.write("-"*80 + "\n")
            f.write(f"Purity Range: {self.df['Purity'].min():.2f} - {self.df['Purity'].max():.2f}\n")
            f.write(f"Mean Purity: {self.df['Purity'].mean():.2f}\n")
            f.write(f"Mean VAF (All) Range: {self.df['Mean_VAF_All'].min():.3f} - {self.df['Mean_VAF_All'].max():.3f}\n")
            f.write(f"Mean VAF (HQ) Range: {self.df['Mean_VAF_HQ'].min():.3f} - {self.df['Mean_VAF_HQ'].max():.3f}\n")
            f.write(f"Overall Mean VAF (HQ): {self.df['Mean_VAF_HQ'].mean():.3f}\n")
            f.write(f"Average Clonal Fraction: {self.df['Clonal_Fraction'].mean():.3f}\n\n")

            if 'Cluster' in self.df.columns:
                f.write("4. CLUSTERING RESULTS\n")
                f.write("-"*80 + "\n")
                for cluster_id in sorted(self.df['Cluster'].unique()):
                    cluster_df = self.df[self.df['Cluster'] == cluster_id]
                    f.write(f"\nCluster {cluster_id}:\n")
                    f.write(f"  Samples: {', '.join(cluster_df['Sample'].tolist())}\n")
                    f.write(f"  Mean Purity: {cluster_df['Purity'].mean():.2f}\n")
                    f.write(f"  Mean Total SVs: {cluster_df[['DEL', 'DUP', 'INV', 'BND', 'INS']].sum(axis=1).mean():.1f}\n")
                    f.write(f"  Mean VAF (HQ): {cluster_df['Mean_VAF_HQ'].mean():.3f}\n")
                    f.write(f"  Mean Clonal Fraction: {cluster_df['Clonal_Fraction'].mean():.3f}\n")

            f.write("\n" + "="*80 + "\n")
            f.write("Analysis complete. See generated plots and tables for details.\n")
            f.write("="*80 + "\n")

        print(f"\nReport saved to: {report_path}")


def main():
    """Main analysis pipeline with SIGNIFICANT VAF focus"""
    # Use paths from CONFIGURATION SECTION at top of script
    sample_file = SAMPLE_FILE
    vcf_dir = VCF_DIR
    output_dir = OUTPUT_DIR

    print("="*80)
    print("STRUCTURAL VARIANT ANALYSIS FOR GBM SAMPLES (VERSION 2)")
    print("FOCUS: SIGNIFICANT VAF METRICS ONLY")
    print("="*80)
    print(f"Sample file: {sample_file}")
    print(f"VCF directory: {vcf_dir}")
    print(f"Output directory: {output_dir}")
    print("="*80 + "\n")

    print("VERSION 2 IMPROVEMENTS:")
    print("  - High-quality variant filtering (PASS + PRECISE + Support ≥ 10)")
    print("  - Clonal VAF analysis (VAF ≥ 0.3)")
    print("  - Weighted VAF (emphasizes high-VAF variants)")
    print("  - Top percentile VAF (top 10%, 25%)")
    print("  - Multiple clustering strategies")
    print("="*80 + "\n")

    # Initialize analyzer
    analyzer = SVAnalyzer(sample_file, vcf_dir, output_dir)

    # Step 1: Analyze all samples
    print("Step 1: Analyzing VCF files with quality filtering...")
    df = analyzer.analyze_all_samples()

    # Step 2: Create summary table
    print("\nStep 2: Creating summary table...")
    analyzer.save_summary_table()

    # Step 3: Perform clustering with DIFFERENT VAF STRATEGIES
    print("\nStep 3: Performing hierarchical clustering with multiple VAF strategies...")

    # Strategy 1: High-Quality VAF (RECOMMENDED)
    print("\n  3a. Clustering with HIGH-QUALITY VAF (RECOMMENDED)...")
    analyzer.perform_clustering(method='ward', n_clusters=3, use_strategy='high_quality')

    # Strategy 2: Clonal VAF
    print("\n  3b. Clustering with CLONAL VAF (VAF ≥ 0.3)...")
    analyzer.perform_clustering(method='ward', n_clusters=3, use_strategy='clonal')

    # Strategy 3: Weighted VAF
    print("\n  3c. Clustering with WEIGHTED VAF...")
    analyzer.perform_clustering(method='ward', n_clusters=3, use_strategy='weighted')

    # Step 4: PCA analysis
    print("\nStep 4: Performing PCA analysis...")
    analyzer.perform_pca_clustering(use_hq_vaf=True)

    # Step 5: Create visualizations
    print("\nStep 5: Creating visualizations...")
    analyzer.create_heatmap()
    analyzer.plot_sv_distribution()
    analyzer.plot_vaf_analysis()

    # Step 6: NEW - SV Type Stacked Histogram
    print("\nStep 6: Creating SV type stacked histogram...")
    analyzer.plot_sv_stacked_histogram()

    # Step 7: NEW - Comprehensive VAF Analysis
    print("\nStep 7: Creating comprehensive VAF analysis...")
    analyzer.plot_comprehensive_vaf_analysis()

    # Step 8: Generate report
    print("\nStep 8: Generating analysis report...")
    analyzer.generate_report()

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print(f"All results saved to: {output_dir}")
    print("="*80)
    print("\nGenerated files:")
    print("  - sv_summary_table.csv (full table with ALL VAF metrics)")
    print("  - sv_summary_display.csv (display table)")
    print("  - clustering_results.csv (samples with cluster assignments)")
    print("\n  Clustering dendrograms (3 strategies):")
    print("    - clustering_dendrogram_high_quality.png (RECOMMENDED)")
    print("    - clustering_dendrogram_clonal.png")
    print("    - clustering_dendrogram_weighted.png")
    print("\n  Visualizations:")
    print("    - pca_clustering.png")
    print("    - sv_heatmap.png")
    print("    - sv_distribution.png")
    print("    - vaf_analysis.png")
    print("    - sv_stacked_histogram.png (NEW - stacked SV type histogram)")
    print("    - comprehensive_vaf_analysis.png (NEW - compares VAF metrics)")
    print("\n  Report:")
    print("    - analysis_report.txt")
    print("="*80)
    print("\nRECOMMENDATION: Use 'high_quality' VAF metrics for clustering")
    print("  - Filters for PASS, PRECISE variants with Support ≥ 10")
    print("  - Focuses on reliable, significant variants")
    print("  - Includes clonal fraction for tumor evolution analysis")
    print("="*80)


if __name__ == "__main__":
    main()
