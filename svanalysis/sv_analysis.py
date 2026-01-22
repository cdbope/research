#!/usr/bin/env python3
"""
Structural Variant Analysis for GBM Samples
Analyzes Sniffles VCF files for SV counts, VAF, and clustering
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

    def parse_vcf(self, vcf_path):
        """
        Parse a VCF file and extract SV information

        Args:
            vcf_path: Path to VCF file (gzipped)

        Returns:
            Dictionary with SV counts and VAF information
        """
        sv_counts = {'DEL': 0, 'DUP': 0, 'INV': 0, 'BND': 0, 'INS': 0}
        vaf_values = []

        open_func = gzip.open if vcf_path.endswith('.gz') else open

        with open_func(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue

                info = fields[7]

                # Extract SV type
                svtype = None
                af = None

                for item in info.split(';'):
                    if item.startswith('SVTYPE='):
                        svtype = item.split('=')[1]
                    elif item.startswith('AF='):
                        try:
                            af = float(item.split('=')[1])
                        except:
                            pass

                # Count SV types
                if svtype in sv_counts:
                    sv_counts[svtype] += 1

                # Collect VAF values
                if af is not None:
                    vaf_values.append(af)

        return {
            'sv_counts': sv_counts,
            'vaf_values': vaf_values,
            'mean_vaf': np.mean(vaf_values) if vaf_values else 0,
            'median_vaf': np.median(vaf_values) if vaf_values else 0
        }

    def analyze_all_samples(self):
        """Analyze all VCF files and create summary table"""
        results = []

        for sample_id, sample_info in self.samples.items():
            vcf_file = os.path.join(self.vcf_dir, f"{sample_id}.wf_sv.vcf.gz")

            if not os.path.exists(vcf_file):
                print(f"Warning: VCF file not found for {sample_id}")
                continue

            print(f"Processing {sample_id}...")
            vcf_data = self.parse_vcf(vcf_file)

            result = {
                'Sample': sample_id,
                'Purity': sample_info['purity'],
                'DEL': vcf_data['sv_counts']['DEL'],
                'DUP': vcf_data['sv_counts']['DUP'],
                'INV': vcf_data['sv_counts']['INV'],
                'BND': vcf_data['sv_counts']['BND'],
                'INS': vcf_data['sv_counts']['INS'],
                'Mean_VAF': vcf_data['mean_vaf'],
                'Median_VAF': vcf_data['median_vaf']
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

    def perform_clustering(self, method='ward', distance_metric='euclidean', n_clusters=None):
        """
        Perform hierarchical clustering optimized for GBM samples

        For GBM samples, we consider:
        - SV burden (total counts)
        - SV type distribution
        - Purity
        - VAF metrics (clonality indicators)

        Args:
            method: Linkage method ('ward', 'average', 'complete')
            distance_metric: Distance metric for clustering
            n_clusters: Number of clusters (if None, will be determined by dendrogram)
        """
        # Prepare features for clustering
        feature_cols = ['Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS', 'Mean_VAF', 'Median_VAF']
        X = self.df[feature_cols].values
        sample_names = self.df['Sample'].values

        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Perform hierarchical clustering
        linkage_matrix = linkage(X_scaled, method=method, metric=distance_metric)

        # Create dendrogram
        plt.figure(figsize=(12, 6))
        dendrogram(linkage_matrix, labels=sample_names, leaf_rotation=45, leaf_font_size=10)
        plt.title(f'Hierarchical Clustering of GBM Samples\n(Method: {method}, Metric: {distance_metric})')
        plt.xlabel('Sample')
        plt.ylabel('Distance')
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'clustering_dendrogram.png'), dpi=300, bbox_inches='tight')
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

    def perform_pca_clustering(self):
        """
        Perform PCA for dimensionality reduction and visualization
        Useful for understanding sample relationships in GBM cohort
        """
        feature_cols = ['Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS', 'Mean_VAF', 'Median_VAF']
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
        plt.title('PCA of GBM Samples (SV Features + Purity + VAF)')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'pca_clustering.png'), dpi=300, bbox_inches='tight')
        plt.close()

        print(f"PCA Explained Variance: PC1={pca.explained_variance_ratio_[0]:.2%}, PC2={pca.explained_variance_ratio_[1]:.2%}")

    def create_heatmap(self):
        """Create heatmap of SV counts and features"""
        # Prepare data for heatmap
        feature_cols = ['Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS', 'Mean_VAF', 'Median_VAF']
        heatmap_data = self.df[feature_cols].copy()
        heatmap_data.index = self.df['Sample']

        # Normalize each column for better visualization
        heatmap_data_norm = (heatmap_data - heatmap_data.min()) / (heatmap_data.max() - heatmap_data.min())

        # Create heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data_norm.T, annot=heatmap_data.T.values,
                   fmt='.2f', cmap='YlOrRd', cbar_kws={'label': 'Normalized Value'},
                   linewidths=0.5)
        plt.title('Heatmap of SV Features and VAF Across GBM Samples')
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
        """Plot VAF analysis"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # VAF vs Purity
        axes[0].scatter(self.df['Purity'], self.df['Mean_VAF'], s=100, alpha=0.7, edgecolors='black')
        for i, sample in enumerate(self.df['Sample']):
            axes[0].annotate(sample,
                           (self.df.iloc[i]['Purity'], self.df.iloc[i]['Mean_VAF']),
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
        axes[0].set_xlabel('Tumor Purity')
        axes[0].set_ylabel('Mean VAF')
        axes[0].set_title('Mean VAF vs Tumor Purity')
        axes[0].grid(True, alpha=0.3)

        # VAF distribution
        axes[1].bar(self.df['Sample'], self.df['Mean_VAF'], alpha=0.7, label='Mean VAF', color='steelblue')
        axes[1].bar(self.df['Sample'], self.df['Median_VAF'], alpha=0.5, label='Median VAF', color='coral')
        axes[1].set_xlabel('Sample')
        axes[1].set_ylabel('VAF')
        axes[1].set_title('VAF Distribution Across Samples')
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
            f.write(f"Mean VAF Range: {self.df['Mean_VAF'].min():.3f} - {self.df['Mean_VAF'].max():.3f}\n")
            f.write(f"Overall Mean VAF: {self.df['Mean_VAF'].mean():.3f}\n\n")

            if 'Cluster' in self.df.columns:
                f.write("4. CLUSTERING RESULTS\n")
                f.write("-"*80 + "\n")
                for cluster_id in sorted(self.df['Cluster'].unique()):
                    cluster_df = self.df[self.df['Cluster'] == cluster_id]
                    f.write(f"\nCluster {cluster_id}:\n")
                    f.write(f"  Samples: {', '.join(cluster_df['Sample'].tolist())}\n")
                    f.write(f"  Mean Purity: {cluster_df['Purity'].mean():.2f}\n")
                    f.write(f"  Mean Total SVs: {cluster_df[['DEL', 'DUP', 'INV', 'BND', 'INS']].sum(axis=1).mean():.1f}\n")
                    f.write(f"  Mean VAF: {cluster_df['Mean_VAF'].mean():.3f}\n")

            f.write("\n" + "="*80 + "\n")
            f.write("Analysis complete. See generated plots and tables for details.\n")
            f.write("="*80 + "\n")

        print(f"\nReport saved to: {report_path}")


def main():
    """Main analysis pipeline"""
    # Configuration
    base_dir = "/home/chbope/extension/script/svanalysis"
    sample_file = os.path.join(base_dir, "sample.txt")
    vcf_dir = base_dir
    output_dir = os.path.join(base_dir, "results")

    print("="*80)
    print("STRUCTURAL VARIANT ANALYSIS FOR GBM SAMPLES")
    print("="*80)
    print(f"Sample file: {sample_file}")
    print(f"VCF directory: {vcf_dir}")
    print(f"Output directory: {output_dir}")
    print("="*80 + "\n")

    # Initialize analyzer
    analyzer = SVAnalyzer(sample_file, vcf_dir, output_dir)

    # Step 1: Analyze all samples
    print("Step 1: Analyzing VCF files...")
    df = analyzer.analyze_all_samples()

    # Step 2: Create summary table
    print("\nStep 2: Creating summary table...")
    analyzer.save_summary_table()

    # Step 3: Perform clustering
    print("\nStep 3: Performing hierarchical clustering...")
    analyzer.perform_clustering(method='ward', n_clusters=3)

    # Step 4: PCA analysis
    print("\nStep 4: Performing PCA analysis...")
    analyzer.perform_pca_clustering()

    # Step 5: Create visualizations
    print("\nStep 5: Creating visualizations...")
    analyzer.create_heatmap()
    analyzer.plot_sv_distribution()
    analyzer.plot_vaf_analysis()

    # Step 6: Generate report
    print("\nStep 6: Generating analysis report...")
    analyzer.generate_report()

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print(f"All results saved to: {output_dir}")
    print("="*80)
    print("\nGenerated files:")
    print("  - sv_summary_table.csv (full table with VAF)")
    print("  - sv_summary_display.csv (display table)")
    print("  - clustering_results.csv (samples with cluster assignments)")
    print("  - clustering_dendrogram.png")
    print("  - pca_clustering.png")
    print("  - sv_heatmap.png")
    print("  - sv_distribution.png")
    print("  - vaf_analysis.png")
    print("  - analysis_report.txt")
    print("="*80)


if __name__ == "__main__":
    main()
