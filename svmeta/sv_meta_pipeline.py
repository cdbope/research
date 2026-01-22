#!/usr/bin/env python3
"""
Structural Variant Meta-Analysis Pipeline for GBM Cohort
Comprehensive analysis of 100 GBM samples with SV recurrence, hotspots, and patterns

WORKFLOW:
1. Normalize & merge SVs across all samples (SURVIVOR)
2. Build sample × SV-event matrix
3. Identify recurrent SV hotspots & genes
4. Cohort-wide pattern analysis (SV signatures)
5. Clinical & molecular associations
6. Pathway-level enrichment
7. External dataset comparison
"""

import os
import sys
import gzip
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE PATHS FOR YOUR SYSTEM
# ============================================================================

# Base directory for meta-analysis
BASE_DIR = "/home/chbope/extension/script/sv_meta_analysis"

# Directory containing normalized VCF files (*.norm.vcf.gz)
VCF_DIR = "/path/to/normalized/vcf/files"

# Reference genome FASTA (for normalization)
REF_FASTA = "/path/to/reference/genome.fa"

# Sample metadata file (sample_id, purity, IDH_status, MGMT_status, age, survival, etc.)
# Format: tab-separated with header
SAMPLE_METADATA = os.path.join(BASE_DIR, "sample_metadata.txt")

# Gene annotation file (BED format: chr, start, end, gene_name, strand)
GENE_ANNOTATION = os.path.join(BASE_DIR, "gencode_genes.bed")

# Known GBM driver genes list
GBM_DRIVERS = os.path.join(BASE_DIR, "gbm_driver_genes.txt")

# Output directory for all results
OUTPUT_DIR = os.path.join(BASE_DIR, "results")

# SURVIVOR binary path (download from: https://github.com/fritzsedlazeck/SURVIVOR)
SURVIVOR_BIN = "SURVIVOR"  # If in PATH, otherwise full path

# bcftools binary path
BCFTOOLS_BIN = "bcftools"  # If in PATH, otherwise full path

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================

# SURVIVOR merge parameters
SURVIVOR_MAX_DIST = 1000  # Maximum distance between SVs to merge (bp)
SURVIVOR_MIN_SUPPORT = 1  # Minimum number of samples supporting merged SV
SURVIVOR_TYPE_MATCH = 1   # Require same SV type (1=yes, 0=no)
SURVIVOR_STRAND_MATCH = 0 # Require same strand (1=yes, 0=no)
SURVIVOR_SIZE_SIMILARITY = 0.3  # Size similarity threshold (0-1)

# Recurrence thresholds
MIN_RECURRENCE = 3  # Minimum number of samples for "recurrent" SV
HOTSPOT_WINDOW = 50000  # Window size for defining hotspots (bp)
MIN_HOTSPOT_SAMPLES = 5  # Minimum samples per hotspot

# SV size bins for pattern analysis (bp)
SIZE_BINS = {
    'tiny': (50, 1000),
    'small': (1000, 10000),
    'medium': (10000, 100000),
    'large': (100000, 1000000),
    'xlarge': (1000000, float('inf'))
}

# Quality filters
MIN_SUPPORT_READS = 10
MIN_VAF = 0.1

# ============================================================================

class SVMetaAnalyzer:
    """Main class for SV meta-analysis"""

    def __init__(self):
        """Initialize the meta-analyzer"""
        self.base_dir = BASE_DIR
        self.vcf_dir = VCF_DIR
        self.ref_fasta = REF_FASTA
        self.output_dir = OUTPUT_DIR

        # Create output directories
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "normalized_vcfs"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "merged"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "matrices"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "hotspots"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "genes"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "patterns"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "plots"), exist_ok=True)

        # Load metadata
        self.sample_metadata = self.load_metadata()
        self.gene_annotation = self.load_gene_annotation()
        self.gbm_drivers = self.load_gbm_drivers()

        print(f"Initialized SVMetaAnalyzer")
        print(f"  Samples: {len(self.sample_metadata)}")
        print(f"  Genes annotated: {len(self.gene_annotation)}")
        print(f"  GBM drivers: {len(self.gbm_drivers)}")

    def load_metadata(self):
        """Load sample metadata"""
        if not os.path.exists(SAMPLE_METADATA):
            print(f"Warning: Metadata file not found: {SAMPLE_METADATA}")
            return pd.DataFrame()

        df = pd.read_csv(SAMPLE_METADATA, sep='\t')
        print(f"Loaded metadata for {len(df)} samples")
        return df

    def load_gene_annotation(self):
        """Load gene annotation from BED file"""
        if not os.path.exists(GENE_ANNOTATION):
            print(f"Warning: Gene annotation not found: {GENE_ANNOTATION}")
            return pd.DataFrame()

        df = pd.read_csv(GENE_ANNOTATION, sep='\t',
                        names=['chr', 'start', 'end', 'gene', 'strand'],
                        comment='#')
        return df

    def load_gbm_drivers(self):
        """Load list of known GBM driver genes"""
        if not os.path.exists(GBM_DRIVERS):
            # Default GBM drivers if file not found
            return {'EGFR', 'PTEN', 'NF1', 'CDKN2A', 'CDKN2B', 'TP53',
                   'PDGFRA', 'MET', 'PIK3CA', 'PIK3R1', 'RB1', 'TERT',
                   'IDH1', 'ATRX', 'MDM2', 'MDM4', 'CDK4', 'CDK6'}

        with open(GBM_DRIVERS) as f:
            genes = {line.strip() for line in f if line.strip()}
        return genes

    # ========================================================================
    # STEP 1: NORMALIZATION & MERGING
    # ========================================================================

    def normalize_vcf(self, input_vcf, output_vcf):
        """
        Normalize VCF using bcftools norm
        - Left-normalize indels
        - Sort variants
        - Apply quality filters
        """
        try:
            # Normalize and sort
            cmd = f"{BCFTOOLS_BIN} norm -f {self.ref_fasta} {input_vcf} | " \
                  f"{BCFTOOLS_BIN} sort -Oz -o {output_vcf}"
            subprocess.run(cmd, shell=True, check=True, capture_output=True)

            # Index
            subprocess.run(f"{BCFTOOLS_BIN} index {output_vcf}",
                          shell=True, check=True, capture_output=True)

            return True
        except subprocess.CalledProcessError as e:
            print(f"Error normalizing {input_vcf}: {e}")
            return False

    def create_vcf_list(self, vcf_files):
        """Create file listing all VCF files for SURVIVOR"""
        list_file = os.path.join(self.output_dir, "merged", "vcf_list.txt")
        with open(list_file, 'w') as f:
            for vcf in vcf_files:
                f.write(f"{vcf}\n")
        return list_file

    def merge_with_survivor(self, vcf_list_file):
        """
        Merge VCFs using SURVIVOR
        Creates a cohort-level merged SV catalog
        """
        output_vcf = os.path.join(self.output_dir, "merged", "merged_SV.vcf")

        cmd = [
            SURVIVOR_BIN, "merge",
            vcf_list_file,
            str(SURVIVOR_MAX_DIST),
            str(SURVIVOR_MIN_SUPPORT),
            str(SURVIVOR_TYPE_MATCH),
            str(SURVIVOR_STRAND_MATCH),
            str(SURVIVOR_SIZE_SIMILARITY),
            "50",  # Minimum size
            output_vcf
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True)
            print(f"✓ SURVIVOR merge complete: {output_vcf}")
            return output_vcf
        except subprocess.CalledProcessError as e:
            print(f"Error running SURVIVOR: {e}")
            return None

    # ========================================================================
    # STEP 2: BUILD SAMPLE × SV-EVENT MATRIX
    # ========================================================================

    def parse_merged_vcf(self, merged_vcf):
        """
        Parse merged VCF and extract:
        - SV events (id, chr, pos, end, type, size)
        - Sample support (which samples have which SVs)
        """
        events = []
        sample_support = []

        open_func = gzip.open if merged_vcf.endswith('.gz') else open

        with open_func(merged_vcf, 'rt') as f:
            samples = []
            for line in f:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    # Get sample names
                    fields = line.strip().split('\t')
                    samples = fields[9:]  # Samples start at column 9
                    continue

                fields = line.strip().split('\t')
                chrom = fields[0]
                pos = int(fields[1])
                svid = fields[2]
                info = fields[7]

                # Parse INFO field
                svtype = None
                svlen = 0
                end = pos

                for item in info.split(';'):
                    if item.startswith('SVTYPE='):
                        svtype = item.split('=')[1]
                    elif item.startswith('SVLEN='):
                        svlen = abs(int(item.split('=')[1]))
                    elif item.startswith('END='):
                        end = int(item.split('=')[1])

                # Get sample support
                sample_genotypes = fields[9:]
                supporting_samples = []
                for i, gt in enumerate(sample_genotypes):
                    # If genotype is not 0/0 or ./., sample has this SV
                    if gt and not gt.startswith('0/0') and not gt.startswith('./.'):
                        supporting_samples.append(samples[i])

                events.append({
                    'sv_id': svid,
                    'chr': chrom,
                    'pos': pos,
                    'end': end,
                    'svtype': svtype,
                    'svlen': svlen,
                    'n_samples': len(supporting_samples)
                })

                sample_support.append({
                    'sv_id': svid,
                    'samples': supporting_samples
                })

        events_df = pd.DataFrame(events)
        return events_df, sample_support, samples

    def build_presence_matrix(self, events_df, sample_support, samples):
        """
        Build sample × SV-event binary matrix
        Rows: samples
        Columns: SV events
        Values: 0/1 (presence/absence)
        """
        # Create matrix
        matrix = pd.DataFrame(0, index=samples, columns=events_df['sv_id'])

        for support in sample_support:
            sv_id = support['sv_id']
            for sample in support['samples']:
                if sample in matrix.index:
                    matrix.loc[sample, sv_id] = 1

        # Save matrix
        matrix_file = os.path.join(self.output_dir, "matrices", "sample_by_sv_matrix.csv")
        matrix.to_csv(matrix_file)
        print(f"✓ Sample × SV matrix saved: {matrix_file}")

        return matrix

    def build_sample_summary(self, matrix, events_df):
        """
        Create sample-level summary features
        - Total SV count
        - Counts by type
        - Number of recurrent SVs
        - Number of SVs in GBM genes
        """
        summary = []

        for sample in matrix.index:
            sample_svs = matrix.loc[sample]
            sv_ids = sample_svs[sample_svs == 1].index.tolist()

            # Get SV types for this sample
            sample_events = events_df[events_df['sv_id'].isin(sv_ids)]

            summary.append({
                'sample': sample,
                'total_svs': len(sv_ids),
                'del_count': len(sample_events[sample_events['svtype'] == 'DEL']),
                'dup_count': len(sample_events[sample_events['svtype'] == 'DUP']),
                'inv_count': len(sample_events[sample_events['svtype'] == 'INV']),
                'bnd_count': len(sample_events[sample_events['svtype'] == 'BND']),
                'ins_count': len(sample_events[sample_events['svtype'] == 'INS']),
                'recurrent_svs': len(sample_events[sample_events['n_samples'] >= MIN_RECURRENCE])
            })

        summary_df = pd.DataFrame(summary)
        summary_file = os.path.join(self.output_dir, "matrices", "sample_summary.csv")
        summary_df.to_csv(summary_file, index=False)
        print(f"✓ Sample summary saved: {summary_file}")

        return summary_df

    # ========================================================================
    # STEP 3: IDENTIFY RECURRENT SV HOTSPOTS & GENES
    # ========================================================================

    def identify_hotspots(self, events_df):
        """
        Identify genomic hotspots where many samples share SVs
        - Sliding window approach
        - Aggregate breakpoints by chromosome and position
        """
        hotspots = []

        # Process each chromosome
        for chrom in events_df['chr'].unique():
            chrom_events = events_df[events_df['chr'] == chrom].sort_values('pos')

            # Sliding window
            positions = chrom_events['pos'].values
            n_samples_list = chrom_events['n_samples'].values

            i = 0
            while i < len(positions):
                window_start = positions[i]
                window_end = window_start + HOTSPOT_WINDOW

                # Find all SVs in this window
                in_window = (positions >= window_start) & (positions < window_end)
                window_samples = n_samples_list[in_window]

                # If enough samples, it's a hotspot
                total_samples = np.sum(window_samples)
                if total_samples >= MIN_HOTSPOT_SAMPLES:
                    hotspots.append({
                        'chr': chrom,
                        'start': window_start,
                        'end': window_end,
                        'n_svs': np.sum(in_window),
                        'total_samples': total_samples,
                        'mean_samples_per_sv': np.mean(window_samples)
                    })

                # Move window
                i += 1

        hotspots_df = pd.DataFrame(hotspots)
        hotspots_file = os.path.join(self.output_dir, "hotspots", "sv_hotspots.csv")
        hotspots_df.to_csv(hotspots_file, index=False)
        print(f"✓ Identified {len(hotspots_df)} hotspots")

        return hotspots_df

    def analyze_gene_recurrence(self, events_df, sample_support):
        """
        Gene-level recurrence analysis
        - Intersect SV breakpoints with genes
        - Count samples with SVs affecting each gene
        """
        gene_recurrence = defaultdict(lambda: {
            'n_samples': 0,
            'n_del': 0,
            'n_dup': 0,
            'n_inv': 0,
            'n_bnd': 0,
            'n_ins': 0,
            'is_driver': False,
            'samples': set()
        })

        # For each SV, check which genes it intersects
        for idx, sv in events_df.iterrows():
            # Find genes overlapping this SV
            genes_hit = self.gene_annotation[
                (self.gene_annotation['chr'] == sv['chr']) &
                (self.gene_annotation['start'] <= sv['end']) &
                (self.gene_annotation['end'] >= sv['pos'])
            ]

            # Get samples for this SV
            sv_samples = [s['samples'] for s in sample_support if s['sv_id'] == sv['sv_id']]
            if sv_samples:
                sv_samples = sv_samples[0]
            else:
                sv_samples = []

            # Update gene counts
            for gene in genes_hit['gene'].unique():
                gene_recurrence[gene]['n_samples'] += len(sv_samples)
                gene_recurrence[gene]['samples'].update(sv_samples)

                if sv['svtype'] == 'DEL':
                    gene_recurrence[gene]['n_del'] += 1
                elif sv['svtype'] == 'DUP':
                    gene_recurrence[gene]['n_dup'] += 1
                elif sv['svtype'] == 'INV':
                    gene_recurrence[gene]['n_inv'] += 1
                elif sv['svtype'] == 'BND':
                    gene_recurrence[gene]['n_bnd'] += 1
                elif sv['svtype'] == 'INS':
                    gene_recurrence[gene]['n_ins'] += 1

                if gene in self.gbm_drivers:
                    gene_recurrence[gene]['is_driver'] = True

        # Convert to DataFrame
        gene_data = []
        for gene, data in gene_recurrence.items():
            gene_data.append({
                'gene': gene,
                'n_samples_affected': len(data['samples']),
                'n_del': data['n_del'],
                'n_dup': data['n_dup'],
                'n_inv': data['n_inv'],
                'n_bnd': data['n_bnd'],
                'n_ins': data['n_ins'],
                'total_svs': sum([data['n_del'], data['n_dup'], data['n_inv'],
                                 data['n_bnd'], data['n_ins']]),
                'is_gbm_driver': data['is_driver']
            })

        gene_df = pd.DataFrame(gene_data).sort_values('n_samples_affected', ascending=False)
        gene_file = os.path.join(self.output_dir, "genes", "gene_recurrence.csv")
        gene_df.to_csv(gene_file, index=False)
        print(f"✓ Gene recurrence analysis complete: {len(gene_df)} genes affected")

        return gene_df

    # ========================================================================
    # STEP 4: SV PATTERN ANALYSIS
    # ========================================================================

    def compute_sv_features(self, sample_summary, events_df, matrix):
        """
        Compute SV-feature vectors for each sample
        - Size-binned counts
        - Intra vs inter-chromosomal
        - Clustered SVs
        """
        features = []

        for sample in matrix.index:
            sample_svs = matrix.loc[sample]
            sv_ids = sample_svs[sample_svs == 1].index.tolist()
            sample_events = events_df[events_df['sv_id'].isin(sv_ids)]

            feature_dict = {'sample': sample}

            # Size-binned counts by type
            for sv_type in ['DEL', 'DUP', 'INV', 'BND', 'INS']:
                type_events = sample_events[sample_events['svtype'] == sv_type]

                for bin_name, (min_size, max_size) in SIZE_BINS.items():
                    count = len(type_events[
                        (type_events['svlen'] >= min_size) &
                        (type_events['svlen'] < max_size)
                    ])
                    feature_dict[f'{sv_type}_{bin_name}'] = count

            # Inter-chromosomal translocations (BND only)
            feature_dict['inter_chrom_bnd'] = len(sample_events[
                (sample_events['svtype'] == 'BND')
            ])

            features.append(feature_dict)

        features_df = pd.DataFrame(features)
        features_file = os.path.join(self.output_dir, "patterns", "sv_features.csv")
        features_df.to_csv(features_file, index=False)
        print(f"✓ SV feature vectors computed")

        return features_df

    def cluster_samples_by_patterns(self, features_df):
        """
        Cluster samples based on SV patterns
        - PCA, UMAP, t-SNE
        - Hierarchical clustering
        """
        # Prepare data
        X = features_df.drop('sample', axis=1).values
        samples = features_df['sample'].values

        # Standardize
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # PCA
        pca = PCA(n_components=min(10, X_scaled.shape[1]))
        X_pca = pca.fit_transform(X_scaled)

        # UMAP
        reducer = umap.UMAP(n_components=2, random_state=42)
        X_umap = reducer.fit_transform(X_scaled)

        # t-SNE
        tsne = TSNE(n_components=2, random_state=42)
        X_tsne = tsne.fit_transform(X_scaled)

        # Save embeddings
        embeddings_df = pd.DataFrame({
            'sample': samples,
            'PC1': X_pca[:, 0],
            'PC2': X_pca[:, 1],
            'UMAP1': X_umap[:, 0],
            'UMAP2': X_umap[:, 1],
            'tSNE1': X_tsne[:, 0],
            'tSNE2': X_tsne[:, 1]
        })
        embeddings_file = os.path.join(self.output_dir, "patterns", "sample_embeddings.csv")
        embeddings_df.to_csv(embeddings_file, index=False)

        # Plot
        self.plot_dimensionality_reduction(embeddings_df)

        # Hierarchical clustering
        linkage_matrix = linkage(X_scaled, method='ward')

        # Plot dendrogram
        self.plot_dendrogram(linkage_matrix, samples)

        return embeddings_df, linkage_matrix

    def plot_dimensionality_reduction(self, embeddings_df):
        """Plot PCA, UMAP, t-SNE"""
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # PCA
        axes[0].scatter(embeddings_df['PC1'], embeddings_df['PC2'], alpha=0.6)
        axes[0].set_xlabel('PC1')
        axes[0].set_ylabel('PC2')
        axes[0].set_title('PCA')

        # UMAP
        axes[1].scatter(embeddings_df['UMAP1'], embeddings_df['UMAP2'], alpha=0.6)
        axes[1].set_xlabel('UMAP1')
        axes[1].set_ylabel('UMAP2')
        axes[1].set_title('UMAP')

        # t-SNE
        axes[2].scatter(embeddings_df['tSNE1'], embeddings_df['tSNE2'], alpha=0.6)
        axes[2].set_xlabel('t-SNE1')
        axes[2].set_ylabel('t-SNE2')
        axes[2].set_title('t-SNE')

        plt.tight_layout()
        plot_file = os.path.join(self.output_dir, "plots", "dimensionality_reduction.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Dimensionality reduction plot saved")

    def plot_dendrogram(self, linkage_matrix, samples):
        """Plot hierarchical clustering dendrogram"""
        fig, ax = plt.subplots(figsize=(15, 8))

        dendrogram(linkage_matrix, labels=samples, ax=ax, leaf_font_size=8)
        ax.set_xlabel('Sample')
        ax.set_ylabel('Distance')
        ax.set_title('Hierarchical Clustering of Samples by SV Patterns')

        plt.tight_layout()
        plot_file = os.path.join(self.output_dir, "plots", "clustering_dendrogram.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Clustering dendrogram saved")


def main():
    """Main meta-analysis pipeline"""
    print("="*80)
    print("SV META-ANALYSIS PIPELINE FOR GBM COHORT")
    print("="*80)
    print()

    # Initialize analyzer
    analyzer = SVMetaAnalyzer()

    # TODO: Add full pipeline steps here
    # This is a framework - users will need to implement specific steps

    print("\n" + "="*80)
    print("Pipeline framework created. See documentation for next steps.")
    print("="*80)


if __name__ == "__main__":
    main()
