#!/usr/bin/env python3
"""
Step 3: Build sample × SV matrix and perform comprehensive analysis

This script:
1. Parses merged VCF
2. Builds presence/absence matrix
3. Identifies recurrent SVs and hotspots
4. Analyzes gene-level recurrence
5. Computes SV pattern features
6. Performs clustering and dimensionality reduction
"""

import os
import sys
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from collections import defaultdict
try:
    import umap
except ImportError:
    print("Warning: umap-learn not installed. UMAP analysis will be skipped.")
    umap = None
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

MERGED_VCF = "/home/chbope/extension/script/svmeta/results/merged/merged_SV.vcf.gz"
GENE_ANNOTATION = "/home/chbope/extension/nWGS_manuscript_data/data/reference/gencode.v48.annotation.gff3"  # UPDATE THIS
GBM_DRIVERS = "/home/chbope/extension/script/svmeta/gbm_driver_genes.txt"
OUTPUT_DIR = "/home/chbope/extension/script/svmeta/results"

MIN_RECURRENCE = 3  # Minimum samples for "recurrent" SV
HOTSPOT_WINDOW = 50000  # bp
MIN_HOTSPOT_SAMPLES = 5

# SV size bins (bp)
SIZE_BINS = {
    'tiny': (50, 1000),
    'small': (1000, 10000),
    'medium': (10000, 100000),
    'large': (100000, 1000000),
    'xlarge': (1000000, float('inf'))
}

# ============================================================================

def parse_merged_vcf(vcf_file):
    """Parse merged VCF and extract SV events and sample support"""
    events = []
    sample_matrix_data = {}

    # Always use gzip.open for .gz files, handle both regular and BGZF
    print(f"DEBUG: vcf_file = {vcf_file}")
    print(f"DEBUG: ends with .gz? {vcf_file.endswith('.gz')}")
    print(f"DEBUG: gzip.open = {gzip.open}")

    if vcf_file.endswith('.gz'):
        print(f"DEBUG: In if branch")
        open_func = gzip.open
        mode = 'rt'
        print(f"DEBUG: Set open_func to {open_func}")
    else:
        print(f"DEBUG: In else branch")
        open_func = open
        mode = 'r'

    print(f"DEBUG: Final open_func = {open_func}")
    print(f"  Opening: {vcf_file}")
    print(f"  Using: {open_func.__name__}")

    with open_func(vcf_file, mode) as f:
        samples = []
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                samples = fields[9:] if len(fields) > 9 else []
                # Initialize matrix
                for sample in samples:
                    sample_matrix_data[sample] = {}
                continue

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            svid = fields[2]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]

            # Parse INFO
            svtype = None
            svlen = 0
            end = pos
            support_samples = 0

            for item in info.split(';'):
                if item.startswith('SVTYPE='):
                    svtype = item.split('=')[1]
                elif item.startswith('SVLEN='):
                    try:
                        svlen = abs(int(item.split('=')[1]))
                    except:
                        pass
                elif item.startswith('END='):
                    try:
                        end = int(item.split('=')[1])
                    except:
                        pass
                elif item.startswith('SUPP='):
                    try:
                        support_samples = int(item.split('=')[1])
                    except:
                        pass

            # Get sample genotypes
            if len(fields) > 9:
                genotypes = fields[9:]
                n_supporting = 0
                for i, gt in enumerate(genotypes):
                    if samples and i < len(samples):
                        sample = samples[i]
                        # Mark presence if not 0/0 or ./.
                        if gt and not gt.startswith('0/0') and not gt.startswith('./.'):
                            sample_matrix_data[sample][svid] = 1
                            n_supporting += 1
                        else:
                            sample_matrix_data[sample][svid] = 0

                # Use actual count from genotypes
                support_samples = n_supporting

            events.append({
                'sv_id': svid,
                'chr': chrom,
                'pos': pos,
                'end': end,
                'svtype': svtype,
                'svlen': svlen if svlen > 0 else (end - pos),
                'n_samples': support_samples
            })

    events_df = pd.DataFrame(events)

    # Build matrix DataFrame
    if samples and sample_matrix_data:
        matrix_df = pd.DataFrame(sample_matrix_data).T
        matrix_df = matrix_df.fillna(0).astype(int)
    else:
        matrix_df = pd.DataFrame()

    return events_df, matrix_df, samples


def save_sample_summary(matrix_df, events_df, output_dir):
    """Create sample-level summary"""
    summary = []

    for sample in matrix_df.index:
        sample_svs = matrix_df.loc[sample]
        sv_ids = sample_svs[sample_svs == 1].index.tolist()
        sample_events = events_df[events_df['sv_id'].isin(sv_ids)]

        # Count by type
        type_counts = sample_events['svtype'].value_counts().to_dict()

        summary.append({
            'sample': sample,
            'total_svs': len(sv_ids),
            'DEL': type_counts.get('DEL', 0),
            'DUP': type_counts.get('DUP', 0),
            'INV': type_counts.get('INV', 0),
            'BND': type_counts.get('BND', 0),
            'INS': type_counts.get('INS', 0),
            'recurrent_svs': len(sample_events[sample_events['n_samples'] >= MIN_RECURRENCE])
        })

    summary_df = pd.DataFrame(summary)
    summary_file = os.path.join(output_dir, "matrices", "sample_summary.csv")
    summary_df.to_csv(summary_file, index=False)
    print(f"✓ Sample summary: {summary_file}")

    return summary_df


def identify_recurrent_svs(events_df, output_dir):
    """Identify recurrent SVs (present in multiple samples)"""
    recurrent = events_df[events_df['n_samples'] >= MIN_RECURRENCE].copy()
    recurrent = recurrent.sort_values('n_samples', ascending=False)

    recurrent_file = os.path.join(output_dir, "hotspots", "recurrent_svs.csv")
    recurrent.to_csv(recurrent_file, index=False)
    print(f"✓ Recurrent SVs ({len(recurrent)}): {recurrent_file}")

    return recurrent


def identify_hotspots(events_df, output_dir):
    """Identify genomic hotspots with many SVs"""
    hotspots = []

    for chrom in sorted(events_df['chr'].unique()):
        chrom_events = events_df[events_df['chr'] == chrom].sort_values('pos')

        # Sliding window
        for i in range(len(chrom_events)):
            window_start = chrom_events.iloc[i]['pos']
            window_end = window_start + HOTSPOT_WINDOW

            in_window = chrom_events[
                (chrom_events['pos'] >= window_start) &
                (chrom_events['pos'] < window_end)
            ]

            total_samples = in_window['n_samples'].sum()
            if total_samples >= MIN_HOTSPOT_SAMPLES:
                hotspots.append({
                    'chr': chrom,
                    'start': window_start,
                    'end': window_end,
                    'n_svs': len(in_window),
                    'total_samples': total_samples,
                    'mean_samples': in_window['n_samples'].mean()
                })

    if hotspots:
        hotspots_df = pd.DataFrame(hotspots)
        # Remove overlapping hotspots
        hotspots_df = hotspots_df.drop_duplicates(['chr', 'start'], keep='first')

        hotspots_file = os.path.join(output_dir, "hotspots", "sv_hotspots.csv")
        hotspots_df.to_csv(hotspots_file, index=False)
        print(f"✓ Hotspots ({len(hotspots_df)}): {hotspots_file}")

        return hotspots_df
    else:
        print("✓ No hotspots identified")
        return pd.DataFrame()




def analyze_genes(events_df, matrix_df, gene_file, driver_genes_file, output_dir):
    """Gene-level recurrence analysis (supports GFF3/GFF or BED-like files)."""
    if not os.path.exists(gene_file):
        print(f"⚠ Gene annotation not found: {gene_file}")
        print("  Skipping gene analysis")
        return None

    # ------------------------------------------------------------------
    # 1. Load gene annotation
    # ------------------------------------------------------------------
    if gene_file.endswith((".gff3", ".gff")):
        # GENCODE GFF3: seqid, source, type, start, end, score, strand, phase, attributes
        raw = pd.read_csv(
            gene_file,
            sep="\t",
            comment="#",
            header=None,
            dtype={0: str, 2: str}
        )

        # Keep only gene features
        genes_raw = raw[raw[2] == "gene"].copy()

        # Helper to parse attributes into a dict
        def parse_attrs(attr):
            d = {}
            for field in str(attr).split(";"):
                field = field.strip()
                if not field:
                    continue
                if "=" in field:
                    k, v = field.split("=", 1)
                    d[k] = v
            return d

        attrs = genes_raw[8].apply(parse_attrs)

        genes_df = pd.DataFrame()
        genes_df["chr"]    = genes_raw[0].astype(str)
        genes_df["start"]  = genes_raw[3].astype(int)
        genes_df["end"]    = genes_raw[4].astype(int)
        genes_df["strand"] = genes_raw[6].astype(str)

        # Extract Ensembl gene_id and HGNC symbol (gene_name)
        genes_df["gene_id"] = attrs.apply(
            lambda d: d.get("gene_id", d.get("ID", ""))
        )
        genes_df["gene_symbol"] = attrs.apply(
            lambda d: d.get("gene_name", d.get("Name", d.get("gene_id", d.get("ID", ""))))
        )

    else:
        # Assume BED-like with 6 columns: chr, start, end, gene_symbol, score, strand
        genes_df = pd.read_csv(
            gene_file,
            sep="\t",
            names=["chr", "start", "end", "gene_symbol", "score", "strand"],
            comment="#",
            dtype={"chr": str}
        )
        genes_df["start"] = pd.to_numeric(genes_df["start"], errors="coerce")
        genes_df["end"]   = pd.to_numeric(genes_df["end"], errors="coerce")
        genes_df = genes_df.dropna(subset=["start", "end"])
        genes_df["start"] = genes_df["start"].astype(int)
        genes_df["end"]   = genes_df["end"].astype(int)
        genes_df["gene_id"] = ""

    # ------------------------------------------------------------------
    # Normalize chr naming to match events_df
    # ------------------------------------------------------------------
    def norm_chr(x):
        x = str(x)
        return x[3:] if x.startswith("chr") else x

    genes_df["chr_norm"] = genes_df["chr"].apply(norm_chr)
    events_df = events_df.copy()
    events_df["chr_norm"] = events_df["chr"].astype(str).apply(norm_chr)

    # ------------------------------------------------------------------
    # 2. Load driver genes (assumed to be gene symbols)
    # ------------------------------------------------------------------
    driver_genes = set()
    if os.path.exists(driver_genes_file):
        with open(driver_genes_file) as f:
            driver_genes = {line.strip() for line in f if line.strip()}

    # ------------------------------------------------------------------
    # 3. Count SVs per gene
    # ------------------------------------------------------------------
    gene_sv_counts = defaultdict(lambda: {
        "n_svs": 0,
        "sv_ids": set(),      # set of sv_ids overlapping this gene
        "DEL": 0,
        "DUP": 0,
        "INV": 0,
        "BND": 0,
        "INS": 0,
        "gene_id": None,
    })

    for _, sv in events_df.iterrows():
        # Find overlapping genes: simple interval overlap
        overlapping = genes_df[
            (genes_df["chr_norm"] == sv["chr_norm"]) &
            (genes_df["start"] <= sv["end"]) &
            (genes_df["end"]   >= sv["pos"])
        ]

        for _, row in overlapping.iterrows():
            symbol = row["gene_symbol"]
            gid    = row.get("gene_id", "")

            gene_sv_counts[symbol]["n_svs"] += 1
            gene_sv_counts[symbol]["sv_ids"].add(sv["sv_id"])

            # Map SVTYPE to our tracked categories
            svtype = sv["svtype"]
            if svtype == "TRA":
                svtype = "BND"  # treat translocations as breakends

            if svtype in ["DEL", "DUP", "INV", "BND", "INS"]:
                gene_sv_counts[symbol][svtype] += 1

            gene_sv_counts[symbol]["gene_id"] = gid

    # ------------------------------------------------------------------
    # 4. Build output DataFrame
    # ------------------------------------------------------------------
    gene_data = []
    for symbol, data in gene_sv_counts.items():
        # Calculate number of unique SAMPLES (not SVs) with SVs in this gene
        sv_ids = list(data["sv_ids"])
        if len(matrix_df) > 0 and len(sv_ids) > 0:
            # Get columns for these SVs (if they exist in matrix)
            existing_sv_ids = [sv_id for sv_id in sv_ids if sv_id in matrix_df.columns]
            if existing_sv_ids:
                # Count samples where at least one of these SVs is present (value = 1)
                gene_sample_matrix = matrix_df[existing_sv_ids]
                n_unique_samples = (gene_sample_matrix.sum(axis=1) > 0).sum()
            else:
                n_unique_samples = 0
        else:
            n_unique_samples = 0

        gene_data.append({
            "gene": symbol,                   # main column = HGNC symbol
            "gene_id": data["gene_id"],       # Ensembl ID
            "n_svs": data["n_svs"],
            "n_samples": n_unique_samples,    # CORRECTED: unique samples, not unique SVs
            "DEL": data["DEL"],
            "DUP": data["DUP"],
            "INV": data["INV"],
            "BND": data["BND"],
            "INS": data["INS"],
            "is_driver": symbol in driver_genes,
        })

    if gene_data:
        gene_df = pd.DataFrame(gene_data).sort_values("n_svs", ascending=False)
        os.makedirs(os.path.join(output_dir, "genes"), exist_ok=True)
        gene_file_out = os.path.join(output_dir, "genes", "gene_sv_counts.csv")
        gene_df.to_csv(gene_file_out, index=False)
        print(f"✓ Gene analysis: {gene_file_out}")

        print("\n  Top 10 genes by SV count:")
        for _, row in gene_df.head(10).iterrows():
            driver_mark = " [DRIVER]" if row["is_driver"] else ""
            print(f"    {row['gene']} ({row['gene_id']}): {row['n_svs']} SVs{driver_mark}")

        return gene_df
    else:
        print("✓ No genes with SVs")
        return None





def compute_sv_features(matrix_df, events_df, output_dir):
    """Compute SV feature vectors per sample"""
    features = []

    for sample in matrix_df.index:
        sample_svs = matrix_df.loc[sample]
        sv_ids = sample_svs[sample_svs == 1].index.tolist()
        sample_events = events_df[events_df['sv_id'].isin(sv_ids)]

        feature_dict = {'sample': sample}

        # Size-binned counts by type
        for svtype in ['DEL', 'DUP', 'INV', 'BND', 'INS']:
            type_events = sample_events[sample_events['svtype'] == svtype]

            for bin_name, (min_size, max_size) in SIZE_BINS.items():
                count = len(type_events[
                    (type_events['svlen'] >= min_size) &
                    (type_events['svlen'] < max_size)
                ])
                feature_dict[f'{svtype}_{bin_name}'] = count

        features.append(feature_dict)

    features_df = pd.DataFrame(features)
    features_file = os.path.join(output_dir, "patterns", "sv_features.csv")
    features_df.to_csv(features_file, index=False)
    print(f"✓ SV features: {features_file}")

    return features_df


def perform_clustering(features_df, output_dir):
    """Cluster samples by SV patterns"""
    X = features_df.drop('sample', axis=1).values
    samples = features_df['sample'].values

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # PCA
    pca = PCA(n_components=min(10, X_scaled.shape[1], X_scaled.shape[0]))
    X_pca = pca.fit_transform(X_scaled)

    # t-SNE
    if X_scaled.shape[0] > 2:
        tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, X_scaled.shape[0]-1))
        X_tsne = tsne.fit_transform(X_scaled)
    else:
        X_tsne = X_scaled[:, :2] if X_scaled.shape[1] >= 2 else np.zeros((X_scaled.shape[0], 2))

    # UMAP
    if umap and X_scaled.shape[0] > 2:
        reducer = umap.UMAP(n_components=2, random_state=42)
        X_umap = reducer.fit_transform(X_scaled)
    else:
        X_umap = X_scaled[:, :2] if X_scaled.shape[1] >= 2 else np.zeros((X_scaled.shape[0], 2))

    # Save embeddings
    embeddings_df = pd.DataFrame({
        'sample': samples,
        'PC1': X_pca[:, 0] if X_pca.shape[1] > 0 else 0,
        'PC2': X_pca[:, 1] if X_pca.shape[1] > 1 else 0,
        'UMAP1': X_umap[:, 0],
        'UMAP2': X_umap[:, 1],
        'tSNE1': X_tsne[:, 0],
        'tSNE2': X_tsne[:, 1]
    })
    embeddings_file = os.path.join(output_dir, "patterns", "sample_embeddings.csv")
    embeddings_df.to_csv(embeddings_file, index=False)
    print(f"✓ Embeddings: {embeddings_file}")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].scatter(embeddings_df['PC1'], embeddings_df['PC2'], alpha=0.6, s=50)
    axes[0].set_xlabel('PC1')
    axes[0].set_ylabel('PC2')
    axes[0].set_title('PCA')
    axes[0].grid(True, alpha=0.3)

    axes[1].scatter(embeddings_df['UMAP1'], embeddings_df['UMAP2'], alpha=0.6, s=50)
    axes[1].set_xlabel('UMAP1')
    axes[1].set_ylabel('UMAP2')
    axes[1].set_title('UMAP')
    axes[1].grid(True, alpha=0.3)

    axes[2].scatter(embeddings_df['tSNE1'], embeddings_df['tSNE2'], alpha=0.6, s=50)
    axes[2].set_xlabel('t-SNE1')
    axes[2].set_ylabel('t-SNE2')
    axes[2].set_title('t-SNE')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plot_file = os.path.join(output_dir, "plots", "dimensionality_reduction.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Dimensionality reduction plot: {plot_file}")

    # Hierarchical clustering
    if X_scaled.shape[0] > 1:
        linkage_matrix = linkage(X_scaled, method='ward')

        fig, ax = plt.subplots(figsize=(max(12, len(samples)*0.3), 6))
        dendrogram(linkage_matrix, labels=samples, ax=ax, leaf_font_size=8)
        ax.set_xlabel('Sample')
        ax.set_ylabel('Distance')
        ax.set_title('Hierarchical Clustering (SV Patterns)')
        plt.xticks(rotation=90)
        plt.tight_layout()

        dendro_file = os.path.join(output_dir, "plots", "clustering_dendrogram.png")
        plt.savefig(dendro_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Dendrogram: {dendro_file}")


def main():
    print("="*80)
    print("STEP 3: MATRIX BUILDING AND ANALYSIS")
    print("="*80)
    print()

    # Create directories
    os.makedirs(os.path.join(OUTPUT_DIR, "matrices"), exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, "hotspots"), exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, "genes"), exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, "patterns"), exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, "plots"), exist_ok=True)

    # Parse merged VCF
    print("Parsing merged VCF...")
    events_df, matrix_df, samples = parse_merged_vcf(MERGED_VCF)
    print(f"  Total SVs: {len(events_df)}")
    print(f"  Samples: {len(samples)}")
    print()

    # Save matrix
    matrix_file = os.path.join(OUTPUT_DIR, "matrices", "sample_by_sv_matrix.csv")
    matrix_df.to_csv(matrix_file)
    print(f"✓ Matrix saved: {matrix_file}")
    print(f"  Shape: {matrix_df.shape}")
    print()

    # Save events
    events_file = os.path.join(OUTPUT_DIR, "matrices", "sv_events.csv")
    events_df.to_csv(events_file, index=False)
    print(f"✓ Events saved: {events_file}")
    print()

    # Sample summary
    print("Computing sample summaries...")
    summary_df = save_sample_summary(matrix_df, events_df, OUTPUT_DIR)
    print()

    # Recurrent SVs
    print("Identifying recurrent SVs...")
    recurrent_df = identify_recurrent_svs(events_df, OUTPUT_DIR)
    print()

    # Hotspots
    print("Identifying hotspots...")
    hotspots_df = identify_hotspots(events_df, OUTPUT_DIR)
    print()

    # Gene analysis
    print("Analyzing gene-level recurrence...")
    gene_df = analyze_genes(events_df, matrix_df, GENE_ANNOTATION, GBM_DRIVERS, OUTPUT_DIR)
    print()

    # SV features
    print("Computing SV feature vectors...")
    features_df = compute_sv_features(matrix_df, events_df, OUTPUT_DIR)
    print()

    # Clustering
    print("Performing clustering and dimensionality reduction...")
    perform_clustering(features_df, OUTPUT_DIR)
    print()

    print("="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()
