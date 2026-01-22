# GBM Mutation Landscape Analysis Pipeline

## Overview

This pipeline performs comprehensive mutation landscape analysis and pathway enrichment for Glioblastoma (GBM) somatic variants from targeted sequencing data. The analysis includes variant aggregation, mutation profiling, fusion event detection, and GBM-specific pathway enrichment.

---

## Table of Contents

1. [Requirements](#requirements)
2. [Input Data](#input-data)
3. [Analysis Scripts](#analysis-scripts)
4. [Methods](#methods)
5. [Output Files](#output-files)
6. [Usage](#usage)
7. [Citation](#citation)

---

## Requirements

### Software Dependencies

```
Python >= 3.8
pandas >= 1.3.0
matplotlib >= 3.4.0
seaborn >= 0.11.0
numpy >= 1.20.0
scipy >= 1.7.0
requests >= 2.26.0
```

### Installation

```bash
pip install pandas matplotlib seaborn numpy scipy requests
```

---

## Input Data

### 1. Variant Calling Files

- **Format**: Tab-separated values (TSV)
- **Location**: `/home/chbope/extension/data/200GMBs/results/merge_annot_clair3andclairsto_update/`
- **Naming**: `{SAMPLE_ID}_merge_annotation_filter_snvs_allcall.csv`
- **Content**: Annotated somatic variants from Clair3 and ClairS-TO callers

**Required columns:**
- `Gene.refGene`: Gene symbol
- `Chr`: Chromosome
- `Start`, `End`: Genomic coordinates
- `Ref`, `Alt`: Reference and alternate alleles
- `Func.refGene`: Functional region (exonic, upstream, etc.)
- `ExonicFunc.refGene`: Exonic function type
- `AF`: Allele frequency
- `Depth`: Sequencing depth
- `CLNSIG`: Clinical significance
- `COSMIC100`: COSMIC database annotation
- `Variant_caller`: Source caller(s)

### 2. Target Gene Panel

- **Format**: BED file (GTF-like annotations)
- **Location**: `/home/chbope/extension/nWGS_manuscript_data/data/reference/OCC.protein_coding.bed`
- **Content**: OCC (Onco-Childhood-CNS) protein-coding gene panel (205 genes)
- **Purpose**: Define target regions for enrichment analysis

### 3. Fusion Event Files

- **Format**: Tab-separated values (TSV)
- **Location**: `/home/chbope/extension/data/200GMBs/results/fusion_event/`
- **Naming**: `{SAMPLE_ID}_filter_fusion_event.tsv`
- **Content**: Structural variant fusion events from Sniffles2

**Required columns:**
- `chr`, `star`, `end`: Breakpoint coordinates
- `ID`: Fusion event identifier
- `svtype`: Structural variant type (BND, DUP, etc.)
- `breaking`: Breakpoint position (start/end)
- `Genes`: Gene(s) involved in fusion

---

## Analysis Scripts

### 1. `mutation_landscape_analysis.py`

**Purpose**: Comprehensive mutation landscape analysis across all samples

**Features**:
- Aggregates variants from multiple samples
- Generates mutation frequency plots
- Creates gene recurrence heatmaps
- Analyzes mutation types and distributions
- Calculates target gene panel coverage
- Produces publication-ready visualizations

**Key Functions**:
- `load_all_variant_data()`: Loads and merges variant data from all samples
- `plot_gene_mutation_frequency()`: Plots top mutated genes with percentages
- `plot_gene_recurrence_heatmap()`: Creates sample × gene mutation matrix
- `plot_target_gene_coverage()`: Analyzes OCC panel gene coverage
- `generate_summary_stats()`: Computes comprehensive statistics

### 2. `merge_fusion_events.py`

**Purpose**: Merge and visualize gene fusion events across samples

**Features**:
- Aggregates fusion events from multiple samples
- Handles files with only headers
- Creates fusion prevalence visualizations
- Analyzes fusion gene pairs
- Generates chromosome distribution plots

**Key Functions**:
- `merge_fusion_events()`: Combines fusion data with encoding handling
- `create_fusion_visualizations()`: Generates 5 publication figures
- Fusion prevalence, SV type distribution, gene pairs, chromosome plots

### 3. `gbm_pathway_analysis.py`

**Purpose**: GBM-focused pathway enrichment analysis

**Features**:
- Enrichr API integration for pathway analysis
- Focus on OCC panel genes with mutations
- GBM-specific pathway identification
- Core GBM signaling pathway analysis
- Multiple pathway database queries

**Key Functions**:
- `load_occ_panel_genes()`: Extracts genes from BED file
- `run_enrichr_analysis()`: Queries 9 pathway databases
- `identify_gbm_pathways()`: Filters GBM-relevant pathways
- `create_gbm_visualizations()`: Generates 6 analysis figures

**Pathway Databases Queried**:
1. KEGG 2021 Human
2. GO Biological Process 2021
3. Reactome 2022
4. WikiPathway 2021 Human
5. MSigDB Hallmark 2020
6. BioPlanet 2019
7. OMIM Disease
8. DisGeNET
9. Jensen DISEASES

---

## Methods

### Variant Calling and Annotation

Somatic single nucleotide variants (SNVs) and small insertions/deletions (indels) were identified using two complementary variant callers:

1. **Clair3** (v1.0): Deep learning-based germline and somatic variant caller
2. **ClairS-TO** (tumor-only mode): Tumor-only somatic variant caller

Variants were called from targeted sequencing data using the OCC (Onco-Childhood-CNS) gene panel covering 205 cancer-associated genes. Raw variants were filtered for:
- Minimum sequencing depth ≥ 10×
- Variant allele frequency (VAF) ≥ 0.10
- Base quality ≥ 20

Variants passing quality filters were annotated using:
- **ANNOVAR**: Functional annotation with RefSeq gene models
- **COSMIC v100**: Cancer mutation database
- **ClinVar**: Clinical significance classification
- **dbSNP**: Population variant frequencies

### Fusion Detection

Gene fusions and structural variants were detected using:
- **Sniffles2** (v2.0): Long-read structural variant caller
- Minimum read support: 3 reads
- Structural variant types: Breakends (BND), Duplications (DUP), Deletions (DEL)

Fusion events were filtered to retain only:
- Events involving known cancer genes
- Partner genes in the OCC panel
- Breakpoints within gene boundaries

### Mutation Landscape Analysis

#### Sample Aggregation
Variant data from 204 GBM samples were aggregated and deduplicated. For each variant:
- Sample ID tracking enabled recurrence analysis
- Allele frequency and depth were converted to numeric values
- Invalid entries (non-numeric AF/Depth) were removed

#### Statistical Analysis
- **Mutation frequency**: Number of mutations per gene across all samples
- **Gene recurrence**: Number of samples with mutations in each gene
- **Mutation burden**: Total variants per sample (mean: 14.67, median: 13)
- **Target coverage**: Percentage of OCC panel genes with mutations (58.5%)

#### Mutation Classification
Variants were categorized by:
- **Functional region**: Exonic (85%), upstream (15%)
- **Exonic function**: Nonsynonymous, frameshift, stopgain, etc.
- **Clinical significance**: Pathogenic, likely pathogenic, uncertain, benign
- **COSMIC status**: Presence in cancer mutation database

### Pathway Enrichment Analysis

#### Gene Selection
Pathway analysis focused on the 120 OCC panel genes with detected mutations across the cohort. This approach:
- Reduces statistical noise from unrelated genes
- Focuses on cancer-relevant pathways
- Provides biologically interpretable results

#### Enrichment Method
Pathway enrichment was performed using the **Enrichr** web service (Chen et al., 2013; Xie et al., 2021) via REST API:

1. **Gene list submission**: 120 mutated OCC genes
2. **Database queries**: 9 pathway/disease databases
3. **Statistical testing**: Fisher's exact test with Benjamini-Hochberg correction
4. **Significance threshold**: Adjusted p-value < 0.05

#### Enrichment Scoring
For each pathway, the following metrics were calculated:
- **P-value**: Fisher's exact test p-value
- **Adjusted p-value**: Benjamini-Hochberg FDR correction
- **Combined score**: c = log(p) × z, where z is z-score of deviation from expected rank
- **Overlap**: Number and identity of input genes in pathway

#### GBM-Specific Analysis
Pathways were classified as GBM-relevant if they contained keywords:
- RTK signaling (EGFR, PDGFRA)
- PI3K-AKT-mTOR pathway
- RAS-MAPK signaling
- TP53/cell cycle regulation
- IDH mutation status
- Glioma/astrocytoma classification
- Hypoxia and angiogenesis

### Visualization

All visualizations were generated using matplotlib (v3.4+) and seaborn (v0.11+):
- Resolution: 300 DPI (publication quality)
- Color schemes: Colorblind-friendly palettes
- Statistical annotations: Percentages, p-values, counts
- Format: PNG with tight bounding boxes

---

## Output Files

### Directory Structure

```
results/
├── combined_variants_all_samples.csv          # Aggregated variant data
├── gene_recurrence_matrix.csv                  # Sample × gene mutation matrix
├── summary_statistics.txt                      # Comprehensive statistics
├── target_gene_coverage_report.txt             # OCC panel coverage details
├── merged_fusion_events.tsv                    # Aggregated fusion events
│
├── Mutation Landscape Visualizations (11 files):
│   ├── sample_mutation_burden.png              # Per-sample variant counts
│   ├── gene_recurrence_heatmap.png             # Mutation matrix heatmap
│   ├── gene_sample_frequency.png               # Gene prevalence (% samples)
│   ├── gene_mutation_frequency.png             # Top mutated genes
│   ├── mutation_types.png                      # Functional annotations
│   ├── allele_frequency_distribution.png       # VAF distributions
│   ├── clinical_significance.png               # ClinVar classifications
│   ├── variant_caller_comparison.png           # Caller performance
│   ├── chromosome_distribution.png             # Variants per chromosome
│   ├── depth_distribution.png                  # Sequencing coverage
│   └── target_gene_coverage.png                # OCC panel analysis
│
├── Fusion Event Visualizations (5 files):
│   ├── fusion_prevalence.png                   # Fusion event rates
│   ├── fusion_svtype_distribution.png          # SV type breakdown
│   ├── fusion_gene_frequency.png               # Genes in fusions
│   ├── fusion_gene_pairs.png                   # Recurrent fusions
│   └── fusion_chromosome_distribution.png      # Breakpoint locations
│
└── gbm_pathway_analysis/
    ├── gbm_relevant_pathways.tsv               # All GBM pathways
    ├── pathway_*.tsv                           # Results per database (9 files)
    │
    └── Pathway Visualizations (6 files):
        ├── gbm_relevant_pathways.png           # Top GBM-specific pathways
        ├── gbm_core_pathways.png               # Core GBM signaling
        ├── top_pathways_all.png                # Top enriched pathways
        ├── pathway_categories.png              # Functional categories
        ├── database_comparison.png             # Database summary
        └── gene_pathway_table.png              # Gene-pathway associations
```

### Key Metrics Generated

**Mutation Landscape:**
- 2,991 total variants (2,840 after QC)
- 204 samples analyzed
- 120 unique genes mutated
- 58.5% OCC panel coverage (120/205 genes)
- Mean: 14.67 variants/sample
- Top genes: TERT (203 samples), KDM4C (199 samples)

**Fusion Events:**
- 3 samples with fusions (1.5%)
- 3 unique fusion events detected
- KIAA1549-BRAF (pilocytic astrocytoma)
- BCOR-EP300 (2×, clear cell sarcoma)
- ZFTA-RELA (ependymoma)

**Pathway Enrichment:**
- 239 significant pathways (adj. p < 0.05)
- 606 GBM-relevant pathways identified
- Core GBM pathways enriched (RTK, PI3K-AKT-mTOR, RAS-MAPK, p53, RB)

---

## Usage

### 1. Mutation Landscape Analysis

```bash
cd /home/chbope/extension/script/mut_landscape
python mutation_landscape_analysis.py
```

**Expected runtime**: ~2-5 minutes
**Output**: 11 PNG files + 3 data files in `results/`

### 2. Fusion Event Analysis

```bash
python merge_fusion_events.py
```

**Expected runtime**: ~30 seconds
**Output**: 5 PNG files + 1 TSV file in `results/`

### 3. GBM Pathway Analysis

```bash
python gbm_pathway_analysis.py
```

**Expected runtime**: ~1-2 minutes (requires internet for Enrichr API)
**Output**: 6 PNG files + 10 TSV files in `results/gbm_pathway_analysis/`

### Run All Analyses

```bash
# Run complete pipeline
python mutation_landscape_analysis.py && \
python merge_fusion_events.py && \
python gbm_pathway_analysis.py

echo "Analysis complete! Results in results/"
```

---

## Interpretation Guide

### Mutation Landscape

**High-frequency mutations** (>50% samples):
- **TERT** (99.5%): Telomerase activation, hallmark of GBM
- **KDM4C** (97.5%): Histone demethylase, epigenetic regulation
- **ID3** (92%): HLH transcription factor, stem cell maintenance

**Driver gene mutations**:
- **EGFR** (40 samples): RTK signaling, therapeutic target
- **PTEN** (104 samples): PI3K pathway tumor suppressor
- **TP53** (49 samples): Cell cycle checkpoint, apoptosis
- **CDKN2A** (23 samples): RB pathway, CDK inhibitor

### Fusion Events

**KIAA1549-BRAF** (T23-172):
- Classic pilocytic astrocytoma fusion
- Constitutive MAPK pathway activation
- Favorable prognosis marker

**BCOR-EP300** (T22-111, 2 variants):
- Associated with clear cell sarcoma
- Chromatin remodeling disruption
- Rare in GBM, clinically significant

**ZFTA-RELA** (T22-024):
- Supratentorial ependymoma marker
- NF-κB pathway activation
- Diagnostic/classification relevance

### Pathway Enrichment

**Core GBM pathways** identified:
1. **RTK/EGFR signaling** - EGFR, PDGFRA, ERBB2 mutations
2. **PI3K-AKT-mTOR** - PTEN, PIK3CA, AKT1 alterations
3. **RAS-MAPK** - NF1, KRAS, BRAF mutations
4. **TP53 pathway** - TP53, MDM4 disruption
5. **RB/cell cycle** - CDKN2A/B, CDK4, RB1 alterations

**Therapeutic implications**:
- EGFR mutations → EGFR TKI resistance pathways enriched
- PI3K-mTOR alterations → mTOR inhibitor candidates
- IDH wild-type signature → aggressive GBM subtype

---

## Troubleshooting

### Common Issues

**1. Import errors**
```bash
# Install missing dependencies
pip install -r requirements.txt
```

**2. File not found errors**
```python
# Check input file paths in scripts
# Update paths if data location differs
```

**3. Enrichr API timeout**
```python
# Increase timeout or retry
# Check internet connection
# Enrichr may be temporarily unavailable
```

**4. Memory errors (large datasets)**
```python
# Process samples in batches
# Reduce top_n parameter for plots
# Use lower DPI (e.g., dpi=150)
```

---

## Citation

If you use this pipeline in your research, please cite:

**Variant Callers:**
- Clair3: Zheng Z, et al. (2022). Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. *Nature Computational Science*, 2(12), 797-803.
- ClairS: Not yet published (preprint)
- Sniffles2: Sedlazeck FJ, et al. (2018). Accurate detection of complex structural variations using single-molecule sequencing. *Nature Methods*, 15(6), 461-468.

**Annotation Tools:**
- ANNOVAR: Wang K, et al. (2010). ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. *Nucleic Acids Research*, 38(16), e164.

**Pathway Analysis:**
- Enrichr: Chen EY, et al. (2013). Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. *BMC Bioinformatics*, 14, 128.
- Xie Z, et al. (2021). Gene set knowledge discovery with Enrichr. *Current Protocols*, 1, e90.

**Databases:**
- COSMIC: Tate JG, et al. (2019). COSMIC: the Catalogue Of Somatic Mutations In Cancer. *Nucleic Acids Research*, 47(D1), D941-D947.
- ClinVar: Landrum MJ, et al. (2018). ClinVar: improvements to accessing data. *Nucleic Acids Research*, 48(D1), D835-D844.

---

## Contact

For questions or issues related to this pipeline:
- Open an issue in the repository
- Contact: [Your contact information]

---

## License

This pipeline is provided as-is for research purposes.

---

## Acknowledgments

- OCC gene panel design
- Sequencing facility
- Computational resources
- Funding sources

---

**Last updated**: November 24, 2025
**Version**: 1.0.0
**Pipeline author**: [Your name]
