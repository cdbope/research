# Complete SV Meta-Analysis Pipeline - Summary

## 🎯 Full Pipeline Overview

This directory contains a **complete end-to-end pipeline** for structural variant (SV) meta-analysis in GBM cohorts.

---

## 📋 Pipeline Steps

### Step 0: Setup (One-time)
```bash
./00_setup_environment.sh  # Automated installation
conda activate svmeta_env
```

### Step 1: Prepare VCFs
```bash
./01_prepare_vcfs.sh
```
- Decompresses VCF.gz files to plain VCF
- Output: `results/prepared_vcfs/`

### Step 2: Merge SVs with SURVIVOR
```bash
./02_merge_with_survivor.sh
```
- Merges SVs across samples into cohort-level catalog
- Output: `results/merged/merged_SV.vcf.gz`

### Step 3: Build Matrix & Analyze
```bash
python 03_build_matrix_and_analyze.py
```
- Creates sample × SV matrix
- Identifies recurrent SVs and hotspots
- Gene-level SV burden analysis
- Sample clustering (PCA, UMAP, t-SNE)
- Output: `results/matrices/`, `results/hotspots/`, `results/genes/`, `results/plots/`

### Step 4: External Dataset Comparison (Optional)
```bash
python 04_external_dataset_comparison.py
```
- Compares cohort with TCGA-GBM and PCAWG
- Identifies shared vs cohort-specific SVs
- Validates against known GBM patterns
- Output: `results/external_comparison/`

### Step 5: Pathway Enrichment (Optional)
```bash
# Install g:Profiler first (optional but recommended)
pip install gprofiler-official

python 05_pathway_enrichment.py
```
- GO, KEGG, Reactome enrichment
- GBM-specific pathway analysis
- Exports for DAVID/Enrichr/GSEA
- Output: `results/pathway_enrichment/`

---

## 📊 Complete Output Structure

```
results/
├── prepared_vcfs/              # Step 1: Decompressed VCFs
│   ├── Sample001.wf_sv.vcf
│   └── ...
│
├── merged/                     # Step 2: SURVIVOR merge
│   ├── merged_SV.vcf.gz
│   ├── merged_SV.vcf.gz.tbi
│   └── vcf_list.txt
│
├── matrices/                   # Step 3: Sample × SV matrices
│   ├── sample_by_sv_matrix.csv
│   ├── sv_events.csv
│   └── sample_summary.csv
│
├── hotspots/                   # Step 3: Recurrence analysis
│   ├── recurrent_svs.csv
│   └── sv_hotspots.csv
│
├── genes/                      # Step 3: Gene-level analysis
│   └── gene_recurrence.csv
│
├── patterns/                   # Step 3: SV features
│   ├── sv_features.csv
│   └── sample_embeddings.csv
│
├── plots/                      # Step 3: Visualizations
│   ├── dimensionality_reduction.png
│   ├── clustering_dendrogram.png
│   ├── sv_type_distribution.png
│   └── sv_size_distribution.png
│
├── external_comparison/        # Step 4: Dataset comparison
│   ├── tcga_comparison.csv
│   ├── pcawg_comparison.csv
│   ├── universal_patterns.csv
│   ├── frequency_comparison_TCGA-GBM.png
│   ├── frequency_comparison_PCAWG.png
│   ├── external_comparison_summary.png
│   └── external_comparison_report.txt
│
└── pathway_enrichment/         # Step 5: Enrichment analysis
    ├── gprofiler_results.csv
    ├── gbm_pathway_overlap.csv
    ├── input_genes.csv
    ├── enrichment_overview.png
    ├── enrichment_bubble.png
    ├── gbm_pathway_overlap.png
    ├── gene_list.txt
    ├── gene_list.gmt
    ├── gene_list_ranked.rnk
    └── enrichment_report.txt
```

---

## 🔬 Key Analyses Performed

### 1. **Recurrent SV Identification**
- SVs present in multiple samples
- Frequency across cohort
- Expected: EGFR amp (~40-50%), CDKN2A del (~50-70%), PTEN del (~30-50%)

### 2. **Genomic Hotspots**
- Regions with clustered SVs
- Sliding window approach (50 kb)
- Identifies breakpoint-prone regions

### 3. **Gene-Level SV Burden**
- Which genes are most affected?
- SV type breakdown per gene
- Flags known GBM driver genes

### 4. **Sample Clustering**
- Groups samples by SV patterns
- May correlate with molecular subtypes
- Identifies outliers (chromothripsis, high burden)

### 5. **External Validation**
- Compares with TCGA-GBM (333 samples)
- Compares with PCAWG (65 samples)
- Identifies novel vs known alterations

### 6. **Pathway Enrichment**
- GO Biological Process enrichment
- KEGG pathway analysis
- Reactome pathway analysis
- GBM-specific pathway validation

---

## 📁 All Scripts & Files

| File | Purpose |
|------|---------|
| **Setup & Installation** | |
| `00_setup_environment.sh` | Automated installation (conda + SURVIVOR) |
| `environment.yml` | Conda environment definition |
| `install_gprofiler.sh` | Install g:Profiler for enrichment |
| **Pipeline Steps** | |
| `01_prepare_vcfs.sh` | Prepare VCFs for SURVIVOR |
| `02_merge_with_survivor.sh` | Merge SVs with SURVIVOR |
| `03_build_matrix_and_analyze.py` | Core analysis script |
| `04_external_dataset_comparison.py` | Compare with TCGA/PCAWG |
| `05_pathway_enrichment.py` | Pathway enrichment analysis |
| **Helper Scripts** | |
| `sv_meta_pipeline.py` | Extended Python framework |
| `download_external_datasets.sh` | External data guide |
| **Reference Data** | |
| `gbm_driver_genes.txt` | Known GBM driver genes |
| `sample_metadata_template.txt` | Sample metadata template |
| `external_datasets/tcga_gbm_sv_summary.csv` | TCGA-GBM reference (41 SVs) |
| `external_datasets/pcawg_gbm_sv_summary.csv` | PCAWG reference (35 SVs) |
| **Documentation** | |
| `README.md` | Complete documentation |
| `QUICKSTART.md` | Quick start guide |
| `PROJECT_SUMMARY.md` | Project overview |
| `INSTALLATION_GUIDE.md` | Detailed installation |
| `EXTERNAL_COMPARISON_GUIDE.md` | External comparison guide |
| `PATHWAY_ENRICHMENT_GUIDE.md` | Enrichment analysis guide |
| `READY_TO_RUN.md` | Quick reference for Step 4 |
| `COMPLETE_PIPELINE_SUMMARY.md` | This file |

---

## ⏱️ Expected Runtime

| Step | Time (8 samples) | Time (100+ samples) |
|------|------------------|---------------------|
| Step 1 | 2-5 min | 10-20 min |
| Step 2 | 5-10 min | 30-60 min |
| Step 3 | 5-15 min | 15-45 min |
| Step 4 | 2-5 min | 5-10 min |
| Step 5 | 2-5 min | 5-10 min |
| **Total** | **16-40 min** | **65-145 min** |

---

## 🎯 Expected GBM Results

### Recurrent SVs
- **EGFR amplifications**: 40-57% of samples
- **CDKN2A/B deletions**: 49-61% of samples
- **PTEN deletions**: 30-41% of samples
- **Chr7 gain**: 50-70% of samples
- **Chr10 loss**: 60-80% of samples

### Top Genes with SVs
1. EGFR (amplifications, rearrangements)
2. CDKN2A, CDKN2B (deletions)
3. PTEN (deletions)
4. TP53 (deletions, mutations)
5. NF1 (deletions, truncations)
6. PDGFRA, CDK4, MDM2 (amplifications)
7. RB1 (deletions)

### Enriched Pathways
- RTK/RAS/PI3K signaling
- Cell cycle regulation
- p53 pathway
- MAPK signaling
- DNA damage response
- Glioma pathway (KEGG)

---

## 📖 Usage Examples

### Complete Workflow
```bash
# One-time setup
./00_setup_environment.sh
conda activate svmeta_env

# Run full pipeline
./01_prepare_vcfs.sh
./02_merge_with_survivor.sh
python 03_build_matrix_and_analyze.py
python 04_external_dataset_comparison.py
python 05_pathway_enrichment.py
```

### Quick Analysis (Steps 3-5 only)
```bash
# If you already have merged VCF
conda activate svmeta_env
python 03_build_matrix_and_analyze.py
python 04_external_dataset_comparison.py
python 05_pathway_enrichment.py
```

### Rerun with Different Parameters
```bash
# Edit parameters in scripts
# Example: Lower recurrence threshold
# Edit 03_build_matrix_and_analyze.py line 48:
MIN_RECURRENCE = 2  # Changed from 3

# Rerun
python 03_build_matrix_and_analyze.py
```

---

## 🔧 Customization

### Adjust SURVIVOR Merging
Edit `02_merge_with_survivor.sh`:
```bash
MAX_DIST=1000         # Distance threshold (bp)
SIZE_THRESHOLD=0.3    # Size similarity (30%)
```

### Change Recurrence Threshold
Edit `03_build_matrix_and_analyze.py`:
```python
MIN_RECURRENCE = 3         # Min samples for "recurrent"
HOTSPOT_WINDOW = 50000     # Hotspot window size
```

### Adjust Enrichment Parameters
Edit `05_pathway_enrichment.py`:
```python
TOP_N_GENES = 100          # Number of genes to analyze
MIN_GENE_RECURRENCE = 3    # Min samples per gene
PVALUE_THRESHOLD = 0.05    # Significance threshold
```

---

## 📚 Key References

### Tools
1. **SURVIVOR**: Jeffares et al. (2017) Bioinformatics
2. **g:Profiler**: Raudvere et al. (2019) Nucleic Acids Research
3. **bcftools**: Danecek et al. (2021) GigaScience

### GBM Genomics
1. **TCGA-GBM**: Brennan et al. (2013) Cell
2. **TCGA Nature 2008**: Comprehensive genomic characterization
3. **PCAWG**: Li et al. (2020) Nature Communications
4. **GBM Subtypes**: Verhaak et al. (2010) Cancer Cell

---

## 💡 Tips & Best Practices

### Quality Control
1. **Check merged VCF**: Verify reasonable SV count (100-500 per sample)
2. **Inspect recurrent SVs**: Should see known drivers (EGFR, CDKN2A, PTEN)
3. **Review clustering**: Look for outliers or unexpected patterns
4. **Validate with external**: 40-60% SVs should match TCGA/PCAWG

### Interpretation
1. **Recurrent SVs**: Focus on those in cancer genes
2. **Cohort-specific**: May be novel or technical - investigate further
3. **Pathway enrichment**: Expected pathways validate your analysis
4. **Sample clustering**: May correlate with clinical/molecular features

### Publication
1. Use visualizations from `plots/` and `external_comparison/`
2. Cite tools and reference datasets properly
3. Report both raw counts and frequencies
4. Mention SURVIVOR parameters used
5. Include pathway enrichment results

---

## 🆘 Troubleshooting

### Common Issues

**"VCF files not found"**
- Check INPUT_VCF_DIR path in `01_prepare_vcfs.sh`

**"SURVIVOR not found"**
- Run `00_setup_environment.sh` to install
- Or compile manually (see INSTALLATION_GUIDE.md)

**"No recurrent SVs found"**
- Lower MIN_RECURRENCE in `03_build_matrix_and_analyze.py`
- Check if merged VCF has reasonable SV counts

**"g:Profiler error"**
- Install: `pip install gprofiler-official`
- Or skip - script works without it

**"Very few external matches"**
- Try relaxing overlap threshold (50% → 30%)
- Different SV callers have different breakpoint precision

---

## ✅ Validation Checklist

- [ ] Environment installed (`conda activate svmeta_env` works)
- [ ] SURVIVOR compiled (`tools/SURVIVOR/Debug/SURVIVOR` exists)
- [ ] VCFs prepared (check `results/prepared_vcfs/`)
- [ ] Merge successful (check `results/merged/merged_SV.vcf.gz`)
- [ ] Matrix created (check `results/matrices/sample_by_sv_matrix.csv`)
- [ ] Recurrent SVs found (check `results/hotspots/recurrent_svs.csv`)
- [ ] Known drivers detected (EGFR, CDKN2A, PTEN in top genes)
- [ ] Visualizations generated (check `results/plots/`)
- [ ] External comparison complete (if running Step 4)
- [ ] Pathway enrichment complete (if running Step 5)

---

## 🎉 You're All Set!

This pipeline is **production-ready** and has everything you need for comprehensive SV meta-analysis of GBM cohorts.

For questions or issues:
1. Check the relevant `*_GUIDE.md` file
2. Review the scripts - they're well-commented
3. See `README.md` for detailed documentation

**Happy analyzing!** 🧬🔬
