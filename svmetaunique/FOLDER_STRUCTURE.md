# Folder Structure: svmetaunique vs svmeta

## Overview

The **svmetaunique** folder contains the conservative counting pipeline that is completely separate from the main **svmeta** pipeline.

---

## Folder Separation

### svmeta/ (Original Pipeline - READ ONLY)
```
svmeta/
├── 01_prepare_vcfs.sh                    # Step 01: Filter VCFs
├── 02_merge_with_survivor.sh             # Step 02: Merge with SURVIVOR
├── 03_build_matrix_and_analyze.py        # Step 03: Standard counting
├── 04_external_dataset_comparison_v2.py  # Step 04: Standard comparison
├── external_datasets/
│   ├── refseq_genes_hg38.bed            # Gene annotations (READ ONLY)
│   ├── tcga_gbm_sv_summary_hg38.csv     # TCGA reference (READ ONLY)
│   └── pcawg_gbm_sv_summary_hg38.csv    # PCAWG reference (READ ONLY)
├── gbm_driver_genes.txt                  # Driver gene list (READ ONLY)
└── results/
    ├── prepared_vcfs/filter_vcf/         # Filtered VCFs (READ ONLY)
    └── merged/
        └── merged_SV.vcf.gz              # Merged VCF (READ ONLY - used as input)
```

**Note**: svmetaunique scripts ONLY READ from svmeta/, never write to it.

---

### svmetaunique/ (Conservative Pipeline - ALL OUTPUTS HERE)
```
svmetaunique/
├── 03_build_matrix_conservative.py       # Step 03: Conservative counting
├── 04_external_dataset_comparison_conservative.py  # Step 04: Conservative comparison
├── 05_conservative_gene_counting.py      # Analysis: Compare approaches
├── run_conservative_pipeline.sh          # Master runner script
├── README.md                              # Documentation
├── FOLDER_STRUCTURE.md                    # This file
├── external_datasets/                     # Created if needed
│   └── refseq_genes_hg38.bed             # Copy of gene annotations (if created)
└── results/                               # ALL outputs go here
    ├── matrices/
    │   └── gene_sample_matrix_conservative.csv
    ├── genes/
    │   └── gene_frequencies_conservative.csv
    ├── external_comparison/
    │   ├── all_genes_comparison_conservative.csv
    │   ├── high_confidence_validated_genes_conservative.csv
    │   ├── top10_genes_conservative_vs_TCGA.csv
    │   ├── top10_genes_conservative_vs_PCAWG.csv
    │   └── COMPARISON_REPORT_CONSERVATIVE.md
    └── conservative_vs_standard_comparison.csv
```

---

## Input/Output Summary

### Inputs (Read-Only from svmeta/)
| File | Purpose | Location |
|------|---------|----------|
| **Merged VCF** | SV calls across 200 samples | `svmeta/results/merged/merged_SV.vcf.gz` |
| **Gene annotations** | RefSeq genes hg38 | `svmeta/external_datasets/refseq_genes_hg38.bed` |
| **TCGA reference** | TCGA-GBM SV frequencies | `svmeta/external_datasets/tcga_gbm_sv_summary_hg38.csv` |
| **PCAWG reference** | PCAWG-GBM SV frequencies | `svmeta/external_datasets/pcawg_gbm_sv_summary_hg38.csv` |
| **Driver genes** | Known GBM driver genes | `svmeta/gbm_driver_genes.txt` |

### Outputs (All Written to svmetaunique/)
| File | Description | Location |
|------|-------------|----------|
| **Gene matrix** | Binary gene×sample matrix | `svmetaunique/results/matrices/` |
| **Gene frequencies** | Conservative frequencies | `svmetaunique/results/genes/` |
| **External comparison** | TCGA/PCAWG comparison | `svmetaunique/results/external_comparison/` |
| **Analysis reports** | Markdown reports | `svmetaunique/results/` |

---

## File Access Guarantees

### ✅ Safe Operations (Read-Only)
- ✅ Reading merged VCF from svmeta/
- ✅ Reading gene annotations from svmeta/
- ✅ Reading external datasets from svmeta/
- ✅ Reading driver gene list from svmeta/

### ❌ Forbidden Operations (Never Performed)
- ❌ Writing to svmeta/results/
- ❌ Modifying files in svmeta/
- ❌ Creating new files in svmeta/
- ❌ Overwriting existing svmeta outputs

### ✅ All Writes Go To
- ✅ svmetaunique/results/
- ✅ svmetaunique/external_datasets/ (if created)

---

## Running the Pipeline

### Prerequisites
```bash
# Complete steps 01-02 in svmeta (if not done)
cd /home/chbope/extension/script/svmeta
bash 01_prepare_vcfs.sh
bash 02_merge_with_survivor.sh
```

### Run Conservative Pipeline
```bash
# Option 1: Run complete pipeline
cd /home/chbope/extension/script/svmetaunique
bash run_conservative_pipeline.sh

# Option 2: Run steps individually
python3 03_build_matrix_conservative.py
python3 04_external_dataset_comparison_conservative.py

# Option 3: Compare conservative vs standard
python3 05_conservative_gene_counting.py
```

---

## Verification

To verify no files were written to svmeta/:

```bash
# Check for files modified in svmeta/results/ (should be empty)
find /home/chbope/extension/script/svmeta/results -type f -mmin -60

# Check all outputs are in svmetaunique/
ls -lh /home/chbope/extension/script/svmetaunique/results/
```

---

## Key Differences

| Aspect | svmeta/ | svmetaunique/ |
|--------|---------|---------------|
| **Counting method** | Standard (count all SVs) | Conservative (max 1 per gene) |
| **Frequencies** | Can exceed 100% | Max 100% per gene |
| **Methodology** | Novel long-read insight | Literature-comparable |
| **Outputs** | Standard results | Conservative results |
| **File writes** | To svmeta/results/ | To svmetaunique/results/ |

---

## Important Notes

1. **No Interference**: svmetaunique does NOT modify any files in svmeta/
2. **Read-Only Access**: svmetaunique only READS from svmeta/ (input data)
3. **Separate Results**: All svmetaunique outputs go to svmetaunique/results/
4. **Can Run Both**: Both pipelines can coexist and run independently
5. **Comparison Possible**: You can compare results from both approaches

---

## Contact

For questions about:
- **svmeta pipeline**: See svmeta/README.md
- **svmetaunique pipeline**: See svmetaunique/README.md
- **Comparison**: Run 05_conservative_gene_counting.py
