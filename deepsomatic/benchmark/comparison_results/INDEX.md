# Complete Index of Benchmark Comparison Files

## 📊 Main Comparison Tables

### 1. **gene_comparison_simple.csv** ⭐ RECOMMENDED
**Simple gene-level comparison table**
- Easy to read in Excel/spreadsheet
- Shows which genes each caller detected
- Highlights common and unique genes per sample

**Columns:**
- `sample_id` - Sample identifier
- `deepsomatic_genes` - Genes detected by DeepSomatic
- `clairsto_genes` - Genes detected by ClairS-TO  
- `common_genes` - Genes found by BOTH
- `missing_in_clairsto` - DeepSomatic-only genes
- `missing_in_deepsomatic` - ClairS-TO-only genes

### 2. **quality_metrics_comparison.csv** ⭐ NEW
**Quality metrics comparison for common variants**
- Compares GQ, Depth, AD, GT, AF for the same variants
- Shows 31 common variants detected by both callers

**Key Finding:** DeepSomatic calls all variants as 1/1 (homozygous), ClairS-TO as 0/1 (heterozygous)

**Columns:**
- `sample_id`, `gene`, `variant`
- `ds_gt`, `clair_gt`, `gt_match`
- `ds_gq`, `clair_gq`, `gq_diff`
- `ds_depth`, `clair_depth`, `depth_diff`
- `ds_ad`, `clair_ad`
- `ds_af`, `clair_af`, `af_diff`

### 3. **per_sample_comparison.csv**
**Statistical summary per sample**
- Variant counts, concordance percentage
- VAF ranges, gene counts
- Pathogenic variants, COSMIC hits

### 4. **variant_comparison_detailed.csv**
**Detailed variant descriptions**
- Full variant information with amino acid changes
- Format: Gene:p.change (chr:pos ref>alt VAF=%)

## 📁 Per-Sample Detailed Tables

**Directory:** `per_sample_tables/`

**13 individual CSV files** (one per sample):
- T18-020_variants_detailed.csv
- T18-021_variants_detailed.csv
- T18-022_variants_detailed.csv
- T18-023_variants_detailed.csv
- T18-024_variants_detailed.csv
- T18-026_variants_detailed.csv
- T19-001_variants_detailed.csv
- T19-003_variants_detailed.csv
- T19-008_variants_detailed.csv
- T19-009_variants_detailed.csv
- T19-011_variants_detailed.csv
- T23-185_variants_detailed.csv
- T23-190_variants_detailed.csv

**Each file contains:**
All variants for that sample from BOTH callers with full details:
- Gene, Nucleotide_alteration, Ref, Alt
- Func (functional region)
- CLNSIG (clinical significance)
- COSMIC100 (cancer database ID)
- GQ, Depth, AD, GT, AF (quality metrics)
- Chr, Start, End
- **Variant_caller** (DeepSomatic or ClairS-TO)

**Columns (16 total):** Gene, Nucleotide_alteration, Ref, Alt, Func, CLNSIG, COSMIC100, GQ, Depth, AD, GT, AF, Chr, Start, End, Variant_caller

## 📄 Summary Reports

### 1. **QUALITY_METRICS_SUMMARY.md** ⭐ NEW
**Comprehensive quality comparison analysis**
- Genotype calling differences (1/1 vs 0/1)
- GQ, Depth, AF comparison statistics
- Interpretation and recommendations
- Clinical reporting guidance

### 2. **VISUAL_BENCHMARK_REPORT.txt**
**Visual report with ASCII charts**
- Bar charts for variant counts
- Gene frequency analysis
- Venn diagrams
- VAF distribution visualizations
- Per-sample gene comparison with visual bars
- Concordance distribution

### 3. **BENCHMARK_REPORT.md**
**Original comprehensive text report**
- Overall statistics
- Winner determination
- Key findings

### 4. **README.md**
**Documentation and usage guide**
- File descriptions
- Usage examples
- Key findings summary

## 🔍 Quick Reference

### Want to see...

**Simple gene comparison?**
→ `gene_comparison_simple.csv`

**Quality metrics for common variants?**
→ `quality_metrics_comparison.csv`
→ `QUALITY_METRICS_SUMMARY.md`

**Detailed variants for specific sample?**
→ `per_sample_tables/{sample_id}_variants_detailed.csv`

**Visual charts and graphs?**
→ `VISUAL_BENCHMARK_REPORT.txt`

**Overall statistics?**
→ `BENCHMARK_REPORT.md` or `README.md`

## 📈 Key Statistics

- **Total samples analyzed:** 13
- **DeepSomatic total variants:** 66
- **ClairS-TO total variants:** 33
- **Common variants:** 31 (detected by both)
- **DeepSomatic unique:** 35 variants
- **ClairS-TO unique:** 2 variants
- **Winner:** DeepSomatic (7-0 across all metrics)

## ⚠️ Important Findings

1. **Genotype Discrepancy:** 
   - DeepSomatic: ALL variants called as 1/1 (homozygous)
   - ClairS-TO: ALL variants called as 0/1 (heterozygous)
   - Agreement: 0% on genotype calls

2. **VAF Concordance:**
   - Average difference: 1.17%
   - High agreement on allele frequencies

3. **Depth Filtering:**
   - Average difference: 27.2 reads
   - Different filtering strategies

4. **Sensitivity:**
   - DeepSomatic detects 100% more variants
   - Better low VAF detection (5.36% vs 13.04%)

## 📊 File Summary

```
comparison_results/
├── gene_comparison_simple.csv              ⭐ Simple gene table
├── quality_metrics_comparison.csv          ⭐ Quality metrics comparison
├── QUALITY_METRICS_SUMMARY.md              ⭐ Quality analysis report
├── VISUAL_BENCHMARK_REPORT.txt             📊 Visual charts
├── BENCHMARK_REPORT.md                     📄 Original report
├── README.md                               📖 Documentation
├── INDEX.md                                📑 This file
├── per_sample_comparison.csv               📊 Statistical summary
├── variant_comparison_detailed.csv         📋 Detailed variants
├── benchmark_summary.txt                   📄 Text summary
└── per_sample_tables/                      📁 13 sample files
    ├── T18-020_variants_detailed.csv
    ├── T18-021_variants_detailed.csv
    ├── ...
    └── T23-190_variants_detailed.csv
```

---

**Last Updated:** 2025-11-13
**Total Files:** 21 (8 summary files + 13 per-sample tables)
