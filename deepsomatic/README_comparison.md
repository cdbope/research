# Variant Caller Comparison Tool

Compare variant calls between DeepSomatic and Clair3/ClairS-TO callers.

## Key Findings for T25-152

The comparison reveals **low concordance (9.5%)** between the two callers:

| Metric | Count | Percentage |
|--------|-------|------------|
| **DeepSomatic variants** | 9 | - |
| **Clair3/ClairS-TO variants** | 14 | - |
| **Shared variants** | 2 | 9.5% |
| **DeepSomatic only** | 7 | 33.3% |
| **Clair3/ClairS-TO only** | 12 | 57.1% |

### Shared Variants (High Confidence)

Only **2 variants** were called by both callers:
1. **chr1:226064434 A>T** - H3-3A (K28M mutation)
   - VAF: DeepSomatic 31.8%, Clair 31.5%
   - Clinical: Likely pathogenic
   - COSMIC: COSV64731746

2. **chr4:54267392 A>G** - PDGFRA (Y288C mutation)
   - VAF: DeepSomatic 76.4%, Clair 74.6%
   - Clinical: Uncertain significance
   - COSMIC: COSV57264810

Both shared variants show **excellent VAF concordance** (100% within ±10%).

## Quick Start

### Option 1: Use Default Paths (T25-152)

```bash
cd /home/chbope/extension/script/deepsomatic
python3 compare_callers.py
```

### Option 2: Use Shell Wrapper

```bash
bash compare_callers.sh
```

### Option 3: Custom Files

```bash
python3 compare_callers.py \
    path/to/deepsomatic.csv \
    path/to/clair.csv \
    output_directory
```

## Output Files

The tool generates 4 files in `comparison_results/`:

1. **shared_variants.tsv** - Variants called by both callers
   - Includes VAF from both callers
   - VAF difference calculation
   - Full annotation

2. **deepsomatic_only.tsv** - Variants unique to DeepSomatic (7 variants)
   - Possible reasons: Higher sensitivity, different calling algorithms

3. **clair_only.tsv** - Variants unique to Clair3/ClairS-TO (12 variants)
   - Includes variants from multiple callers (Pileup, Merged, ClairS_TO)

4. **comparison_summary.txt** - Quick stats summary

## Analysis Metrics

### 1. Concordance Rate
Percentage of shared variants out of total unique variants.

**For T25-152**: 9.52%
- This is relatively low, suggesting different sensitivity/specificity

### 2. VAF Comparison
For shared variants, compares Variant Allele Frequency (VAF):
- Mean difference: 1.06%
- 100% within ±10% threshold

### 3. Gene-Level Overlap
- Shared genes with variants: 2 (H3-3A, PDGFRA)
- DeepSomatic unique genes: 5
- Clair unique genes: 8

### 4. Functional Impact
Both callers primarily detected exonic variants:
- DeepSomatic: 9/9 exonic (100%)
- Clair: 13/14 exonic (93%)

### 5. Clinical Significance
- Pathogenic: 0 (DS), 1 (Clair)
- Likely pathogenic: 1 (DS), 1 (Clair) - **Shared: H3-3A**
- Uncertain: 2 (DS), 1 (Clair)

## Interpreting Results

### Why Low Concordance?

Several factors may contribute:

1. **Different Filtering Criteria**
   - DeepSomatic: PASS filter + custom exonic/nonsynonymous filter
   - Clair3/ClairS-TO: Multiple caller merge (Pileup, Merged, ClairS_TO)

2. **Different Sensitivity/Specificity Trade-offs**
   - DeepSomatic (tumor-only): May be more conservative
   - Clair3 (germline caller): Higher sensitivity, may call more variants

3. **Model Training**
   - DeepSomatic: Trained on somatic variants
   - Clair3: Trained on germline variants

4. **Pre-filtering Steps**
   - Your DeepSomatic: Filtered for exonic nonsynonymous (excludes benign)
   - Your Clair: Includes all passing variants

### High Confidence Variants

Variants called by **both callers** are highest confidence:
- chr1:226064434 - H3-3A K28M (known cancer driver)
- chr4:54267392 - PDGFRA Y288C

### Investigating Discordant Variants

#### DeepSomatic-Only Variants (7)
Check:
- Are they low VAF? (May be below Clair's threshold)
- Are they in difficult regions? (homopolymers, low complexity)
- Do they have low depth?

View with:
```bash
less comparison_results/deepsomatic_only.tsv
```

#### Clair-Only Variants (12)
Check:
- Are they germline variants? (DeepSomatic is tumor-only focused)
- Do they pass your custom filters?
- Which caller detected them? (check Variant_caller column)

View with:
```bash
less comparison_results/clair_only.tsv
```

## Recommendations

### For High-Confidence Variant List
Use **intersection** (shared variants only):
```bash
# These 2 variants have highest confidence
cat comparison_results/shared_variants.tsv
```

### For Comprehensive Screening
Use **union** (all variants from both):
```bash
# Combine all unique variants
# Review DeepSomatic-only and Clair-only manually
```

### For Clinical Reporting
1. **High confidence**: Report shared variants
2. **Manual review**: Check caller-specific variants for:
   - Known cancer genes
   - High clinical significance
   - Literature support

## Advanced Usage

### Compare Different Samples

```bash
# For sample T001
python3 compare_callers.py \
    output/T001_annotateandfilter_deep_somatic.csv \
    /path/to/T001_clair.csv \
    comparison_results/T001
```

### Batch Comparison

Create a script to compare multiple samples:

```bash
#!/bin/bash
for sample in T25-152 T001 T002; do
    echo "Comparing $sample..."
    python3 compare_callers.py \
        output/${sample}_annotateandfilter_deep_somatic.csv \
        /path/to/${sample}_clair.csv \
        comparison_results/${sample}
done
```

### Export for Visualization

The TSV files can be imported into:
- Excel/LibreOffice for manual review
- R/Python for custom plots
- IGV for visual inspection

```R
# Example R code
library(tidyverse)
library(VennDiagram)

shared <- read_tsv("comparison_results/shared_variants.tsv")
ds_only <- read_tsv("comparison_results/deepsomatic_only.tsv")
clair_only <- read_tsv("comparison_results/clair_only.tsv")

# Create Venn diagram
venn.plot <- draw.pairwise.venn(
    area1 = 9,
    area2 = 14,
    cross.area = 2,
    category = c("DeepSomatic", "Clair3/ClairS-TO")
)
```

## Dependencies

- Python 3
- pandas library

Install pandas:
```bash
pip3 install pandas
```

## File Format Requirements

### DeepSomatic File
Expected columns:
- Chr, Start, End, Ref, Alt
- Gene.refGene, ExonicFunc.refGene, AAChange.refGene
- CLNSIG, COSMIC100
- Otherinfo10 (contains VAF)

### Clair File
Expected columns:
- Chr, Start, End, Ref, Alt
- Gene.refGene, ExonicFunc.refGene, AAChange.refGene
- CLNSIG, COSMIC100
- AF (Variant Allele Frequency)
- Variant_caller

## Troubleshooting

### "File not found" Error
```bash
# Check if files exist
ls -lh /home/chbope/extension/script/deepsomatic/output/T25-152_annotateandfilter_deep_somatic.csv
ls -lh /home/chbope/extension/data/t25-152/merge_annot_clair3andclairsto/T25-152_merge_annotation_filter_snvs_allcall.csv
```

### "pandas not found" Error
```bash
pip3 install pandas
```

### Empty Output Files
- Check if input files have variants
- Verify file formats match expected structure

## Understanding the Numbers

### What is Normal Concordance?

Concordance between callers varies widely:
- **Germline variants**: 95-99% (high agreement)
- **Somatic variants**: 30-70% (moderate agreement)
- **Low VAF somatic**: 10-40% (low agreement)

Your **9.5% concordance** is low, but may be explained by:
1. Different filtering strategies (you filtered more stringently)
2. Tumor-only vs germline-trained models
3. Low VAF variants (harder to call consistently)

### VAF Concordance: 100% ✓

For the 2 shared variants, VAF agreement is excellent (within 1-2%), indicating:
- Both callers accurately estimate VAF
- These are real variants, not artifacts

## Further Investigation

### Visual Inspection in IGV

For discordant variants:
```bash
# Extract coordinates
awk -F'\t' '{print $2":"$3"-"$4}' comparison_results/deepsomatic_only.tsv

# Load in IGV with both BAM files
# Check:
# - Read support
# - Strand bias
# - Mapping quality
# - Base quality
```

### Database Lookup

Check variants in:
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [COSMIC](https://cancer.sanger.ac.uk/cosmic)
- [gnomAD](https://gnomad.broadinstitute.org/) (for germline filtering)

## References

- DeepSomatic: https://github.com/google/deepvariant/blob/r1.9/docs/deepsomatic-details.md
- Clair3: https://github.com/HKU-BAL/Clair3
- ClairS-TO: https://github.com/HKU-BAL/ClairS-TO
