# Quality Metrics Comparison: DeepSomatic vs ClairS-TO

## Overview

Analysis of 31 common variants (detected by both callers) across 13 samples to compare quality metrics.

## Critical Finding: Genotype Calling Difference

### ⚠️ Major Discrepancy

**ALL 31 common variants show different genotype calls:**

- **DeepSomatic**: Calls ALL variants as `1/1` (homozygous alternative)
- **ClairS-TO**: Calls ALL variants as `0/1` (heterozygous)

**Genotype Agreement: 0/31 (0.0%)**

### What This Means

This is a **fundamental difference** in how the callers interpret variant calls:

1. **For tumor-only samples**, `1/1` calls suggest:
   - DeepSomatic interprets these as homozygous variants in tumor cells
   - Could represent clonal variants or copy-neutral LOH

2. **ClairS-TO's `0/1` calls** suggest:
   - Heterozygous variant calls
   - More conservative interpretation
   - May reflect subclonal populations or normal cell contamination

### Implication

**This is NOT a minor difference!** The genotype call affects:
- Interpretation of variant clonality
- Estimation of tumor purity
- Clinical reporting
- Downstream analysis

## Quality Metrics Comparison

For the same 31 variants, comparing quality scores:

### 1. Genotype Quality (GQ)

```
Average Difference: 6.4
Range: 0 - 25
```

- **Relatively similar** quality scores
- Most variants differ by <10 quality points
- Both callers have high confidence in their calls

### 2. Read Depth (DP)

```
Average Difference: 27.2 reads
Range: 0 - 806 reads
```

**Key Observations:**
- Large variation in depth between callers
- Suggests **different read filtering strategies**
- One extreme case: PDGFRA variant
  - DeepSomatic: 355 reads
  - ClairS-TO: 1,161 reads (3.3× more)

**Likely causes:**
- Different base quality filters
- Different mapping quality thresholds
- Different duplicate marking strategies
- Different handling of supplementary/secondary alignments

### 3. Allele Frequency (AF/VAF)

```
Average Difference: 1.17%
Range: 0.00% - 4.12%
Variants with >2% difference: 5/31 (16%)
```

**Good Agreement:**
- Despite different depths, VAF estimates are **very similar**
- Average difference only 1.17%
- 84% of variants differ by <2%

**Top 5 Largest VAF Differences:**

| Sample | Gene | Variant | DS AF | Clair AF | Diff |
|--------|------|---------|-------|----------|------|
| T23-190 | PTEN | chr10:87933148 G>A | 86.4% | 90.5% | 4.1% |
| T19-008 | TERT | chr5:1295113 G>A | 38.8% | 42.9% | 4.1% |
| T19-011 | TP53 | chr17:7675074 C>T | 73.1% | 76.0% | 2.9% |
| T19-001 | TERT | chr5:1295135 G>A | 18.2% | 21.1% | 2.9% |
| T23-190 | TERT | chr5:1295113 G>A | 29.4% | 32.3% | 2.9% |

### 4. Allelic Depth (AD)

Both callers report similar allelic depth ratios despite different total depths, supporting the similar VAF observations.

## Example Variant Comparison

**Sample T19-009, TP53 chr17:7675185 C>T:**

| Metric | DeepSomatic | ClairS-TO | Difference |
|--------|-------------|-----------|------------|
| GT | 1/1 | 0/1 | ❌ Different |
| GQ | 37 | 32 | -5 |
| Depth | 32 | 31 | -1 |
| AD | 13,19 | 12,19 | Similar |
| AF | 59.4% | 61.3% | +1.9% |

## Recommendations

### 1. For Tumor-Only Analysis

**Critical Decision:** Understand genotype interpretation

- **DeepSomatic `1/1`**: May indicate clonal tumor variants
- **ClairS-TO `0/1`**: May indicate heterozygous or subclonal variants

⚠️ **Do NOT directly compare genotype calls between callers**

### 2. For Variant Validation

**Use VAF, not genotype:**
- VAF values are highly concordant (avg 1.17% difference)
- VAF is more biologically meaningful for tumor samples
- Both callers agree on allele frequencies

### 3. For Quality Filtering

**Genotype Quality (GQ):**
- Both callers provide high GQ scores (typically 25-40)
- Similar confidence despite different genotypes
- Average difference only 6.4 points

**Read Depth:**
- Be aware of different filtering strategies
- Don't use absolute depth cutoffs
- Consider depth relative to each caller's distribution

### 4. For Clinical Reporting

**Report both:**
1. The variant call (position, ref>alt, gene)
2. The VAF (which is concordant)
3. Note that genotype interpretation differs by caller

## Key Takeaways

✅ **High VAF concordance** (1.17% avg difference)
✅ **Similar quality scores** (GQ diff: 6.4)
❌ **Complete genotype discordance** (0% agreement on GT)
⚠️ **Different depth filtering** (27.2 reads avg difference)

## Conclusion

The two callers **agree on WHAT variants are present** (31 common variants) and **their frequencies** (1.17% VAF difference), but **disagree fundamentally on HOW to genotype them** (1/1 vs 0/1).

For tumor-only variant calling:
- **Focus on VAF**, not genotype
- **Both callers are valid** but interpret genotypes differently
- **DeepSomatic** is more permissive (66 vs 33 total variants)
- **ClairS-TO** is more conservative

The choice depends on your analysis goals:
- **Maximum sensitivity**: DeepSomatic
- **High specificity**: ClairS-TO
- **Best of both**: Use intersection (31 variants with high confidence)

---

**Generated:** 2025-11-13
**Analysis:** 31 common variants across 13 tumor samples
**File:** `quality_metrics_comparison.csv`
