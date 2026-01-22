# Guide to VAF Analysis for Structural Variants in GBM

## The Problem with Simple Mean/Median VAF

### Why Simple Averages Are Misleading

When you calculate mean or median VAF across ALL structural variants, you're including:

1. **Low-quality variants**
   - IMPRECISE calls (fuzzy breakpoints)
   - Low read support (SUPPORT < 10)
   - Failed quality filters (FILTER != PASS)
   - These have unreliable VAF estimates

2. **Artifacts and noise**
   - Chromosome-specific sequencing biases
   - Repetitive region artifacts
   - Alignment errors
   - These contribute false low-VAF signals

3. **Passenger mutations**
   - Neutral variants with low VAF
   - Late subclonal events of no biological significance
   - These dilute the driver signal

4. **Result**: Mean VAF that doesn't represent true tumor biology

### Example from Your Data

```
Sample Sample1:
- ALL variants (36,070): Mean VAF = 0.632
- HIGH-QUALITY variants (25,908): Mean VAF = 0.701
- CLONAL variants (31,569): Mean VAF = 0.691

Difference: 0.069 (6.9 percentage points!)
```

**Interpretation**: The "all variants" mean is pulled down by ~10,000 low-quality variants with inflated low-VAF estimates.

## Solution: Version 2's Multi-Tier VAF Approach

### Tier 1: High-Quality VAF (RECOMMENDED)

**Filtering Criteria:**
```
✓ FILTER = "PASS"
✓ PRECISE breakpoints (not IMPRECISE)
✓ SUPPORT ≥ 10 reads
✓ VAF ≥ 0.1 (removes noise)
```

**What this gives you:**
- Reliable, well-supported variants only
- Accurate VAF estimates
- Suitable for clustering and statistical analysis
- Best for publication-quality results

**Metrics:**
- `Mean_VAF_HQ`: Average VAF of high-quality variants
- `Median_VAF_HQ`: Robust central tendency
- `HQ_Variant_Count`: How many passed filters

### Tier 2: Clonal VAF

**Definition:**
```
Clonal variant: VAF ≥ 0.3
Subclonal variant: 0.1 ≤ VAF < 0.3
```

**Why VAF ≥ 0.3?**
- Theoretical: Diploid cell with heterozygous variant in 60% of cells = VAF 0.3
- Practical: Distinguishes early (clonal) from late (subclonal) events
- For GBM: Identifies likely driver mutations

**Metrics:**
- `Mean_VAF_Clonal`: Average VAF of clonal variants (VAF ≥ 0.3)
- `Clonal_Fraction`: Proportion of variants that are clonal
- `Clonal_Count`: Number of clonal variants

**Interpretation:**
- High Clonal_Fraction (>0.8): Homogeneous, clonal tumor
- Low Clonal_Fraction (<0.5): Heterogeneous, subclonal diversity
- For GBM: Most should have high clonal fraction (aggressive monoclonal expansion)

### Tier 3: Weighted VAF

**Concept:**
Instead of equal weighting, use VAF itself as the weight:

```
Weighted_VAF = Σ(VAF_i × VAF_i) / Σ(VAF_i)
```

**Effect:**
- High-VAF variants contribute more to the average
- Low-VAF variants contribute less
- Emphasizes biologically significant events

**When to use:**
- Focus on driver events
- De-emphasize subclonal noise
- Homogeneous tumors

### Tier 4: Top Percentile VAF

**Definition:**
```
Top10_VAF = mean of top 10% highest VAF variants
Top25_VAF = mean of top 25% highest VAF variants
```

**What this tells you:**
- Presence of homozygous deletions (VAF → 1.0)
- High-confidence clonal events
- Least affected by artifacts

**Interpretation:**
- Top10_VAF ≈ 1.0: Homozygous deletions/duplications present
- Top10_VAF < 0.7: Mostly heterozygous events
- For GBM: Common to see Top10_VAF = 1.0 (homozygous CDKN2A deletion, etc.)

## Choosing the Right VAF Metric for Different Analyses

### For Clustering GBM Samples

**RECOMMENDED: High-Quality VAF**

```python
Features: ['Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS',
          'Mean_VAF_HQ', 'Median_VAF_HQ', 'Clonal_Fraction']
```

**Why:**
- Robust to technical variation
- Captures both SV burden AND clonality
- Balances sensitivity and specificity
- Best for identifying genomic subtypes

**Alternative: Clonal VAF** (for evolution-focused analysis)

```python
Features: ['Purity', 'DEL', 'DUP', 'INV', 'BND', 'INS',
          'Mean_VAF_Clonal', 'Median_VAF_Clonal', 'Clonal_Fraction']
```

### For Tumor Purity Correlation

**RECOMMENDED: Mean_VAF_HQ or Mean_VAF_Clonal**

Higher purity should correlate with higher clonal VAF (assuming monoclonal tumor).

### For Tumor Evolution Analysis

**RECOMMENDED: Clonal_Fraction**

- Track proportion of clonal vs subclonal events
- Infer timing of events (early = clonal, late = subclonal)
- Assess intratumoral heterogeneity

### For Driver Event Identification

**RECOMMENDED: Top10_VAF or Mean_VAF_Clonal**

- Focus on highest-VAF variants (most likely drivers)
- Filter out passenger mutations
- Identify homozygous deletions

## Practical Workflow for Your GBM Data

### Step 1: Run Both Versions

```bash
# Version 1 (basic, all VAF)
python sv_analysis.py

# Version 2 (enhanced, significant VAF)
python sv_analysis_v2.py
```

### Step 2: Compare Results

Look at `comprehensive_vaf_analysis.png` (Version 2):
- Panel 1: Are "All VAF" and "HQ VAF" similar? If very different, many low-quality variants.
- Panel 2: Does Clonal_Fraction correlate with Purity? Should correlate for monoclonal GBM.
- Panel 4: Clonal vs Subclonal counts - GBM should be mostly clonal.

### Step 3: Choose Primary Metric

**For publication**: Use **Mean_VAF_HQ** and **Clonal_Fraction**

Report like this:
```
"We analyzed high-quality structural variants (PASS, PRECISE, Support ≥10)
and calculated mean VAF for each sample. Clonal fraction (proportion of
variants with VAF ≥0.3) was used to assess tumor clonality."
```

### Step 4: Quality Control

Check each sample:

```
Good quality GBM sample:
✓ HQ_Variant_Count > 20,000
✓ Clonal_Fraction > 0.8
✓ Mean_VAF_HQ > 0.65
✓ Top10_VAF ≈ 1.0

Potential issues:
✗ Clonal_Fraction < 0.5 (very heterogeneous or contamination)
✗ Mean_VAF_HQ < 0.4 (low purity or quality issues)
✗ Large gap between Mean_VAF_All and Mean_VAF_HQ (many artifacts)
```

### Step 5: Clustering

Use **high_quality strategy** (default in Version 2):

```
Results: results_v2/clustering_dendrogram_high_quality.png
```

Compare with **clonal strategy** to see if clusters are stable.

## Real Example Interpretation

### Sample: Sample7

```
Purity: 0.52 (relatively low)
Mean_VAF_All: 0.649
Mean_VAF_HQ: 0.725
Mean_VAF_Clonal: 0.703
Clonal_Fraction: 0.889
Clonal_Count: 30,402
Subclonal_Count: 3,798
Weighted_VAF: 0.803
Top10_VAF: 1.0
```

**Interpretation:**
1. **Purity**: Lower (0.52) but still usable
2. **VAF Gap**: 0.725 - 0.649 = 0.076 (7.6% difference)
   - Indicates some low-quality variants affecting "All VAF"
3. **Clonal Fraction**: 0.889 (88.9%)
   - Despite lower purity, tumor is highly clonal
   - Suggests aggressive monoclonal expansion (typical GBM)
4. **Top10_VAF**: 1.0
   - Homozygous deletions present
   - Likely CDKN2A or other GBM drivers
5. **Conclusion**: Good quality sample, highly clonal GBM

### Clustering Assignment

- **High-Quality Strategy**: Cluster 3 (alone)
- **Clonal Strategy**: Cluster 1 (with Sample2, Sample3)

**Why different?**
- High-quality considers overall VAF distribution
- Clonal focuses on clonal architecture
- Sample7 has unique high-quality VAF profile but similar clonal structure to others

**Which is correct?** Both are valid - depends on biological question:
- Genomic subtyping: Use high-quality
- Evolutionary analysis: Use clonal

## Common Questions

### Q: Should I always use high-quality VAF?

**A:** For GBM analysis, YES. It's the most robust metric. Only use "All VAF" for:
- Initial exploration
- Comparison with published studies using all variants
- Sensitivity analysis

### Q: What if my clonal fraction is low (<0.5)?

**A:** Could indicate:
1. Tumor heterogeneity (multiple subclones)
2. Lower purity than estimated
3. Technical issues (contamination, quality)

For GBM, expect high clonal fraction (>0.8). Investigate low values.

### Q: Can I adjust the clonality threshold (0.3)?

**A:** Yes, but justify it:
- VAF ≥ 0.3 is standard for diploid genome
- VAF ≥ 0.5 for more stringent "clonal" definition
- VAF ≥ 0.2 for permissive definition (more heterogeneous tumors)

For GBM: 0.3 is appropriate.

### Q: Which clustering strategy should I use?

**A:** For GBM genomic subtyping:
1. **Primary**: high_quality strategy
2. **Validation**: Compare with clonal and weighted
3. **Report**: If all three agree, clusters are robust

### Q: What about chromosome-specific artifacts?

**A:** Version 2 filters many artifacts through:
1. PASS filter (removes flagged variants)
2. PRECISE filter (removes imprecise calls common in repeats)
3. Support ≥ 10 (removes low-confidence calls)

Remaining artifacts should be rare in high-quality set.

## Summary Recommendations

### ✅ DO:
- Use **Version 2** for analysis
- Report **Mean_VAF_HQ** as primary VAF metric
- Include **Clonal_Fraction** for clonality assessment
- Use **high_quality clustering strategy**
- Check **comprehensive_vaf_analysis.png** for QC

### ❌ DON'T:
- Report only "All VAF" (mean of all variants)
- Ignore quality filtering
- Use single clustering strategy without validation
- Overlook samples with low clonal fraction

### 📊 FOR PUBLICATION:
**Methods:**
"Structural variants were called using Sniffles2. High-quality variants (FILTER=PASS, PRECISE, read support ≥10) were used for downstream analysis. Variant allele frequency (VAF) was calculated for high-quality variants only. Clonal fraction was defined as the proportion of variants with VAF ≥0.3. Hierarchical clustering was performed using purity, SV counts, mean high-quality VAF, median high-quality VAF, and clonal fraction."

**Results:**
"Samples showed high clonal fractions (mean: 0.90, range: 0.89-0.93), consistent with monoclonal GBM tumors. Clustering based on high-quality VAF metrics identified three distinct genomic subtypes."

## File Reference

- **Version 1**: `/home/chbope/extension/script/svanalysis/sv_analysis.py`
- **Version 2**: `/home/chbope/extension/script/svanalysis/sv_analysis_v2.py`
- **Documentation**:
  - `README.md` (Version 1)
  - `README_V2.md` (Version 2)
  - `VAF_ANALYSIS_GUIDE.md` (this file)
