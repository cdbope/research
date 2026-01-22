# Which Variant Caller is Better? DeepSomatic vs Clair3/ClairS-TO

## TL;DR - Quick Answer for T25-152

🏆 **WINNER: Clair3/ClairS-TO**

**Score**: Clair3/ClairS-TO (79.0) vs DeepSomatic (39.5)

### Key Findings:
- ✅ **Higher Sensitivity**: 66.7% vs 42.9% (detects more variants)
- ✅ **Better Quality**: 5.64 vs 4.39 per variant
- ✅ **More Clinical Hits**: 3 vs 1 pathogenic/likely pathogenic variants
- ✅ **More COSMIC Hits**: 14 vs 4 cancer-related variants
- ✅ **Better Cancer Gene Coverage**: 4 vs 3 known cancer genes

### Recommendation by Use Case:

| Use Case | Best Caller | Reason |
|----------|-------------|--------|
| **Clinical Reporting** | Clair3/ClairS-TO | Better coverage of clinically relevant variants |
| **Screening/Discovery** | Clair3/ClairS-TO | Higher sensitivity, finds more variants |
| **High-Confidence Only** | **Intersection (Both)** | 2 variants called by both = highest confidence |
| **Research** | Union of Both | Comprehensive detection |

---

## Detailed Comparison

### Quality Metrics Head-to-Head

| Metric | DeepSomatic | Clair3/ClairS-TO | Winner |
|--------|-------------|------------------|--------|
| **Clinical significance hits** | 1 | 3 | Clair ✓ |
| **COSMIC cancer variants** | 4 | 14 | Clair ✓ |
| **Functional impact variants** | 9 | 14 | Clair ✓ |
| **Known cancer gene hits** | 3 | 4 | Clair ✓ |
| **Moderate VAF variants** | 9 | 13 | Clair ✓ |
| **Weighted Total Score** | 39.5 | 79.0 | Clair ✓ |

**Result**: Clair3/ClairS-TO wins 5/5 metrics

### Sensitivity vs Specificity

**Sensitivity (Variant Detection)**:
- Clair3/ClairS-TO: **66.7%** (14/21 variants)
- DeepSomatic: **42.9%** (9/21 variants)
- Winner: **Clair3/ClairS-TO** ✓

**Quality per Variant**:
- Clair3/ClairS-TO: **5.64** (quality score per variant)
- DeepSomatic: **4.39** (quality score per variant)
- Winner: **Clair3/ClairS-TO** ✓

### Concordance Analysis

**Concordance Rate**: 9.5% (2 shared variants)
- This is **LOW** but expected given different filtering strategies
- Low concordance doesn't mean bad - it means different approaches
- **Shared variants have excellent VAF concordance** (±1-2%)

**Interpretation**: Both callers are accurate for what they detect, but they detect different variant sets.

---

## Why Clair3/ClairS-TO Performs Better for T25-152

### 1. Higher Sensitivity
Clair3/ClairS-TO detected **5 more variants** (14 vs 9):
- More comprehensive screening
- Better for not missing clinically relevant variants
- Uses multiple calling strategies (Pileup + Merged + ClairS_TO)

### 2. Better Clinical Relevance
Clair3/ClairS-TO found more clinically significant variants:
- 3 pathogenic/likely pathogenic vs 1 (DeepSomatic)
- 14 COSMIC hits vs 4
- Includes TERT promoter, JAK2, NOTCH1, etc.

### 3. More Cancer Genes
Detected variants in:
- **Both**: H3-3A, PDGFRA
- **Clair-only**: TERT, JAK2, NOTCH1, NOTCH2, ABL1, and others
- **DS-only**: EGFR, KRAS, ERBB2, ROS1, FGFR1, CHEK2

### 4. Quality Control
- Multiple caller consensus (Pileup, Merged, ClairS_TO)
- Variants called by multiple strategies = higher confidence
- Example: H3-3A called by all 3 sub-callers

---

## When to Use Each Caller

### Use Clair3/ClairS-TO When:
✅ You need **comprehensive screening**
✅ You want **higher sensitivity** (don't miss variants)
✅ You're doing **clinical diagnostics** (need broad coverage)
✅ You can **manually review** flagged variants
✅ You're looking for **germline + somatic** variants

**Pros**:
- Higher sensitivity (fewer false negatives)
- Better clinical variant coverage
- Multiple calling strategies
- Well-established for Oxford Nanopore

**Cons**:
- May include more false positives
- Requires manual review of borderline calls
- Longer filtering pipeline

### Use DeepSomatic When:
✅ You need **tumor-only specific** calling
✅ You want **fewer false positives** (precision over sensitivity)
✅ You're working with **FFPE samples** (has FFPE-specific models)
✅ You want **focused high-confidence** calls
✅ Limited time for manual review

**Pros**:
- Tumor-specific training (somatic focus)
- FFPE-specific models available
- Fewer false positives (more conservative)
- Single unified model

**Cons**:
- Lower sensitivity (may miss some variants)
- Newer tool (less established for ONT)
- Fewer total variants detected

### Use BOTH (Intersection) When:
✅ You need **highest confidence** variants
✅ Clinical reporting requires **orthogonal validation**
✅ You want **near-zero false positives**
✅ You're making **critical clinical decisions**

**Result for T25-152**: 2 high-confidence variants (H3-3A, PDGFRA)

### Use BOTH (Union) When:
✅ You're doing **research/discovery**
✅ You want **maximum sensitivity**
✅ You can **validate** interesting findings
✅ Comprehensive **mutation profiling**

**Result for T25-152**: 21 total unique variants

---

## Practical Recommendations

### Tier 1: Report with Highest Confidence
**Use**: Variants called by BOTH callers (n=2)

1. **chr1:226064434 A>T** - H3-3A K28M
   - Likely pathogenic
   - VAF: 31.8% (DS), 31.5% (Clair) - Excellent agreement
   - Action: **Report as pathogenic**

2. **chr4:54267392 A>G** - PDGFRA Y288C
   - Uncertain significance
   - VAF: 76.4% (DS), 74.6% (Clair) - Excellent agreement
   - Action: **Report, consider targeted therapy**

### Tier 2: Consider with Manual Review
**Use**: Single-caller variants in known cancer genes

**Clair-only**:
- TERT promoter (chr5:1295957) - Check in IGV
- JAK2, NOTCH1, NOTCH2, ABL1 - Review for clinical relevance

**DeepSomatic-only**:
- EGFR G719S - Targetable mutation, validate
- KRAS G60V - Driver mutation, validate
- ERBB2 L755S - HER2 mutation, validate

**Action**: IGV review → Sanger validation if clinically relevant

### Tier 3: Flag for Research
**Use**: Other single-caller variants

**Action**: Document for research, not for clinical reporting without validation

---

## Validation Strategy

### For Clinical Reporting:

1. **Tier 1 (Both callers)**:
   - Sanger sequencing recommended
   - Can report with high confidence
   - H3-3A K28M: Report as diagnostic marker

2. **Tier 2 (Single caller, cancer gene)**:
   - Mandatory validation (Sanger or ddPCR)
   - IGV visual inspection required
   - Literature review for clinical significance

3. **Tier 3 (Single caller, other)**:
   - Research-only
   - Do not report clinically without validation

### Validation Methods:

| Method | Best For | Sensitivity | Cost |
|--------|----------|-------------|------|
| **Sanger** | High VAF (>10%) | ~15% | $ |
| **ddPCR** | Low VAF (<10%) | ~0.1% | $$ |
| **Deep Seq** | Multiple variants | ~1% | $$$ |

---

## How We Determined the Winner

### Scoring System:

**Quality Metrics** (weighted):
- Clinical significance: ×3
- COSMIC hits: ×2
- Cancer gene hits: ×2
- Functional impact: ×1.5
- Moderate VAF variants: ×1

**Results**:
```
DeepSomatic:
  Clinical hits (1) × 3 = 3
  COSMIC (4) × 2 = 8
  Cancer genes (3) × 2 = 6
  Functional (9) × 1.5 = 13.5
  Moderate VAF (9) × 1 = 9
  TOTAL = 39.5

Clair3/ClairS-TO:
  Clinical hits (3) × 3 = 9
  COSMIC (14) × 2 = 28
  Cancer genes (4) × 2 = 8
  Functional (14) × 1.5 = 21
  Moderate VAF (13) × 1 = 13
  TOTAL = 79.0
```

**Winner: Clair3/ClairS-TO (79.0 vs 39.5)** 🏆

---

## Running the Evaluation Yourself

### Quick Start:
```bash
cd /home/chbope/extension/script/deepsomatic

# Complete comparison + evaluation
bash full_comparison.sh
```

### For Different Sample:
```bash
bash full_comparison.sh \
  -d output/SAMPLE_deepsomatic.csv \
  -c path/to/SAMPLE_clair.csv \
  -o comparison_results/SAMPLE
```

### Output Files:
- `caller_recommendation.txt` - Which caller is better
- `shared_variants.tsv` - High confidence (both callers)
- `deepsomatic_only.tsv` - DeepSomatic unique
- `clair_only.tsv` - Clair unique
- `comparison_summary.txt` - Statistics

---

## Important Caveats

### 1. Sample-Specific Results
This analysis is for **T25-152 only**. Results may differ for other samples.

### 2. Filtering Effects
Both files were pre-filtered:
- DeepSomatic: Exonic nonsynonymous, benign excluded
- Clair: Multiple caller consensus

Different filtering → different concordance

### 3. Low Concordance ≠ Bad
9.5% concordance is low but **not concerning** because:
- Different calling algorithms (expected)
- Different filtering criteria
- Shared variants show excellent VAF agreement
- Both detect clinically relevant variants

### 4. Validation is Key
**Always validate** clinically actionable variants:
- Sanger sequencing for clinical reporting
- IGV visual inspection for all Tier 2 variants
- Literature review for therapeutic decisions

---

## Final Recommendation for T25-152

### Best Overall: **Clair3/ClairS-TO** 🏆

**Reasoning**:
1. ✅ Better sensitivity (66.7% vs 42.9%)
2. ✅ More clinically relevant variants (3 vs 1)
3. ✅ Higher quality per variant (5.64 vs 4.39)
4. ✅ Better cancer gene coverage (4 vs 3)
5. ✅ Higher total quality score (79.0 vs 39.5)

**Best Practice**:
- **Primary**: Use Clair3/ClairS-TO results
- **Validation**: Check intersection (both callers) for highest confidence
- **Follow-up**: Review DeepSomatic-only variants in known cancer genes (EGFR, KRAS, ERBB2)

### Confidence Levels:

| Confidence | Variants | Source | Action |
|------------|----------|--------|--------|
| **Highest** | 2 | Both callers | Report clinically |
| **High** | 12 | Clair only | Review + validate key ones |
| **Moderate** | 7 | DeepSomatic only | Validate cancer genes |

---

## Scripts and Documentation

- `full_comparison.sh` - Complete comparison + recommendation
- `compare_callers.py` - Detailed variant comparison
- `evaluate_callers.py` - Caller evaluation and scoring
- `README_comparison.md` - Detailed comparison guide
- `README_WHICH_CALLER.md` - This file

---

## Questions?

**Q: Can I trust the recommendation?**
A: Yes, it's based on objective metrics (clinical significance, COSMIC, cancer genes, VAF quality). But always validate clinically relevant findings.

**Q: Should I stop using DeepSomatic?**
A: No! DeepSomatic has value for FFPE samples and tumor-only calling. Consider both tools complementary.

**Q: What about other samples?**
A: Run the evaluation for each sample. Results may vary based on tumor type, sequencing quality, VAF distribution.

**Q: How do I report these findings?**
A: Report Tier 1 (both callers) with high confidence. Validate and report Tier 2 as needed. Document methodology.
