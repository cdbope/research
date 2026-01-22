# Enhanced Pathway Enrichment Analysis - New Features

## Overview

I've merged advanced features from `/home/chbope/extension/script/mut_landscape/gbm_pathway_analysis.py` into the SV pathway enrichment script [05_pathway_enrichment.py](05_pathway_enrichment.py).

## New Features Added

### 1. **Enrichr API Integration** ✨

- **Function**: `run_enrichr_analysis(gene_list)`
- **What it does**: Queries multiple pathway databases via Enrichr API
- **Databases included**:
  - KEGG_2021_Human
  - GO_Biological_Process_2021
  - Reactome_2022
  - WikiPathway_2021_Human
  - MSigDB_Hallmark_2020
- **Benefit**: Provides alternative/complementary enrichment results to g:Profiler

### 2. **Expanded GBM Pathway Definitions** 📚

#### Core GBM Signaling Pathways (`CORE_GBM_PATHWAYS`):
- **RTK/EGFR Signaling**: 9 genes (EGFR, PDGFRA, MET, FGFR1/2/3, ERBB2/3, etc.)
- **PI3K-AKT-mTOR**: 10 genes (PIK3CA, PIK3R1/R2, AKT1/2/3, MTOR, PTEN, TSC1/2)
- **RAS-MAPK**: 10 genes (KRAS, NRAS, HRAS, BRAF, RAF1, MAP2K1/2, MAPK1/3, NF1)
- **TP53 Pathway**: 8 genes (TP53, MDM2/4, CDKN2A, ATM/R, CHEK1/2)
- **RB/Cell Cycle**: 10 genes (RB1, CDKN2A/B, CDKN1A/B, CDK4/6, CCND1/2/3)
- **Angiogenesis**: 8 genes (VEGFA, VEGFR1/2, HIF1A, EPAS1, VHL, ANGPT1/2)

#### GBM Pathway Keywords (`GBM_PATHWAY_KEYWORDS`):
Comprehensive list covering **19 canonical GBM pathways**:
1. RTK–RAS–MAPK signaling
2. PI3K–AKT–mTOR pathway
3. p53 signaling
4. RB / G1-S cell-cycle axis
5. Telomere maintenance (TERT/ATRX)
6. Chromatin remodeling (SWI/SNF, epigenetic)
7. IDH pathway (oncometabolite)
8. MGMT & TMZ resistance
9. DNA repair (HR, NHEJ, MMR, BER)
10. NF-κB signaling
11. JAK–STAT3 signaling
12. Notch / Wnt / SHH stemness pathways
13. Apoptosis & necroptosis
14. Autophagy regulation
15. Hypoxia / HIF signaling
16. Angiogenesis / VEGF
17. Immune evasion (PD-L1, CD47, TAM/M2)
18. Invasion / EMT / mesenchymal transition
19. Metabolic reprogramming (glycolysis, glutamine)

### 3. **New Analysis Functions** 🔬

#### `identify_gbm_pathways(enrichr_results)`
- Filters Enrichr results for GBM-relevant pathways
- Uses the 19 canonical GBM pathway keywords
- Returns ranked dataframe of GBM-specific enrichments

#### `analyze_core_gbm_pathways(gene_list)`
- Analyzes overlap with the 6 core GBM signaling pathways
- Calculates enrichment ratios
- Shows which key GBM pathways are affected by SVs

### 4. **Improved Bubble Plot** 🎨

The bubble plot (enrichment_bubble.png) now includes:
- **Term name labels** on each bubble
- **Gene labels** showing first 5 genes in each enriched term
- Larger figure size (16×12) for better readability
- Color-coded by source (GO:BP, GO:CC, REAC, TF, etc.)
- Reduced from top 30 to top 20 terms to avoid overcrowding

## Usage

### Configuration

Edit these parameters in the script:

```python
# Line 67-71
TOP_N_GENES = 100              # Number of genes to analyze
MIN_GENE_RECURRENCE = 3        # Min samples with SV
PVALUE_THRESHOLD = 0.05        # Significance threshold
USE_ENRICHR = True             # Enable/disable Enrichr API
```

### Running the Analysis

```bash
conda activate svmeta_env
python 05_pathway_enrichment.py
```

The script will now:
1. Load gene SV counts
2. Filter to top genes with recurrent SVs
3. Run g:Profiler enrichment (if available)
4. Run Enrichr enrichment (if enabled)
5. Identify GBM-relevant pathways from both sources
6. Analyze core GBM signaling pathway overlap
7. Generate enhanced visualizations
8. Export results

## Output Files

### New/Enhanced Files:

1. **enrichment_bubble.png** (Enhanced)
   - Now includes term and gene labels
   - Top 20 most significant terms
   - Color-coded by database source

2. **enrichr_results.csv** (New)
   - Full Enrichr enrichment results
   - All databases combined

3. **gbm_relevant_pathways_enrichr.csv** (New)
   - GBM-specific pathways from Enrichr
   - Filtered by canonical GBM keywords

4. **core_gbm_pathways.csv** (New)
   - Overlap with 6 core GBM signaling pathways
   - Enrichment ratios and gene lists

5. **core_gbm_pathways.png** (New)
   - Visual representation of core pathway enrichment
   - Bar plot showing enrichment ratios

### Existing Files (Unchanged):
- gprofiler_results.csv
- gbm_pathway_overlap.csv
- enrichment_overview.png
- gene_list.txt / .gmt / _ranked.rnk
- enrichment_report.txt

## Expected Results for GBM

### g:Profiler Results:
- 100-150 significant terms (GO:BP, GO:CC, REAC, TF, GO:MF)
- Expected GO terms: cell cycle, DNA repair, signal transduction

### Enrichr Results:
- 50-100 additional pathways from 5 databases
- GBM-specific: ~20-40 pathways matching canonical GBM keywords

### Core GBM Pathway Overlap:
You should see high enrichment for:
- **RTK/EGFR Signaling**: If EGFR, PDGFRA, MET have SVs
- **PI3K-AKT-mTOR**: If PTEN, PIK3CA have SVs (very common in GBM)
- **TP53 Pathway**: If TP53, MDM2, CDKN2A have SVs
- **RB/Cell Cycle**: If CDKN2A/B, CDK4, RB1 have SVs (very common)

Typical enrichment ratios for GBM cohorts:
- RTK/EGFR: 30-50% (3-5 out of 9 genes)
- PI3K-AKT-mTOR: 40-60% (4-6 out of 10 genes)
- TP53 Pathway: 40-60% (3-5 out of 8 genes)
- RB/Cell Cycle: 50-70% (5-7 out of 10 genes)

## Benefits Over Original Script

1. **Dual enrichment sources**: g:Profiler + Enrichr for comprehensive coverage
2. **GBM-focused filtering**: Automatic identification of GBM-relevant pathways
3. **Core pathway analysis**: Specific focus on the 6 key GBM signaling pathways
4. **Better visualizations**: Enhanced labels make plots publication-ready
5. **Expanded pathway database**: 19 canonical GBM pathways vs. 5 original

## Next Steps for Further Enhancement

If you want even more features, consider adding:

1. **Core GBM pathway visualization** (from mut_landscape script):
   - Side-by-side comparison of avg enrichment scores
   - Count of pathways per category

2. **Pathway category pie chart**:
   - Signal transduction, cell cycle, DNA repair, etc.
   - Shows distribution of functional categories

3. **Gene-pathway network table**:
   - Top pathways with their associated genes
   - Helps identify which genes drive which enrichments

4. **Database comparison**:
   - Bar plot comparing pathway counts across databases
   - Shows which databases are most informative

Would you like me to add any of these additional visualizations? Let me know and I can continue enhancing the script!

## References

### Merged from:
- [/home/chbope/extension/script/mut_landscape/gbm_pathway_analysis.py](../mut_landscape/gbm_pathway_analysis.py)

### GBM Pathway Classifications:
- Brennan et al. (2013). The Somatic Genomic Landscape of Glioblastoma. Cell.
- TCGA Research Network (2008). Comprehensive genomic characterization of GBM. Nature.
- Verhaak et al. (2010). Integrated genomic analysis of GBM subtypes. Cancer Cell.

### Enrichment Tools:
- g:Profiler: Raudvere et al. (2019). Nucleic Acids Research.
- Enrichr: Chen et al. (2013). BMC Bioinformatics.
