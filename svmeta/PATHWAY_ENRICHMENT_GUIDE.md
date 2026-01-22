## Summary

✅ **Pathway Enrichment Module Complete!**

I've created a comprehensive pathway enrichment analysis script that:

### What It Does:

1. **Automated Enrichment Analysis**
   - Uses g:Profiler API for GO, KEGG, Reactome enrichment
   - Analyzes known GBM-specific pathways
   - Generates publication-ready visualizations

2. **Multiple Data Sources**
   - Gene Ontology (Biological Process, Molecular Function, Cellular Component)
   - KEGG pathways
   - Reactome pathways
   - Transcription factor targets
   - Disease associations

3. **Exports for External Tools**
   - Plain text list (for DAVID, Enrichr)
   - GMT format (for GSEA)
   - Ranked list (by SV recurrence)

### Quick Start:

```bash
# Install g:Profiler (recommended but optional)
pip install gprofiler-official

# Run enrichment analysis
conda activate svmeta_env
python 05_pathway_enrichment.py
```

### Output Files:

**In `results/pathway_enrichment/`:**

1. **gprofiler_results.csv** - Full enrichment results
   - All significant GO terms, pathways, etc.
   - P-values, FDR, gene lists

2. **gbm_pathway_overlap.csv** - GBM-specific pathway analysis
   - RTK/RAS/PI3K pathway
   - TP53 pathway
   - Cell cycle pathway
   - DNA repair
   - Chromatin remodeling

3. **Visualizations:**
   - `enrichment_overview.png` - Top enriched terms by source
   - `enrichment_bubble.png` - Bubble plot of significant terms
   - `gbm_pathway_overlap.png` - GBM pathway bar chart

4. **Export Files:**
   - `gene_list.txt` - Simple list for DAVID/Enrichr
   - `gene_list.gmt` - GMT format for GSEA
   - `gene_list_ranked.rnk` - Ranked by recurrence

5. **enrichment_report.txt** - Comprehensive text report

### Expected Results for GBM:

You should see enrichment for:

**Top Pathways:**
- PI3K-Akt signaling pathway
- MAPK signaling pathway
- Cell cycle
- p53 signaling pathway
- RTK signaling
- Glioma pathway

**Top GO Terms:**
- Negative regulation of cell proliferation
- Regulation of cell cycle
- DNA damage response
- Apoptotic process
- Signal transduction

**GBM Pathway Overlap:**
- RTK/RAS/PI3K: 60-80% genes (EGFR, PTEN, PIK3CA, NF1, etc.)
- TP53 pathway: 40-60% genes (TP53, MDM2, CDKN2A)
- Cell cycle: 50-70% genes (CDKN2A, CDK4, RB1)

### Using External Tools:

#### DAVID (https://david.ncifcrf.gov/)
```bash
# Upload gene_list.txt
# Select: "Official Gene Symbol" as identifier
# Run Functional Annotation Clustering
```

#### Enrichr (https://maayanlab.cloud/Enrichr/)
```bash
# Paste contents of gene_list.txt
# Explore multiple enrichment libraries:
#   - KEGG 2021 Human
#   - GO Biological Process 2021
#   - WikiPathways 2021 Human
#   - DisGeNET
```

#### g:Profiler Web (https://biit.cs.ut.ee/gprofiler/)
```bash
# Paste gene list
# Organism: Homo sapiens
# Statistical domain scope: Only annotated genes
# Run query - get interactive results
```

### Troubleshooting:

**"gprofiler-official not installed"**
- Script will still run using local GBM pathway analysis
- Install with: `pip install gprofiler-official` for full analysis

**"No significant enrichments found"**
- Try lowering PVALUE_THRESHOLD in script (line 54)
- Check if TOP_N_GENES is too restrictive (line 49)
- Verify gene symbols are correct (HUGO names)

**"Gene recurrence file not found"**
- Run Step 3 first: `python 03_build_matrix_and_analyze.py`

### Advanced Usage:

**Adjust parameters in script:**

```python
# Line 49-54 - Analysis parameters
TOP_N_GENES = 100  # Increase for more genes
MIN_GENE_RECURRENCE = 3  # Lower for less strict filtering
PVALUE_THRESHOLD = 0.05  # Adjust significance
```

**Focus on specific gene set:**

```python
# After loading genes, filter manually
gene_list = ['EGFR', 'PTEN', 'CDKN2A', 'TP53', 'NF1']  # Custom list
enrichment_df = run_gprofiler_enrichment(gene_list)
```

### Integration with Other Steps:

```bash
# Complete workflow
./01_prepare_vcfs.sh
./02_merge_with_survivor.sh
python 03_build_matrix_and_analyze.py
python 04_external_dataset_comparison.py  # Optional
python 05_pathway_enrichment.py          # ← This step
```

### Publication-Ready Results:

The visualizations and tables are ready for:
- Supplementary figures
- Methods sections (cite g:Profiler)
- Results (pathway enrichment paragraph)
- Discussion (biological interpretation)

**Citation for g:Profiler:**
Raudvere et al. (2019). g:Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update). Nucleic Acids Research.

---

**You're all set!** Run the script and explore your pathway enrichment results! 🎉
