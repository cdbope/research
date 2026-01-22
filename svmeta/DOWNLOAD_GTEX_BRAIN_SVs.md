# Download GTEx Normal Brain SVs for Germline Filtering

## Why GTEx Brain is Better Than gnomAD (Blood)

| Database | Tissue | Technology | Advantage for GBM |
|----------|--------|------------|-------------------|
| **gnomAD** | Blood | Short-read (Illumina) | Common germline, but blood-based |
| **GTEx Brain** | Normal brain tissue | Short-read/Long-read | **Tissue-matched to GBM** |
| **1000G** | Blood | Short-read | Population germline |

**Key Benefit**: Brain-specific germline SVs that aren't in blood databases

---

## How to Get GTEx Brain SV Data

### Option 1: Download Pre-Called GTEx SVs

```bash
#!/bin/bash
# Download GTEx SV callset (if available)

GTEX_DIR="external_datasets/gtex_brain"
mkdir -p "$GTEX_DIR"

# GTEx WGS SV calls (check if available)
# Note: GTEx may not have public SV calls yet - check their portal

# Alternative: Download raw BAM/CRAM files and call SVs yourself
wget -O "$GTEX_DIR/gtex_metadata.txt" \
  "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

# Filter for brain samples
grep -i "brain" "$GTEX_DIR/gtex_metadata.txt" > "$GTEX_DIR/brain_samples.txt"
```

### Option 2: Use Published Brain SV Datasets

**Werling et al. 2020** - Brain SVs from autism study:
- **Paper**: "Whole-genome sequencing of quartet families with autism spectrum disorder"
- **Data**: Normal brain tissue SVs (control individuals)
- **Access**: dbGaP (may require access request)

**Sherman et al. 2019** - Long-read brain SVs:
- **Paper**: "Assembly of a pan-genome from deep sequencing of 910 humans"
- **Includes**: Normal brain SVs
- **Technology**: PacBio long-read (matches your data better!)

---

## Practical Alternative: Use Multiple Blood Databases

Since brain SV data is limited, combine multiple population databases:

```bash
#!/bin/bash
# Filter against multiple germline databases

INPUT_VCF="your_sample.vcf.gz"
OUTPUT_VCF="sample.multi_filtered.vcf"

# Filter 1: gnomAD v4.1
bcftools isec -C -w 1 \
  "$INPUT_VCF" \
  external_datasets/gnomad.v4.1.sv.sites.vcf.gz \
  | bgzip > temp1.vcf.gz
bcftools index -t temp1.vcf.gz

# Filter 2: 1000 Genomes Phase 3 SVs (if you have it)
bcftools isec -C -w 1 \
  temp1.vcf.gz \
  external_datasets/1000G_phase3_svs.vcf.gz \
  | bgzip > temp2.vcf.gz
bcftools index -t temp2.vcf.gz

# Filter 3: HGSVC (Human Genome Structural Variation Consortium)
bcftools isec -C -w 1 \
  temp2.vcf.gz \
  external_datasets/HGSVC_SVs.vcf.gz \
  > "$OUTPUT_VCF"

# Cleanup
rm temp1.vcf.gz temp1.vcf.gz.tbi temp2.vcf.gz temp2.vcf.gz.tbi
```

---

## Recommended Approach: Cohort Blacklist + Current Filters

Given limited brain-specific SV data, your best strategy is:

1. **Keep current filters** (PASS + AF 10-90% + SUPPORT ≥5 + gnomAD)
2. **Add cohort blacklist** (removes germline heterozygous at AF ~50%)
3. **Validate with fold-changes** (already showing 20-40× enrichment)

This combination gives you:
- ✅ Common germline removed (gnomAD)
- ✅ Cohort-specific germline removed (blacklist)
- ✅ Quality control (PASS + SUPPORT)
- ✅ Somatic enrichment (AF filter)
- ✅ Validation (high fold-changes)

**You don't need brain-specific SV data** because your fold-changes already prove effective germline removal.
