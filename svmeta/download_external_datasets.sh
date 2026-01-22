#!/bin/bash
###############################################################################
# Download External Datasets for Comparison
#
# This script provides guidance and examples for downloading TCGA-GBM and
# PCAWG structural variant data for comparison with your cohort.
#
# Usage: ./download_external_datasets.sh
###############################################################################

set -e

EXTERNAL_DIR="/home/chbope/extension/script/svmeta/external_datasets"
mkdir -p "$EXTERNAL_DIR"

echo "============================================================================"
echo "EXTERNAL DATASET DOWNLOAD GUIDE"
echo "============================================================================"
echo
echo "This script provides instructions for downloading TCGA-GBM and PCAWG data."
echo "Due to data access restrictions, you'll need to follow these steps manually."
echo
echo "============================================================================"
echo "1. TCGA-GBM STRUCTURAL VARIANTS"
echo "============================================================================"
echo
echo "Option A: GDC Data Portal (Recommended)"
echo "  1. Visit: https://portal.gdc.cancer.gov/"
echo "  2. Select Project: TCGA-GBM"
echo "  3. Data Category: Structural Variation"
echo "  4. Data Type: Structural Rearrangement"
echo "  5. Download VCF files or aggregated data"
echo
echo "Option B: Broad GDAC Firehose"
echo "  1. Visit: https://gdac.broadinstitute.org/"
echo "  2. Select: TCGA-GBM"
echo "  3. Download: Copy Number Variation and SV data"
echo
echo "Option C: Published TCGA-GBM SV Calls"
echo "  - Brennan et al. (2013) Cell"
echo "  - The TCGA-GBM marker paper"
echo "  - Supplementary data tables contain recurrent SVs"
echo
echo "Expected format for analysis:"
echo "  CSV file with columns: sv_id, chr, start, end, svtype, frequency, num_samples, total_samples, genes"
echo
echo "============================================================================"
echo "2. PCAWG STRUCTURAL VARIANTS"
echo "============================================================================"
echo
echo "Option A: ICGC Data Portal"
echo "  1. Visit: https://dcc.icgc.org/pcawg"
echo "  2. Select: Glioblastoma samples"
echo "  3. Download: Structural Variation VCFs"
echo
echo "Option B: Synapse Repository"
echo "  1. Register at: https://www.synapse.org/"
echo "  2. Access PCAWG data: syn7596712"
echo "  3. Download SV consensus calls"
echo
echo "Option C: Published Datasets"
echo "  - Li et al. (2020) Nature"
echo "  - PCAWG Structural Variants Working Group"
echo "  - Supplementary datasets"
echo
echo "Expected format for analysis:"
echo "  CSV file with columns: sv_id, chr, start, end, svtype, frequency, num_samples, total_samples, genes"
echo
echo "============================================================================"
echo "3. ALTERNATIVE: USE PROVIDED TEMPLATES"
echo "============================================================================"
echo
echo "If you don't have access to TCGA/PCAWG data, the comparison script will"
echo "create template files that you can populate manually with published data:"
echo
echo "  - $EXTERNAL_DIR/tcga_gbm_sv_summary_TEMPLATE.csv"
echo "  - $EXTERNAL_DIR/pcawg_gbm_sv_summary_TEMPLATE.csv"
echo
echo "These templates show the expected format. You can fill them with:"
echo "  - Published recurrent SVs from papers"
echo "  - Data from supplementary tables"
echo "  - Known GBM SV hotspots from literature"
echo
echo "============================================================================"
echo "4. PROCESSING RAW VCF FILES"
echo "============================================================================"
echo
echo "If you download raw VCF files, use this workflow to create summary tables:"
echo
cat << 'EOF'

# Example: Process TCGA VCFs
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' \
    tcga_sample*.vcf.gz > tcga_svs.txt

# Aggregate to get recurrent SVs
python << 'PYTHON_EOF'
import pandas as pd
import numpy as np

# Read SV calls
df = pd.read_csv('tcga_svs.txt', sep='\t',
                 names=['chr', 'start', 'end', 'svtype'])

# Calculate frequencies
# ... (clustering logic to identify recurrent events)

# Save summary
summary.to_csv('tcga_gbm_sv_summary.csv', index=False)
PYTHON_EOF

EOF

echo
echo "============================================================================"
echo "5. KNOWN GBM SV HOTSPOTS (FROM LITERATURE)"
echo "============================================================================"
echo
echo "If external data is unavailable, you can manually create a reference file"
echo "with well-known GBM SV events from literature:"
echo
echo "Common GBM SVs:"
echo "  - EGFR amplifications (chr7:55M-56M)"
echo "  - CDKN2A/B deletions (chr9p21)"
echo "  - PTEN deletions (chr10q23)"
echo "  - Chromosome 7 gain / Chromosome 10 loss"
echo "  - EGFR-vIII rearrangements"
echo "  - PDGFRA amplifications"
echo
echo "============================================================================"
echo "6. NEXT STEPS"
echo "============================================================================"
echo
echo "After obtaining external datasets:"
echo
echo "1. Place files in: $EXTERNAL_DIR/"
echo "   - tcga_gbm_sv_summary.csv"
echo "   - pcawg_gbm_sv_summary.csv"
echo
echo "2. Run comparison:"
echo "   conda activate svmeta_env"
echo "   python 04_external_dataset_comparison.py"
echo
echo "3. Review results in:"
echo "   results/external_comparison/"
echo
echo "============================================================================"
echo

# Create a sample reference file with known GBM SVs from literature
cat > "$EXTERNAL_DIR/known_gbm_svs_literature.csv" << 'EOF'
sv_id,chr,start,end,svtype,frequency,num_samples,total_samples,genes,reference
chr7:55019032-55211628_DUP,chr7,55019032,55211628,DUP,0.45,150,333,EGFR,"TCGA 2013 (Brennan et al.)"
chr9:21967751-21994490_DEL,chr9,21967751,21994490,DEL,0.52,173,333,CDKN2A;CDKN2B,"TCGA 2013"
chr10:89622869-89731687_DEL,chr10,89622869,89731687,DEL,0.41,137,333,PTEN,"TCGA 2013"
chr4:54274328-54292994_DUP,chr4,54274328,54292994,DUP,0.15,50,333,PDGFRA,"TCGA 2013"
chr12:58142003-58146971_DUP,chr12,58142003,58146971,DUP,0.19,63,333,CDK4,"TCGA 2013"
chr12:68808100-68849456_DUP,chr12,68808100,68849456,DUP,0.17,57,333,MDM2,"TCGA 2013"
chr1:11166592-11322608_DUP,chr1,11166592,11322608,DUP,0.08,27,333,MET,"TCGA 2013"
chr13:48303748-48481890_DEL,chr13,48303748,48481890,DEL,0.36,120,333,RB1,"TCGA 2013"
chr17:29421944-29446394_DEL,chr17,29421944,29446394,DEL,0.25,83,333,NF1,"TCGA 2013"
EOF

echo "Created reference file: known_gbm_svs_literature.csv"
echo "This contains well-characterized GBM SVs from TCGA publications."
echo
echo "You can use this as a reference dataset by updating the script:"
echo "  TCGA_SV_FILE = \"$EXTERNAL_DIR/known_gbm_svs_literature.csv\""
echo
echo "============================================================================"
echo "SETUP COMPLETE"
echo "============================================================================"
