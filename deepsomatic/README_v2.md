# DeepSomatic Pipeline v2 with ANNOVAR Annotation

This enhanced version of the DeepSomatic pipeline includes automated variant annotation using ANNOVAR.

## Pipeline Overview

The pipeline performs the following steps:

1. **Variant Calling**: Runs DeepSomatic for somatic variant detection
2. **Filtering**: Filters VCF to keep only PASS variants
3. **Format Conversion**: Converts VCF to ANNOVAR input format
4. **Annotation**: Annotates variants using ANNOVAR with:
   - RefGene (gene-based annotation)
   - ClinVar (clinical significance)
   - COSMIC (cancer mutations database)
5. **Post-filtering**: Filters annotated variants based on:
   - Exonic nonsynonymous variants (excluding benign)
   - Upstream variants
   - TERT promoter mutations
   - Distance-based filtering (excludes dist=166)

## Prerequisites

1. **Docker** - For running DeepSomatic
2. **bcftools** - For VCF filtering
3. **ANNOVAR** - For variant annotation
   - Download from: http://www.openbioinformatics.org/annovar/
   - Required databases: refGene, clinvar_20240611, cosmic100coding2024

## Setup Instructions

### 1. Install ANNOVAR (if not already installed)

```bash
# Download ANNOVAR (requires registration)
# Visit: http://www.openbioinformatics.org/annovar/annovar_download_form.php

# Extract ANNOVAR
tar -xzvf annovar.latest.tar.gz

# Download required databases
cd annovar
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240611 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic100coding2024 humandb/
```

### 2. Configure the Pipeline

Edit [config.sh](config.sh) and update the following variables:

```bash
# Update these paths to match your system
ANNOVAR_DIR="/path/to/annovar"
HUMANDB_DIR="/path/to/annovar/humandb"

# Update sample information
SAMPLE_ID="your_sample_id"
BAM_FILE="your_sample.bam"

# Optional: Update other parameters as needed
MODEL_TYPE="ONT_TUMOR_ONLY"  # or FFPE_WGS_TUMOR_ONLY, FFPE_WES_TUMOR_ONLY
```

### 3. Install bcftools (if not already installed)

```bash
# Ubuntu/Debian
sudo apt-get install bcftools

# Or conda
conda install -c bioconda bcftools
```

## Usage

### Basic Usage

```bash
cd /home/chbope/extension/script/deepsomatic
bash deepsomatic_v2.sh
```

### Running with Custom Configuration

```bash
# Edit config.sh first
nano config.sh

# Then run the pipeline
bash deepsomatic_v2.sh
```

## Output Files

The pipeline generates the following output files in `OUTPUT_DIR`:

1. **Raw VCF**: `${SAMPLE_ID}_ont_deepsomatic_output.vcf.gz`
   - All variants called by DeepSomatic

2. **Filtered VCF**: `${SAMPLE_ID}_ont_deepsomatic_output.pass.vcf.gz`
   - Only PASS-quality variants

3. **ANNOVAR Input**: `deepsomatic_To_snv_avinput`
   - Converted format for ANNOVAR

4. **Full Annotation**: `deepsomatic_To_snv_avinput_snv.hg38_multianno.txt`
   - Complete ANNOVAR annotation results

5. **Final Filtered CSV**: `${SAMPLE_ID}_annotateandfilter_deep_somatic.csv`
   - Filtered and prioritized variants (columns 1-16, 25, 26)

## Filtering Criteria

The final output applies the following filters:

- **Include**:
  - Exonic nonsynonymous variants (excluding benign)
  - Upstream variants
  - TERT promoter region variants
  - Exonic variants

- **Exclude**:
  - Benign variants
  - Variants with dist=166

## Model Types

Available DeepSomatic models (set in `config.sh`):

- `ONT_TUMOR_ONLY` - Oxford Nanopore tumor-only sequencing
- `FFPE_WGS_TUMOR_ONLY` - FFPE whole genome sequencing (tumor only)
- `FFPE_WES_TUMOR_ONLY` - FFPE whole exome sequencing (tumor only)
- `WGS_TUMOR_ONLY` - Fresh frozen WGS (tumor only)
- `WES_TUMOR_ONLY` - Fresh frozen WES (tumor only)

## Troubleshooting

### Error: "Could not find TensorRT"
- This is a warning, not an error. DeepSomatic will run on CPU without GPU acceleration.

### Error: "could not load fasta and/or fai"
- Ensure the reference directory is properly mounted
- Check that the .fai index file exists alongside the FASTA file

### Error: "ANNOVAR directory not found"
- Update `ANNOVAR_DIR` in config.sh to point to your ANNOVAR installation

### Error: "bcftools: command not found"
- Install bcftools: `sudo apt-get install bcftools`

## Version History

- **v2**: Added ANNOVAR annotation pipeline
- **v1**: Basic DeepSomatic variant calling

## References

- DeepSomatic: https://github.com/google/deepvariant/blob/r1.9/docs/deepsomatic-details.md
- ANNOVAR: http://annovar.openbioinformatics.org/
