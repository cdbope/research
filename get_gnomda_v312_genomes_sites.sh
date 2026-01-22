#!/usr/bin/env bash
set -euo pipefail

# ------------------ CONFIG ------------------
# Public HTTP URL for gnomAD v3.1.2 genomes (VCF, HGDP+TGP subset)
# Using HTTP URL instead of GCS for direct downloads
HTTP_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes"

# Output filenames
OUT_VCF="genomes-v3.1.2.sites.vcf.bgz"

# Chromosomes to fetch (adjust if you also want chrMT)
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
        chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

# Working directory for per-chrom files (created if missing)
WORKDIR="gnomad_v3.1.2_bychr"
# --------------------------------------------

# ---- sanity checks
command -v wget >/dev/null 2>&1 || { echo "ERROR: wget not found"; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix not found"; exit 1; }

mkdir -p "${WORKDIR}"

echo "==> Downloading per-chromosome VCFs to ${WORKDIR}"
for CHR in "${CHROMS[@]}"; do
  FILE="gnomad.genomes.v3.1.2.hgdp_tgp.${CHR}.vcf.bgz"
  IDX="${FILE}.tbi"
  SRC_VCF="${HTTP_BASE}/${FILE}"
  SRC_TBI="${HTTP_BASE}/${IDX}"
  
  # Skip if file already exists and is non-empty
  if [[ -s "${WORKDIR}/${FILE}" ]]; then
    echo "Skipping ${FILE} (already exists)"
  else
    echo "Downloading ${FILE}..."
    wget -q --show-progress -O "${WORKDIR}/${FILE}" "${SRC_VCF}" || { echo "Failed: ${SRC_VCF}"; exit 1; }
  fi
  
  if [[ -s "${WORKDIR}/${IDX}" ]]; then
    echo "Skipping ${IDX} (already exists)"
  else
    echo "Downloading ${IDX}..."
    wget -q --show-progress -O "${WORKDIR}/${IDX}" "${SRC_TBI}" || { echo "Failed: ${SRC_TBI}"; exit 1; }
  fi
done

echo "==> Verifying all files present"
missing=0
for CHR in "${CHROMS[@]}"; do
  FILE="${WORKDIR}/gnomad.genomes.v3.1.2.hgdp_tgp.${CHR}.vcf.bgz"
  IDX="${FILE}.tbi"
  [[ -s "${FILE}" && -s "${IDX}" ]] || { echo "Missing: ${FILE} or ${IDX}"; missing=1; }
done
[[ $missing -eq 0 ]] || { echo "ERROR: Missing files. Aborting."; exit 1; }

echo "==> Merging chromosomes into one VCF (bgzipped)"
# bcftools concat keeps headers properly and ensures valid output
# -Oz: bgzip output; --threads: parallelize compression if supported
BCFTOOLS_THREADS=${BCFTOOLS_THREADS:-4}

# Build list of input VCFs in correct order
INPUTS=()
for CHR in "${CHROMS[@]}"; do
  INPUTS+=("${WORKDIR}/gnomad.genomes.v3.1.2.hgdp_tgp.${CHR}.vcf.bgz")
done

bcftools concat \
  --threads "${BCFTOOLS_THREADS}" \
  -Oz -o "${OUT_VCF}" \
  "${INPUTS[@]}"

echo "==> Indexing merged VCF with tabix"
tabix -f -p vcf "${OUT_VCF}"

# Optional: also create CSI index (better for huge files)
# bcftools index --threads "${BCFTOOLS_THREADS}" -f -c "${OUT_VCF}"

echo "==> Done."
echo "Merged VCF: ${OUT_VCF}"
echo "Index:      ${OUT_VCF}.tbi"

