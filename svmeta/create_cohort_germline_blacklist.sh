#!/bin/bash
###############################################################################
# Create Cohort-Specific Germline Blacklist
#
# Logic: True somatic SVs are sample-specific (chromothripsis is unique).
#        Germline SVs appear in many samples at similar positions.
#        SVs present in >30% of cohort are likely germline heterozygous.
#
# Usage: ./create_cohort_germline_blacklist.sh
###############################################################################

set -euo pipefail

# Configuration
INPUT_VCF_DIR="results/prepared_vcfs/filter_vcf"
OUTPUT_DIR="results/cohort_germline_analysis"
SURVIVOR="/path/to/SURVIVOR"  # Update with your SURVIVOR path
MIN_SAMPLES=60  # 30% of 200 samples (adjust as needed)

mkdir -p "$OUTPUT_DIR"

echo "============================================================================"
echo "COHORT-SPECIFIC GERMLINE BLACKLIST GENERATION"
echo "============================================================================"
echo "Input directory: $INPUT_VCF_DIR"
echo "Minimum samples for blacklist: $MIN_SAMPLES (30% of cohort)"
echo "============================================================================"
echo

# Step 1: Create sample list for SURVIVOR
echo "Step 1: Creating sample list..."
ls "$INPUT_VCF_DIR"/*.vcf > "$OUTPUT_DIR/sample_files.txt"
n_samples=$(wc -l < "$OUTPUT_DIR/sample_files.txt")
echo "  Found $n_samples VCF files"

# Step 2: Merge all samples to find recurrent SVs
echo "Step 2: Merging samples with SURVIVOR..."
$SURVIVOR merge \
  "$OUTPUT_DIR/sample_files.txt" \
  1000 \
  2 \
  1 \
  1 \
  0 \
  50 \
  "$OUTPUT_DIR/cohort_merged.vcf"

echo "  ✓ Merged cohort VCF created"

# Step 3: Analyze frequency of each SV across cohort
echo "Step 3: Identifying recurrent SVs (present in >${MIN_SAMPLES} samples)..."

python3 - <<'EOF'
import sys

vcf_file = "results/cohort_germline_analysis/cohort_merged.vcf"
min_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 60
output_bed = "results/cohort_germline_analysis/germline_blacklist.bed"
output_vcf = "results/cohort_germline_analysis/germline_blacklist.vcf"

germline_svs = []
total_svs = 0

with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('##'):
            continue

        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            n_total_samples = len(header) - 9
            print(f"Total samples in merged VCF: {n_total_samples}")
            continue

        total_svs += 1
        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]

        # Count samples with this SV (non-reference genotype)
        genotypes = fields[9:]
        n_samples_with_sv = sum(1 for gt in genotypes if not gt.startswith('0/0') and not gt.startswith('./.'))

        freq = n_samples_with_sv / n_total_samples if n_total_samples > 0 else 0

        # Flag as likely germline if present in many samples
        if n_samples_with_sv >= min_samples:
            # Get SV end position from INFO field
            info = fields[7]
            end_pos = pos + 1000  # Default
            for item in info.split(';'):
                if item.startswith('END='):
                    end_pos = int(item.split('=')[1])

            germline_svs.append({
                'chrom': chrom,
                'start': pos,
                'end': end_pos,
                'ref': ref,
                'alt': alt,
                'n_samples': n_samples_with_sv,
                'freq': freq,
                'line': line.strip()
            })

# Write BED file for filtering
with open(output_bed, 'w') as f:
    f.write("# Germline blacklist: SVs present in ≥{} samples ({:.1f}%)\n".format(
        min_samples, 100 * min_samples / n_total_samples))
    for sv in germline_svs:
        f.write(f"{sv['chrom']}\t{sv['start']}\t{sv['end']}\t"
                f"n={sv['n_samples']};freq={sv['freq']:.2f}\n")

# Write VCF file for inspection
with open(output_vcf, 'w') as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write(f"##source=cohort_germline_blacklist\n")
    f.write(f"##INFO=<ID=NSAMP,Number=1,Type=Integer,Description=\"Number of samples with this SV\">\n")
    f.write(f"##INFO=<ID=FREQ,Number=1,Type=Float,Description=\"Frequency across cohort\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for sv in germline_svs:
        f.write(f"{sv['chrom']}\t{sv['start']}\t.\t{sv['ref']}\t{sv['alt']}\t.\t.\t"
                f"NSAMP={sv['n_samples']};FREQ={sv['freq']:.3f}\n")

print(f"\n{'='*80}")
print(f"RESULTS:")
print(f"{'='*80}")
print(f"Total SVs in merged cohort: {total_svs}")
print(f"Likely germline SVs (≥{min_samples} samples): {len(germline_svs)}")
print(f"Percentage flagged as germline: {100 * len(germline_svs) / total_svs:.1f}%")
print(f"\nOutput files:")
print(f"  BED blacklist: {output_bed}")
print(f"  VCF for inspection: {output_vcf}")
print(f"{'='*80}")

# Show top 10 most frequent
print(f"\nTop 10 most frequent SVs (likely germline):")
germline_svs_sorted = sorted(germline_svs, key=lambda x: x['n_samples'], reverse=True)
for i, sv in enumerate(germline_svs_sorted[:10], 1):
    print(f"  {i}. {sv['chrom']}:{sv['start']}-{sv['end']} | "
          f"{sv['n_samples']} samples ({100*sv['freq']:.1f}%)")

EOF

echo
echo "============================================================================"
echo "Cohort germline blacklist generation complete!"
echo "============================================================================"
echo "Next step: Filter individual samples with apply_cohort_blacklist.sh"
