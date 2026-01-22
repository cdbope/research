#!/bin/bash
###############################################################################
# Identify Outlier Samples with High SV Counts
#
# Logic: Samples with abnormally high SV counts may have:
#        - Technical issues (poor sequencing quality)
#        - High germline contamination
#        - True hypermutator phenotype (rare)
#
# Usage: ./identify_outlier_samples.sh
###############################################################################

set -euo pipefail

INPUT_VCF_DIR="results/prepared_vcfs/filter_vcf"
OUTPUT_DIR="results/qc"

mkdir -p "$OUTPUT_DIR"

echo "============================================================================"
echo "SAMPLE-LEVEL SV COUNT ANALYSIS"
echo "============================================================================"

python3 - <<'EOF'
import os
import statistics

vcf_dir = "results/prepared_vcfs/filter_vcf"
output_file = "results/qc/sample_sv_counts.txt"

# Count SVs per sample
sample_counts = {}
for vcf_file in os.listdir(vcf_dir):
    if not vcf_file.endswith('.vcf'):
        continue

    sample_name = vcf_file.replace('.wf_sv.vcf', '')
    vcf_path = os.path.join(vcf_dir, vcf_file)

    # Count variants
    count = 0
    with open(vcf_path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                count += 1

    sample_counts[sample_name] = count

# Calculate statistics
counts = list(sample_counts.values())
mean_count = statistics.mean(counts)
median_count = statistics.median(counts)
stdev_count = statistics.stdev(counts) if len(counts) > 1 else 0

# Identify outliers (>2 std dev above mean)
outlier_threshold = mean_count + (2 * stdev_count)

# Write results
with open(output_file, 'w') as f:
    f.write("Sample\tSV_Count\tStatus\n")
    for sample, count in sorted(sample_counts.items(), key=lambda x: x[1], reverse=True):
        status = "OUTLIER" if count > outlier_threshold else "NORMAL"
        f.write(f"{sample}\t{count}\t{status}\n")

# Print summary
print(f"Sample SV Count Statistics:")
print(f"  Total samples: {len(counts)}")
print(f"  Mean SV count: {mean_count:.1f}")
print(f"  Median SV count: {median_count:.1f}")
print(f"  Std deviation: {stdev_count:.1f}")
print(f"  Outlier threshold (mean + 2×SD): {outlier_threshold:.1f}")
print()

outliers = [(s, c) for s, c in sample_counts.items() if c > outlier_threshold]
if outliers:
    print(f"Outlier samples ({len(outliers)}):")
    for sample, count in sorted(outliers, key=lambda x: x[1], reverse=True):
        print(f"  {sample}: {count} SVs ({count/mean_count:.1f}× mean)")
    print()
    print("⚠️  Consider excluding these samples or investigating for:")
    print("    - Sequencing quality issues")
    print("    - Germline contamination")
    print("    - True hypermutator phenotype (rare)")
else:
    print("✓ No outlier samples detected")

print(f"\nFull results: {output_file}")

EOF

echo "============================================================================"
echo "Sample QC analysis complete!"
echo "Results: results/qc/sample_sv_counts.txt"
echo "============================================================================"
