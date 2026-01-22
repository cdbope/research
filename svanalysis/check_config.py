#!/usr/bin/env python3
"""
Quick configuration check script
Run this to verify your configuration before running the full analysis
"""

import os
import sys

# Import configuration from sv_analysis_v2
try:
    import sv_analysis_v2 as sv2
    print("✓ sv_analysis_v2.py imported successfully\n")
except Exception as e:
    print(f"✗ Error importing sv_analysis_v2.py: {e}")
    sys.exit(1)

print("="*70)
print("CONFIGURATION CHECK")
print("="*70)

# Check paths
print("\n📁 Path Configuration:")
print(f"  SAMPLE_FILE:  {sv2.SAMPLE_FILE}")
print(f"  VCF_DIR:      {sv2.VCF_DIR}")
print(f"  OUTPUT_DIR:   {sv2.OUTPUT_DIR}")

# Check if files/directories exist
print("\n🔍 Existence Check:")
errors = []

if os.path.exists(sv2.SAMPLE_FILE):
    print(f"  ✓ SAMPLE_FILE exists")
else:
    print(f"  ✗ SAMPLE_FILE not found")
    errors.append("SAMPLE_FILE")

if os.path.exists(sv2.VCF_DIR):
    print(f"  ✓ VCF_DIR exists")
    vcf_count = len([f for f in os.listdir(sv2.VCF_DIR) if f.endswith('.vcf.gz')])
    print(f"    Found {vcf_count} VCF files (.vcf.gz)")
else:
    print(f"  ✗ VCF_DIR not found")
    errors.append("VCF_DIR")

if os.path.exists(sv2.OUTPUT_DIR):
    print(f"  ✓ OUTPUT_DIR exists")
else:
    print(f"  ⚠ OUTPUT_DIR doesn't exist (will be created)")

# Check parameters
print("\n⚙️ Quality Filtering Parameters:")
print(f"  MIN_SUPPORT:      {sv2.MIN_SUPPORT}")
print(f"  MIN_VAF:          {sv2.MIN_VAF}")
print(f"  CLONAL_THRESHOLD: {sv2.CLONAL_THRESHOLD}")

# Final status
print("\n" + "="*70)
if errors:
    print("❌ CONFIGURATION HAS ERRORS")
    print(f"   Please fix: {', '.join(errors)}")
    print("   See CONFIGURATION_GUIDE.md for help")
    sys.exit(1)
else:
    print("✅ CONFIGURATION LOOKS GOOD!")
    print("   Ready to run: python sv_analysis_v2.py")
print("="*70)
