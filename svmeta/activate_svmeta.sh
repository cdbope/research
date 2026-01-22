#!/bin/bash
# Helper script to activate svmeta_env

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate svmeta_env

echo "svmeta_env activated!"
echo "SURVIVOR location: $(dirname $(dirname $(readlink -f "$0")))/tools/SURVIVOR/Debug/SURVIVOR"
