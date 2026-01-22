#!/bin/bash

# Clean up intermediate files
echo "Cleaning up intermediate files..."
sudo rm -rf /home/chbope/extension/data/200GMBs/pcgr/starsigndna/refit/tmp/intermediate_results_dir/*
sudo rm -rf /home/chbope/extension/data/200GMBs/pcgr/starsigndna/refit/tmp/logs/*

# Run deepsomatic
echo "Running DeepSomatic..."
cd /home/chbope/extension/script/deepsomatic
bash deepsomatic.sh
