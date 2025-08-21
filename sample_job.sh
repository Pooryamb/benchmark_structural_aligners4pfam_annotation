#!/bin/bash
#SBATCH --time=02:00:00

module load StdEnv scipy-stack


touch tmp/timestamps/started_binning

mkdir -p data/processed/evalue_bins
file_paths=$(find ./tmp/alis/sample_pf/ -type f -name "*_B*.tsv")
echo "$file_paths" | parallel "python scripts/find_evalue_stratified_labels.py --input {}"

touch tmp/timestamps/ended_binning
