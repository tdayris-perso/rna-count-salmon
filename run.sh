#!/bin/bash
set -eui

# Install conda environment if not installed before
echo "Installing environment if needed."
$(conda info --envs | grep "rna-count-salmon" --quiet) && echo "Pipeline already installed! What a chance!" || conda env create --force -f "/mnt/beegfs/pipelines/rna-count-salmon/pipeline/rna-count-salmon/envs/workflow_flamingo.yaml"

# Check on environment variables: if env are missing
echo "Loading environment"
conda activate rna-count-salmon

# then installation process did not work properly
echo "Running pipeline if possible"
$(export -p | grep "RNA_COUNT_LAUNCHER" --quiet) && python3 ${RNA_COUNT_LAUNCHER} flamingo || echo "Magic did not work. :-("
