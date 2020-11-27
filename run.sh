#!/bin/bash

# Install conda environment if not installed before
$(conda info --envs | grep "rna-count-salmon" --quiet) && echo "Pipeline already installed! What a chance!" || conda env create --force -f "/mnt/beegfs/pipelines/rna-count-salmon/pipeline/rna-count-salmon/envs/workflow_flamingo.yaml"

# Check on environment variables: if env are missing,
# then installation process did not work properly
conda activate rna-count-salmon
$(export -p | grep "RNA_COUNT_LAUNCHER" --quiet) && python3 ${RNA_COUNT_LAUNCHER} flamingo || echo "Magic did not work. :-("
