#!/bin/bash

set -e

python3.7 ../scripts/prepare_config.py genomes/transcriptome.fasta --salmon-index-extra ' --kmerLen 5 ' --salmon-quant-extra ' --noBiasLengthThreshold --minAssignedFrags 1 --noEffectiveLengthCorrection --noLengthCorrection --fasterMapping --noFragLengthDist --allowDovetail --numPreAuxModelSamples 0 --numAuxModelSamples 0 ' --aggregate --libType "ISF"

python3.7 ../scripts/prepare_design.py reads/

snakemake -s ../Snakefile --use-conda -j 4 --force

snakemake -s ../Snakefile --use-conda --report -j 4 --force
