#!/bin/bash

set -e

python3.7 ../scripts/prepare_config.py genome/transcriptome.fasta --salmon_index_extra '\-\-kmerLen 5' --salmon_quant_extra '\-\-noBiasLengthThreshold \-\-minAssignedFrags 1 \-\-noEffectiveLengthCorrection \-\-noLengthCorrection \-\-fasterMapping \-\-noFragLengthDist \-\-allowDovetail \-\-numPreAuxModelSamples 0 \-\-numAuxModelSamples 0' --aggregate --libType "ISF"

python3.7 ../scripts/prepare_design.py reads/

snakemake -s ../Snakefile --use-conda -j 4

snakemake -s ../Snakefile --use-conda --report -j 4
