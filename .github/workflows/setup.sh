#!/bin/bash
set -euo pipefail

export SINGULARITY_VER=3.5.3

if type conda > /dev/null; then exit 0; fi
wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh -O anaconda.sh
bash anaconda.sh -b -p anaconda3
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda activate
conda install --channels conda-forge conda=4.9.2
conda env create --file envs/workflow.yaml
