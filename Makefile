SHELL := bash
.ONESHELL:
.SHELLFLAGS := -euio pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

### Variables ###
# Tools
PYTEST           = pytest
BASH             = bash
CONDA            = conda
PYTHON           = python3.8
SNAKEMAKE        = snakemake
CONDA_ACTIVATE   = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

# Paths
TEST_CONFIG      = scripts/prepare_config.py
TEST_DESIGN      = scripts/prepare_design.py
TEST_COMMON      = scripts/common_script_rna_count_salmon.py
SNAKE_FILE       = Snakefile
ENV_YAML         = envs/workflow.yaml
TRANSCRIPT_PATH  = tests/genome/transcriptome.fasta
GTF_PATH         = tests/genome/annot.gtf
READS_PATH       = tests/reads/

# Arguments
PYTEST_ARGS      = -vv
ENV_NAME         = rna-count-salmon
SNAKE_THREADS    = 1
SAINDEX_ARGS     = ' --kmerLen 5 '
SAQUANT_ARGS     = ' --noBiasLengthThreshold --minAssignedFrags 1 --noEffectiveLengthCorrection --noLengthCorrection --fasterMapping --noFragLengthDist --allowDovetail --numPreAuxModelSamples 0 --numAuxModelSamples 0 '

# Recipes
default: all-unit-tests

# Environment building through conda
conda-tests:
	${CONDA_ACTIVATE} base && \
	${CONDA} env create --file ${ENV_YAML} --force && \
	${CONDA} activate ${ENV_NAME}
.PHONY: conda-tests


### UNIT TESTS ###
# Running all unit-tests (one for each python scripts)
all-unit-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} ${PYTEST_ARGS} ${TEST_CONFIG} ${TEST_DESIGN} ${TEST_COMMON}
.PHONY: all-unit-tests


# Running all unit test (on prepare_config.py only)
config-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} ${PYTEST_ARGS} ${TEST_CONFIG} && \
	${PYTHON} ${TEST_CONFIG} ${TRANSCRIPT_PATH} --salmon-index-extra ${SAINDEX_ARGS} --salmon-quant-extra ${SAQUANT_ARGS} --aggregate --libType ISF --workdir tests --debug
.PHONY: config-tests


# Running all unit test (on prepare_design.py only)
design-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} ${PYTEST_ARGS} ${TEST_DESIGN} && \
	${PYTHON} ${TEST_DESIGN} ${READS_PATH} -o tests/design.tsv --debug
.PHONY: design-tests


# Running all unit test (on common_scrit_rna_count_salmon.py only)
common-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} ${PYTEST_ARGS} ${TEST_COMMON}
.PHONY: common-tests


### Continuous Integration Tests ###
# Running snakemake on test datasets
test-conda-report.html:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTHON} ${TEST_DESIGN} ${READS_PATH} -o ${PWD}/tests/design.tsv --debug && \
	${PYTHON} ${TEST_CONFIG} ${TRANSCRIPT_PATH} ${GTF_PATH} --salmon-index-extra ${SAINDEX_ARGS} --salmon-quant-extra ${SAQUANT_ARGS} --aggregate --libType ISF --workdir ${PWD}/tests --design ${PWD}/tests/design.tsv --threads ${SNAKE_THREADS} --debug && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --configfile ${PWD}/tests/config.yaml --forceall --printshellcmds --reason --directory ${PWD}/tests && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --configfile ${PWD}/tests/config.yaml --directory ${PWD}/tests --report test-conda-report.html


# Running snakemake on test datasets with singularity flag raised on
test-singularity-report.html:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTHON} ${TEST_DESIGN} ${READS_PATH} -o ${PWD}/tests/design.tsv --debug && \
	${PYTHON} ${TEST_CONFIG} ${TRANSCRIPT_PATH} ${GTF_PATH} --salmon-index-extra ${SAINDEX_ARGS} --salmon-quant-extra ${SAQUANT_ARGS} --aggregate --libType ISF --workdir ${PWD}/tests --design ${PWD}/tests/design.tsv --threads ${SNAKE_THREADS} --debug && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --configfile ${PWD}/tests/config.yaml --forceall --printshellcmds --reason --directory ${PWD}/tests --use-singularity && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --configfile ${PWD}/tests/config.yaml --directory ${PWD}/tests --report test-singularity-report.html

# Cleaning Snakemake outputs
clean:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --force --configfile ${PWD}/tests/config.yaml --directory ${PWD}/tests --delete-all-output && \
	rm -r ${PWD}/tests/*-report.html
.PHONY: clean

# Display pipeline graph
workflow.png:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --force --configfile ${PWD}/tests/config.yaml --directory ${PWD}/tests --rulegraph | dot -T png > workflow.png

example.png:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --force --configfile ${PWD}/tests/config.yaml --directory ${PWD}/tests --dag | dot -T png > example.png
