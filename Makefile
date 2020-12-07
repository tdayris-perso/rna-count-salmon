SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eic
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
CONDA_ACTIVATE   = source "$$(conda info --base)/etc/profile.d/conda.sh" ; conda activate ; conda activate

# Paths
TEST_CONFIG      = scripts/prepare_config.py
TEST_DESIGN      = scripts/prepare_design.py
TEST_COMMON      = scripts/common_script_rna_count_salmon.py
TEST_CALLER      = rna-count-salmon.py
SNAKE_FILE       = Snakefile
ENV_YAML         = envs/workflow.yaml
ENV_FLAMINGO     = envs/workflow_flamingo.yaml
ENV_LOCAL        = envs/workflow_local.yaml
TRANSCRIPT_PATH  = tests/genome/transcriptome.fasta
GTF_PATH         = tests/genome/annot.gtf
READS_PATH       = tests/reads/

CONFIG_CALL      = ${PYTHON} rna-count-salmon.py config
DESIGN_CALL      = ${PYTHON} rna-count-salmon.py design
SNAKEF_CALL      = ${PYTHON} rna-count-salmon.py snakemake
REPORT_CALL      = ${PYTHON} rna-count-salmon.py report

# Arguments
PYTEST_ARGS      = -vv
ENV_NAME         = rna-count-salmon
SNAKE_THREADS    = 1
SAINDEX_ARGS     = ' --kmerLen 5 '
SAQUANT_ARGS     = ' --noBiasLengthThreshold --noFragLengthDist --noSingleFragProb --minAssignedFrags 1 --noEffectiveLengthCorrection --noLengthCorrection '

# Recipes
default: quantification-report.html


conda-install-flamingo:
	${CONDA_ACTIVATE} base && \
	${CONDA} env create --file ${ENV_FLAMINGO} --force && \
	${CONDA} activate ${ENV_NAME}
.PHONY: conda-tests


conda-install-local:
	${CONDA_ACTIVATE} base && \
	${CONDA} env create --file ${ENV_LOCAL} --force && \
	${CONDA} activate ${ENV_NAME}
.PHONY: conda-tests


config.yaml:
	${CONDA_ACTIVATE} ${ENV_NAME} && ${CONFIG_CALL} ${FASTA} ${GTF} --debug


design.tsv:
	${CONDA_ACTIVATE} ${ENV_NAME} && ${DESIGN_CALL} . --recursive --debug


quantification-report.html: config.yaml design.tsv
	${CONDA_ACTIVATE} ${ENV_NAME} && ${SNAKEF_CALL} && ${REPORT_CALL}


### UNIT TESTS ###
# Environment building through conda
conda-tests:
	${CONDA_ACTIVATE} base && \
	${CONDA} env create --file ${ENV_YAML} --force && \
	${CONDA} activate ${ENV_NAME}
.PHONY: conda-tests


# Running all unit-tests (one for each python scripts)
all-unit-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} ${PYTEST_ARGS} ${TEST_CONFIG} ${TEST_DESIGN} ${TEST_COMMON} ${TEST_CALLER}
.PHONY: all-unit-tests


# Running all unit test (on prepare_config.py only)
config-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} ${PYTEST_ARGS} ${TEST_CONFIG} && \
	${CONFIG_CALL} ${TRANSCRIPT_PATH} ${GTF_PATH} --salmon-index-extra ${SAINDEX_ARGS} --salmon-quant-extra ${SAQUANT_ARGS} --aggregate --libType ISF --workdir tests --debug -o tests/config.yaml
.PHONY: config-tests


# Running all unit test (on prepare_design.py only)
design-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} ${PYTEST_ARGS} ${TEST_DESIGN} && \
	${DESIGN_CALL} ${READS_PATH} -o tests/design.tsv --debug
.PHONY: design-tests


# Running all unit test (on common_scrit_rna_count_salmon.py only)
common-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} ${PYTEST_ARGS} ${TEST_COMMON}
.PHONY: common-tests


### Continuous Integration Tests ###
# Running snakemake on test datasets
test-cli-wrapper-report.html:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	declare -x SNAKEMAKE_OUTPUT_CACHE="${PWD}/tests/snakemake/cache" && \
	declare -x SNAKEFILE="${PWD}/Snakefile" && \
	declare -x PROFILE="${PWD}/.igr/profile/local" && \
	declare -x PREPARE_CONFIG="${PWD}/scripts/prepare_config.py" && \
	declare -x PREPARE_DESIGN="${PWD}/scripts/prepare_design.py" && \
	declare -x FASTA="${PWD}/tests/genome/transcriptome.fasta" && \
	declare -x GTF="${PWD}/tests/genome/annot.gtf" && \
	declare -x RNA_COUNT_LAUNCHER="${PWD}/rna-count-salmon.py" && \
	export SNAKEMAKE_OUTPUT_CACHE SNAKEFILE PROFILE PREPARE_CONFIG PREPARE_DESIGN FASTA GTF RNA_COUNT_LAUNCHER && \
	${DESIGN_CALL} ${READS_PATH} -o ${PWD}/tests/design.tsv --debug && \
	${CONFIG_CALL} ${TRANSCRIPT_PATH} ${GTF_PATH} \
	--salmon-index-extra ${SAINDEX_ARGS} \
	--salmon-quant-extra ${SAQUANT_ARGS} \
	--aggregate --libType ISF --workdir ${PWD}/tests \
	--design ${PWD}/tests/design.tsv --threads ${SNAKE_THREADS} --debug && \
	${SNAKEF_CALL} --snakemake-args "--configfile tests/config.yaml" && \
	${REPORT_CALL} --snakemake-args "--configfile tests/config.yaml"


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
