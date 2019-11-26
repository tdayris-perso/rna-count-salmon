### Variables ###
# Tools
PYTEST           = pytest
BASH             = bash
CONDA            = conda
PYTHON           = python3.7
SNAKEMAKE        = snakemake
CONDA_ACTIVATE   = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

# Paths
TEST_CONFIG      = scripts/prepare_config.py
TEST_DESIGN      = scripts/prepare_design.py
TEST_AGGREGATION = scripts/aggregate_samples.py
SNAKE_FILE       = Snakefile
ENV_YAML         = envs/workflow.yaml
TRANSCRIPT_PATH  = genomes/transcriptome.fasta
READS_PATH       = tests/reads/

# Arguments
ENV_NAME         = rna-count-salmon
SNAKE_THREADS    = 1
SAINDEX_ARGS     = ' --kmerLen 5 '
SAQUANT_ARGS     = ' --noBiasLengthThreshold --minAssignedFrags 1 --noEffectiveLengthCorrection --noLengthCorrection --fasterMapping --noFragLengthDist --allowDovetail --numPreAuxModelSamples 0 --numAuxModelSamples 0 '

# Recipes
default: all-unit-tests


# Environment building through conda
conda-tests: SHELL:=$(BASH) -i
conda-tests:
	$(CONDA_ACTIVATE) base && \
	$(CONDA) env create --file $(ENV_YAML) --force && \
	$(CONDA) activate $(ENV_NAME)

### UNIT TESTS ###
# Running all tests
all-unit-tests: SHELL:=$(BASH) -i
all-unit-tests:
	$(CONDA_ACTIVATE) $(ENV_NAME) && \
	$(PYTEST) -v $(TEST_CONFIG) $(TEST_DESIGN) $(TEST_AGGREGATION)

# Running tests on configuration only
config-tests: SHELL:=$(BASH) -i
config-tests:
	$(CONDA_ACTIVATE) $(ENV_NAME) && \
	$(PYTEST) -v $(TEST_CONFIG) && \
	$(PYTHON) $(TEST_CONFIG) $(TRANSCRIPT_PATH) --salmon-index-extra $(SAINDEX_ARGS) --salmon-quant-extra $(SAQUANT_ARGS) --aggregate --libType "ISF" --workdir tests --debug

# Running tests on design only
design-tests: SHELL:=$(BASH) -i
design-tests:
	$(CONDA_ACTIVATE) $(ENV_NAME) && \
	$(PYTEST) -v $(TEST_DESIGN) && \
	$(PYTHON) $(TEST_DESIGN) $(READS_PATH) -o tests/design.tsv --debug

# Running tests on aggregation only
aggregation-tests: SHELL:=$(BASH) -i
aggregation-tests:
	$(CONDA_ACTIVATE) $(ENV_NAME) && \
	$(PYTEST) -v $(TEST_AGGREGATION)

### Continuous Integration Tests ###
# Running snakemake on test datasets
ci-tests: SHELL:=$(BASH) -i
ci-tests:
	ls && pwd && \
	$(CONDA_ACTIVATE) $(ENV_NAME) && \
	$(PYTHON) $(TEST_DESIGN) $(READS_PATH) -o ${PWD}/tests/design.tsv --debug && \
	$(PYTHON) $(TEST_CONFIG) $(TRANSCRIPT_PATH) --salmon-index-extra $(SAINDEX_ARGS) --salmon-quant-extra $(SAQUANT_ARGS) --aggregate --libType "ISF" --workdir ${PWD}/tests --design ${PWD}/tests/design.tsv --threads $(SNAKE_THREADS) --debug && \
	$(SNAKEMAKE) -s $(SNAKE_FILE) --use-conda -j $(SNAKE_THREADS) --forceall --printshellcmds --reason --directory ${PWD}/tests && \
	$(SNAKEMAKE) -s $(SNAKE_FILE) --use-conda -j $(SNAKE_THREADS) --directory ${PWD}/tests --report


# Running within Singularity
singularity-tests: SHELL:=$(BASH) -i
singularity-tests:
	$(CONDA_ACTIVATE) $(ENV_NAME) && \
	$(PYTHON) $(TEST_DESIGN) $(READS_PATH) -o ${PWD}/tests/design.tsv --debug && \
	$(PYTHON) $(TEST_CONFIG) $(TRANSCRIPT_PATH) --salmon-index-extra $(SAINDEX_ARGS) --salmon-quant-extra $(SAQUANT_ARGS) --aggregate --libType "ISF" --workdir ${PWD}/tests --design ${PWD}/tests/design.tsv --threads $(SNAKE_THREADS) --debug && \
	$(SNAKEMAKE) -s $(SNAKE_FILE) --use-conda -j $(SNAKE_THREADS) --forceall --printshellcmds --reason --directory ${PWD}/tests --use-singularity && \
	$(SNAKEMAKE) -s $(SNAKE_FILE) --use-conda -j $(SNAKE_THREADS) --directory ${PWD}/tests --report


# Cleaning
clean: SHELL:=$(BASH) -i
clean:
	$(CONDA_ACTIVATE) $(ENV_NAME) && \
	$(SNAKEMAKE) -s $(SNAKE_FILE) --use-conda -j $(SNAKE_THREADS) --force --configfile ${PWD}/tests/config.yaml --use-singularity --directory ${PWD}/tests --delete-all-output
