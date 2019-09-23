### Variables ###
# Tools
PYTEST           = pytest
BASH             = /bin/bash
CONDA            = conda
PYTHON           = python3.7
SNAKEMAKE        = snakemake

# Paths
TEST_CONFIG      = scripts/prepare_config.py
TEST_DESIGN      = scripts/prepare_design.py
TEST_AGGREGATION = scripts/aggregate_samples.py
SNAKE_FILE       = Snakefile
ENV_YAML         = envs/workflow.yaml
TRANSCRIPT_PATH  = genomes/transcriptome.fasta
READS_PATH       = reads/

# Arguments
ENV_NAME         = rna-count-salmon
SNAKE_THREADS    = 1
SAINDEX_ARGS     = ' --kmerLen 5 '
SAQUANT_ARGS     = ' --noBiasLengthThreshold --minAssignedFrags 1 --noEffectiveLengthCorrection --noLengthCorrection --fasterMapping --noFragLengthDist --allowDovetail --numPreAuxModelSamples 0 --numAuxModelSamples 0 '

# Recipes
default: all-unit-tests

### UNIT TESTS ###
# Running all tests
all-unit-tests: SHELL:=$(BASH)
all-unit-tests:
	"$([ -n "$CONDA_EXE" ] && echo "$CONDA_EXE" || which conda)" info
	$(CONDA) activate $(ENV_NAME)
	$(PYTEST) -v $(TEST_CONFIG) $(TEST_DESIGN) $(TEST_AGGREGATION)

# Running tests on configuration only
config-tests: SHELL:=$(BASH) -i
config-tests:
	$(CONDA) activate $(ENV_NAME)
	$(PYTEST) -v $(TEST_CONFIG)
	$(PYTHON) $(TEST_CONFIG) $(TRANSCRIPT_PATH) --salmon-index-extra $(SAINDEX_ARGS) --salmon-quant-extra $(SAQUANT_ARGS) --aggregate --libType "ISF" --workdir tests -o tests/config.yaml

# Running tests on design only
design-tests: SHELL:=$(BASH) -i
design-tests:
	$(CONDA) activate $(ENV_NAME)
	$(PYTEST) -v $(TEST_DESIGN)
	$(PYTHON) $(TEST_DESIGN) $(READS_PATH) -o tests/design.tsv

# Running tests on aggregation only
aggregation-tests: SHELL:=$(BASH) -i
aggregation-tests:
	$(CONDA) activate $(ENV_NAME)
	$(PYTEST) -v $(TEST_AGGREGATION)

### Continuous Integration Tests ###
# Running snakemake on test datasets
ci-tests: SHELL:=$(BASH) -i
ci-tests:
	$(CONDA) activate $(ENV_NAME)
	$(PYTHON) $(TEST_DESIGN) $(READS_PATH) -o tests/design.tsv
	$(PYTHON) $(TEST_CONFIG) $(TRANSCRIPT_PATH) --salmon-index-extra $(SAINDEX_ARGS) --salmon-quant-extra $(SAQUANT_ARGS) --aggregate --libType "ISF" --workdir tests -o tests/config.yaml -d tests/design.tsv
	$(SNAKEMAKE) -s $(SNAKE_FILE) --use-conda -j $(SNAKE_THREADS) --force
	$(SNAKEMAKE) -s $(SNAKE_FILE) --use-conda -j $(SNAKE_THREADS) --report

# Environment building through conda
conda-tests: SHELL:=$(BASH) -i
conda-tests:
	$(CONDA) env create --file $(ENV_YAML) --force
	$(CONDA) activate $(ENV_NAME)

### Cleaning tests results ###
clean:
	rm -rf .snakemake salmon_index/ qc/ raw_data/ pseudo_mapping/ logs/ genome aggregated_salmon_counts/ report.html config.yaml design.tsv
