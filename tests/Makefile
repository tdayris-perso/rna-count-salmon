### Variables ###
# Pytest
PYTEST           = pytest
TEST_CONFIG      = ../scripts/prepare_config.py
TEST_DESIGN      = ../scripts/prepare_design.py
TEST_AGGREGATION = ../scripts/aggregate_samples.py

# Snakemake
BASH             = bash
TEST_PIPELINE    = test.sh


### UNIT TESTS ###
# Running all tests
all-unit-tests:
	$(PYTEST) -v $(TEST_CONFIG) $(TEST_DESIGN) $(TEST_AGGREGATION)

# Running tests on configuration only
config-tests:
	$(PYTEST) -v $(TEST_CONFIG)

# Running tests on design only
design-tests:
	$(PYTEST) -v $(TEST_DESIGN)

# Running tests on aggregation only
alggregation-tests:
	$(PYTEST) -v $(TEST_AGGREGATION)

### Continuous Integration Tests ###
# Running snakemake on test datasets
ci-tests:
	$(BASH) $(TEST_PIPELINE)

### Cleaning tests results ###
clean:
	rm -rf salmon_index/ qc/ raw_data/ pseudo_mapping/ logs/ genome aggregated_salmon_counts/ report.html config.yaml design.tsv
