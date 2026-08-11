# PIEZO1 dynamic structural simulator
#
# One-command reproduction (aim A5): a fresh clone plus `make env` plus
# `make reproduce` rebuilds the entire working state and then checks that every
# number the documentation asserts still comes out of the code.
#
# Everything runs inside the `piezo1` conda environment. CONDA_RUN wraps each
# target so the Makefile works without the environment being pre-activated —
# `make test` from a bare shell does the right thing.

ENV_NAME ?= piezo1
CONDA_RUN := conda run --no-capture-output -n $(ENV_NAME)
PY := $(CONDA_RUN) python

.DEFAULT_GOAL := help
.PHONY: help env lock fetch resources test lint gui reproduce verify quick \
        figures validate clean-derived sizes params audit \
        notebooks coldclone provenance
# `notebooks` and `coldclone` name a directory and a script; without .PHONY
# make sees the target as already built and does nothing.

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
	  | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-16s\033[0m %s\n", $$1, $$2}'

env:  ## Create the conda environment
	bash scripts/create_env.sh $(ENV_NAME)

lock:  ## Regenerate the pinned environment locks
	conda env export -n $(ENV_NAME) --no-builds | grep -v '^prefix:' > environment.lock.yml
	$(CONDA_RUN) python -m pip freeze | grep -viE '^-e |file://' > requirements.lock.txt
	@echo "wrote environment.lock.yml and requirements.lock.txt"

fetch:  ## Download structures, sequences, ligands and the CHAP grid
	$(PY) -m piezo1.io.fetch

resources:  ## Rebuild the curated annotation resources
	$(PY) scripts/build_parameters.py
	$(PY) scripts/build_uniprot_annotations.py
	$(PY) scripts/build_domains.py
	$(PY) scripts/build_functional_residues.py
	$(PY) scripts/build_variants.py
	$(PY) scripts/build_structure_registry.py

test:  ## Run the test suite
	$(PY) -m pytest

lint:  ## Static checks
	$(CONDA_RUN) ruff check piezo1 scripts tests

params:  ## Rebuild the parameter registry and check every citation resolves
	$(PY) scripts/build_parameters.py

audit:  ## Fail if any number in the science modules is unregistered
	$(PY) -m piezo1.parameter_audit

sizes:  ## Fail if any file exceeds the project's 500-line limit
	@find piezo1 scripts tests -name '*.py' -exec wc -l {} + \
	  | awk '$$1 > 500 && $$2 != "total" { print "TOO LONG: " $$2 " (" $$1 " lines)"; bad=1 } \
	         END { if (bad) exit 1; print "all files within 500 lines" }'

gui:  ## Launch the application
	$(PY) -m piezo1

notebooks:  ## Rebuild the example notebooks, running every cell first
	$(PY) scripts/build_notebooks.py

figures:  ## Regenerate documentation figures and screenshots
	$(PY) scripts/make_figures.py
	$(PY) scripts/make_guo2017_figures.py
	$(PY) scripts/screenshot_app.py --structure 8YEZ --analysis

validate:  ## Re-run both pre-registered variant tests
	$(PY) scripts/run_validation.py
	$(PY) scripts/run_validation_round22.py

verify:  ## Check that every documented number still comes out of the code
	$(PY) scripts/reproduce.py --verify

coldclone:  ## Run the suite from an empty clone; a failure is a reproducibility bug
	$(PY) scripts/cold_clone_check.py

provenance:  ## Walk each documented number back to its file, parameters and commit
	$(PY) -m piezo1.analysis.provenance_chain

quick:  ## Fetch, test and verify, skipping the slow steps
	$(PY) scripts/reproduce.py --quick

reproduce:  ## Rebuild everything from scratch and verify the documentation
	$(PY) scripts/reproduce.py

clean-derived:  ## Remove regenerable artefacts (never touches ref/ downloads)
	rm -rf data/derived/reports
	@echo "derived reports removed; re-run 'make validate' and 'make figures'"
