RSCRIPT ?= Rscript

.PHONY: all data analyze test check clean

all: data analyze test

data: data/barley.csv

data/barley.csv: scripts/prepare_data.R
	$(RSCRIPT) --vanilla scripts/prepare_data.R --output=data/barley.csv

analyze: data/barley.csv scripts/run_analysis.R R/analysis_helpers.R
	$(RSCRIPT) --vanilla scripts/run_analysis.R --data=data/barley.csv --results=results --figures=figures

test: analyze
	$(RSCRIPT) --vanilla tests/test_helpers.R
	$(RSCRIPT) --vanilla tests/test_outputs.R

check:
	$(RSCRIPT) --vanilla scripts/check_environment.R

clean:
	rm -f data/barley.csv data/barley.csv.md5 REPORT.md results/*.csv results/*.json results/*.txt figures/*.png
