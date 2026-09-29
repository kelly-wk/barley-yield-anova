# Public Data and Analysis Provenance

## Canonical data source

The workflow loads `lattice::barley` from the installed R package. The dataset contains mean yields for 10 barley varieties at six Minnesota sites in 1931 and 1932, giving a balanced 2 × 10 × 6 table of 120 observations.

[`scripts/prepare_data.R`](../scripts/prepare_data.R) converts factor labels to explicit strings, orders all cells deterministically, validates complete factorial coverage, and writes a local CSV plus MD5 record. It does not impute, normalize, filter, or alter yield values. The generated CSV is intentionally not committed to the public repository.

## Analysis decisions

1. Raw pairwise t tests are replaced by contrasts from the fitted factorial model.
2. Model-derived marginal means average equally across the balanced design.
3. Holm adjustment is applied within each factor family.
4. Confidence intervals, standardized effects, diagnostics, influence checks, and HC3 sensitivity are reported.
5. Model-selection sensitivity includes AICc, BIC, and leave-one-out RMSE rather than AIC alone.
6. The Morris-site year-label concern in the official documentation is tested explicitly as a sensitivity analysis.

## Publication boundary

Original analysis code, tests, reports, and figures are covered by the project license. The `lattice` package and its dataset retain their upstream terms; see [`THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md). The dataset contains no personal or sensitive fields.
