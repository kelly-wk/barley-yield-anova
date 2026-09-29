# Data card: Minnesota barley yield

## Identity

- R object: `lattice::barley`
- Local snapshot: `data/barley.csv`
- Shape: 120 rows, four columns
- Unit of observation: one year × variety × site cell
- Outcome: `yield`, bushels per acre
- Predictors: two years, 10 varieties, six sites
- Missing values: none

The official R documentation describes the values as yields from 1/40-acre plots. Ten varieties were grown in three randomized blocks at six Minnesota stations in each of two years. The published `barley` values are averages across the three blocks, not the individual plot measurements.

## Provenance and citation

Canonical documentation:

- R/CRAN help: <https://search.r-project.org/CRAN/refmans/lattice/html/barley.html>
- CRAN package record: <https://CRAN.R-project.org/package=lattice>

Primary source cited by the R documentation:

> Immer, R. F., Hayes, H. K., and Powers, LeRoy (1934). Statistical Determination of Barley Varietal Adaptation. *Journal of the American Society of Agronomy*, 26, 403–419.

The R documentation also cites Kevin Wright (2013), “Revisiting Immer's Barley Data,” *The American Statistician*, 67(3), 129–133.

## Preparation

`scripts/prepare_data.R` loads the object from the installed `lattice` package, converts the factors to explicit strings, orders the 120 cells deterministically, validates the complete 2 × 10 × 6 factorial coverage, and writes the CSV plus an MD5 integrity record. It does not impute, normalize, filter, or alter yield values.

## Known issue and sensitivity analysis

The official documentation notes inconsistent year labels in the historical source and reports Cleveland's suggestion that the Morris-site years may have been switched. The primary analysis uses the R package exactly as published. `scripts/run_analysis.R` separately swaps only the two Morris year labels and reports how conclusions change; it never silently “corrects” the canonical data.

## Limitations

- The three randomized blocks have already been averaged, so plot-level variance and block effects cannot be estimated.
- Years used different land, according to the package documentation. A two-year fixed contrast should not be generalized as a climate trend.
- Sites and varieties are treated as fixed levels. The study does not estimate a population distribution of sites or genotypes.
- With one aggregated observation per full factorial cell, the three-way interaction cannot be separated from residual error.
- Historical transcription concerns remain, especially at Morris.

## License and privacy

The snapshot is generated locally from the `lattice` package, licensed GPL-2 or GPL-3, and is not committed to this public repository. Consult the current CRAN record before redistributing it independently. The public experimental data contain no personal or sensitive fields. See `THIRD_PARTY_NOTICES.md` for attribution and the upstream boundary.
