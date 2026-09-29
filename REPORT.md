# 大麦収量の要因効果分析 / Factorial Analysis of Barley Yield

## Executive summary

The analysis uses all 120 published observations and a pre-specified fixed-effects model, `yield ~ year + variety + site + variety:site`. The model-estimated 1932 - 1931 contrast was -5.31 bushels/acre (95% CI -7.91 to -2.72; p <0.001).

The variety-by-site interaction was not detectable (F(45, 59) = 0.53, p 0.985). Holm-adjusted model contrasts identified 0 of 45 variety pairs and 10 of 15 site pairs at alpha = 0.05. These comparisons use the fitted factorial model; they are not pooled raw pairwise t tests.

AICc ranked `year_by_site` first (AICc 736.95), so the interaction model is retained as the pre-specified inferential model rather than claimed as the uniquely best predictive model. Full rankings and leave-one-out errors are machine-readable.

## Data and design

The public `lattice::barley` data contain mean yield in bushels per acre for 10 varieties at six Minnesota sites in 1931 and 1932. Every year-by-variety-by-site combination appears once. Each published value is already averaged across three randomized blocks. The original plot-level replicates are unavailable, so block variance and the original randomized-block model cannot be recovered.

The official R documentation is the source of record: <https://search.r-project.org/CRAN/refmans/lattice/html/barley.html>. It cites Immer, Hayes, and Powers (1934) and Wright (2013).

## Methods

The primary model includes fixed main effects for year, variety, and site plus variety-by-site interaction. Because the data are a complete balanced factorial, the reported sequential ANOVA terms are orthogonal for this specification. The three-way interaction is not estimable separately from error after block averaging and is not fitted.

Estimated marginal means are linear combinations of the model coefficients, averaged equally over the other design factors. Pairwise contrasts use the model residual variance and degrees of freedom. P-values are Holm-adjusted separately within the year, variety, and site families. Reported 95% contrast intervals are unadjusted model intervals; the standardized intervals divide those limits by the residual standard deviation and are approximate because uncertainty in that standard deviation is not propagated.

Model-selection sensitivity compares nested scientific candidates with AIC, small-sample AICc, BIC, adjusted R-squared, and exact linear-model leave-one-out RMSE. Diagnostics include residual plots, Shapiro-Wilk, a fitted-value Koenker-Breusch-Pagan auxiliary regression, leverage, externally studentized residuals, Cook's distance, and HC3 covariance sensitivity.

## Factorial ANOVA and effect magnitudes

| Term | df | F | p | partial eta2 | omega2 |
| --- | --- | --- | --- | --- | --- |
| year | 1.000 | 16.829 | 0.000 | 0.222 | 0.062 |
| variety | 9.000 | 2.323 | 0.026 | 0.262 | 0.047 |
| site | 5.000 | 26.353 | 0.000 | 0.691 | 0.500 |
| variety:site | 45.000 | 0.532 | 0.985 | 0.289 | 0.000 |

Partial eta-squared describes each term relative to that term plus residual sum of squares. Omega-squared is a less biased variance-share estimate; both should be interpreted descriptively because sites and varieties are fixed and the data are aggregated.

## Year effect and heteroskedasticity sensitivity

| Covariance | Estimate | SE | 95% low | 95% high | p |
| --- | --- | --- | --- | --- | --- |
| OLS covariance | -5.314 | 1.295 | -7.907 | -2.722 | 0.000 |
| HC3 covariance | -5.314 | 1.848 | -9.011 | -1.618 | 0.006 |

The fitted-value heteroskedasticity check gave p 0.128. The Shapiro-Wilk residual check gave p 0.992. The largest Cook's distance was 0.124; 20 observations exceeded the screening threshold 4/n = 0.033. These thresholds are diagnostics, not automatic deletion rules; no observations were removed.

## Model-selection sensitivity

| Model | k | AICc | delta AICc | weight | BIC | LOOCV RMSE |
| --- | --- | --- | --- | --- | --- | --- |
| year_by_site | 22.000 | 736.947 | 0.000 | 1.000 | 787.839 | 5.039 |
| additive | 17.000 | 806.505 | 69.558 | 0.000 | 847.893 | 6.807 |
| year_site | 8.000 | 810.776 | 73.828 | 0.000 | 831.778 | 7.010 |
| variety_site | 16.000 | 825.954 | 89.007 | 0.000 | 865.272 | 7.394 |
| site | 7.000 | 826.500 | 89.553 | 0.000 | 845.012 | 7.490 |
| year_by_variety | 26.000 | 827.418 | 90.471 | 0.000 | 884.796 | 7.262 |

The leading model contains a year-by-site interaction. Its joint interaction test was F(5, 99) = 20.07 (p <0.001). Site-specific model contrasts are below; Holm adjustment covers the six site contrasts.

| Site | 1932 - 1931 | 95% low | 95% high | Holm p | Holm significant |
| --- | --- | --- | --- | --- | --- |
| Grand Rapids | -8.243 | -12.305 | -4.182 | 0.000 | TRUE |
| Duluth | -4.593 | -8.655 | -0.532 | 0.027 | TRUE |
| University Farm | -6.320 | -10.382 | -2.258 | 0.005 | TRUE |
| Morris | 12.227 | 8.165 | 16.288 | 0.000 | TRUE |
| Crookston | -12.480 | -16.542 | -8.418 | 0.000 | TRUE |
| Waseca | -12.477 | -16.538 | -8.415 | 0.000 | TRUE |

AIC-based selection addresses relative expected predictive fit, whereas ANOVA tests address a pre-specified effect. They need not return the same wording of a conclusion. Accordingly, the report keeps the variety-by-site model for the planned interaction test and treats the year-by-site model as strong sensitivity evidence. Results from both models are exported.

## Morris year-label sensitivity

The package documentation notes Cleveland's suggestion that Morris may have its two years switched. The published data give a 1932 - 1931 estimate of -5.31 (95% CI -7.91 to -2.72). Swapping only the Morris year labels gives -9.39 (95% CI -11.02 to -7.76). The interaction p-value changes from 0.985 to 0.140. This scenario is a sensitivity analysis, not a correction to the canonical data.

The year-by-site model has exactly the same likelihood after a within-Morris year-label swap: it can absorb that relabeling by reversing the Morris simple effect. By contrast, the additive model improves sharply after the swap and comes within about 2 AICc units of the year-by-site model. This shows that the published Morris pattern drives much, but not necessarily all, of the apparent year-by-site heterogeneity.

## Conclusions

Year and site differences are large in this aggregated historical experiment, while variety differences are smaller. The model explicitly asks whether varieties behave differently by site; the joint interaction result and all model-based contrasts are reported without replacing the design with independent raw t tests. Conclusions are limited to these observed years, sites, and varieties and should not be presented as a modern agronomic forecast.

## Reproducibility record

Run `make all` from this directory. It rebuilds the data snapshot from `lattice`, runs the analysis, regenerates this report and every table/figure, and executes the tests. The verified runtime is recorded in `results/session_info.txt`; exact R package pins are in `renv.lock`.

