# Supplementary Note 10. Boundary proximity and parameter-measure sensitivity

## Small denominators and absolute growth differences

The 5.4-fold median ratio in the primary robustness study could be inflated by tiny full-model
growth rates. We therefore reuse all 40000 frozen primary inputs and examine the 14949
full-model-unstable points with three complementary statistics: the growth ratio,
the absolute gap $\lambda_{DNA}-\lambda_F$, and the bounded symmetric difference
$2(\lambda_{DNA}-\lambda_F)/(\lambda_{DNA}+\lambda_F)$. We also define a relative growth
margin $\eta=\lambda_F/\max_j|s_j(J_{-3})|$ and apply prespecified cuts in $\eta$ and raw
$\lambda_F$. The largest K-pole magnitude supplies a local characteristic rate, not a
new instability criterion.

**Table S12. Growth differences after excluding points near the boundary.**

| quantity cut | minimum value | points | growth ratio median | absolute gap median | symmetric error median | ratio Q25 | ratio Q75 |
|---|---|---|---|---|---|---|---|
| relative_margin | 0 | 14949 | 5.41192 | 0.834161 | 1.37616 | 2.71085 | 12.8954 |
| relative_margin | 0.0001 | 14477 | 5.42728 | 0.871598 | 1.37765 | 2.7476 | 12.8104 |
| relative_margin | 0.001 | 12615 | 5.10149 | 0.929899 | 1.34442 | 2.69919 | 11.4626 |
| relative_margin | 0.01 | 8014 | 4.26796 | 0.920589 | 1.24069 | 2.46579 | 8.5689 |
| relative_margin | 0.05 | 4048 | 3.45319 | 0.823129 | 1.10177 | 2.17015 | 5.86694 |
| lambda_pattern | 0.001 | 14891 | 5.41142 | 0.839149 | 1.37611 | 2.71363 | 12.8791 |
| lambda_pattern | 0.01 | 14421 | 5.35189 | 0.880427 | 1.37027 | 2.72058 | 12.5535 |
| lambda_pattern | 0.1 | 10059 | 4.46656 | 1.19609 | 1.26828 | 2.57695 | 9.1841 |

The ratio median decreases from 5.4119 to 3.4532 when $\eta\geq0.05$ (4048 points).
The median absolute gap in that subset is 0.8231 in the model's inverse dimensionless time.
Thus proximity to onset affects the magnitude but does not by itself explain the sampled
growth discrepancy. A ratio is never used as evidence of experimental fold accuracy.

## Prior sensitivity by importance reweighting

We retain the original support and all stored inputs; no new draws are made. In addition to
the reference independent log-uniform rates and uniform Hill exponents, we apply two
illustrative density ratios: (i) independent triangular Hill densities on [1,6] with mode 3;
(ii) a Gaussian copula with correlation 0.75 between log proneural production and log
proneural loss, preserving their original log-uniform marginals. The latter makes a strong
production/turnover relationship explicit without changing marginal ranges. Neither prior
is estimated from biology.

Each seed's instability fraction is a self-normalised weighted estimate. We average the ten
seed estimates and use a Student-t interval with nine degrees of freedom for Monte Carlo
variation across these designs. This is not uncertainty about biological prevalence.
The pooled effective sample size is $(\sum w)^2/\sum w^2$; we report it because reweighting
reduces information. Growth medians use pooled weights only among full-unstable points.

**Table S13. Prior-sensitive instability and growth summaries.**

| prior | unstable fraction | MC 95% lower | MC 95% upper | ESS | mean weight | weighted growth-ratio median |
|---|---|---|---|---|---|---|
| reference | 37.372% | 37.073% | 37.672% | 40000 | 1 | 5.41192 |
| triangular_Hill | 39.129% | 38.612% | 39.645% | 22497.8 | 1.00144 | 5.24507 |
| correlated_S_turnover | 38.322% | 37.331% | 39.313% | 16699.6 | 1.00256 | 5.96058 |

Instability estimates vary from 37.37% to 39.13% and the weighted median growth ratio from
5.245 to 5.961 for these measures. The bounds study in Note 2 produces a larger fraction
change (to 28.89%). These percentages are properties of specified measures; they must not
be called the fraction of all possible biological systems. Importance sampling does not
test behaviour outside the reference support, and these two alternatives do not exhaust
reasonable priors.

## Reproduction

`scripts/reviewer_linear.py` writes `growth_error_points.csv.gz`,
`growth_boundary_strata.csv`, `prior_replicates.csv` and `prior_sensitivity.csv` from the
checksum-verified primary robustness cases. Use the report and figure commands in Note 5.

The design was recorded in `figures/data/reviewer_extensions/protocol.json` before the
production calculations. All selected cases and counterexamples are retained. Source hashes
and case checksums are verified by `scripts/write_reviewer_reports.py`. The primary package
versions and numerical reproduction contract are given in `docs/REPRODUCIBILITY.md`.

![**Figure S12. Growth and prior sensitivity. (A) Reweighted fractions with ten-design Monte Carlo intervals. (B) Ratio medians and IQRs after relative-margin cuts. (C) Absolute growth gap. (D) Bounded symmetric growth difference; all unstable primary points are retained.**](../figures/supplementary/FigS12_growth_and_prior_sensitivity.png){width=100%}
