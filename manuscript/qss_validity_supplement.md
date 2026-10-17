# Supplementary Note 5. Accuracy of the matched quasi-steady-state reduction

## Question and matched comparison

When does eliminating proneural relaxation preserve growth accurately? We use 50 distinct
full-model-unstable cases selected at evenly spaced ranks of the original sample's sorted
growth rates, plus P2 and P3. This is a coverage design, not a new random sample. In each
family, $(\lambda_S,d_S)\mapsto\alpha(\lambda_S,d_S)$ preserves every equilibrium and every
static loop gain. We evaluate 61 logarithmically spaced values $10^{-2}\leq\alpha\leq10^4$
and all 54 distinct graph modes, including the uniform mode: 3172 family evaluations.
The matched DNA model is unchanged by this compensation.

The relative growth error is $e=|\lambda_{DNA}/\lambda_F-1|$ where $\lambda_F>10^{-8}$;
we retain a scaled absolute error separately. The diagnostic
$q=|s_{\rm lead}|/(\alpha d_S)$ measures the leading pole relative to downstream relaxation.
At that pole, replacing the exact downstream factor
$\alpha d_S/(s+\alpha d_S)$ by one gives the DC-normalised difference
$|1-\alpha d_S/(s+\alpha d_S)|=|s/(s+\alpha d_S)|$. These are mode-specific diagnostics,
not error bounds for an eigenvalue.
The sensitivity of a root of the complete characteristic equation matters as well.

## Results and practical use

For each family and tolerance, we report the first sampled $\alpha$ after which **every later
sampled value** meets the tolerance. These are grid-based empirical thresholds, not certified
continuous-interval bounds. All 52 families reach both tolerances by $10^4$.
The 10% thresholds range from 1 to 3162.28 (median 125.89); the 5% thresholds range from
1.995 to 6309.57 (median 251.19). At P2, the respective thresholds are 199.53 and 398.11;
at P3 they are 398.11 and 794.33. There are no full-faster-than-DNA cases in this sweep.

A small $q$ alone is insufficient for relative growth accuracy: 27 of 677 evaluated cases
with $q\leq0.01$ still exceed 10% growth error. These are dependent parameter-family
evaluations, not independent Bernoulli trials. At original sample 3471 and $\alpha=0.01$,
$q=0.00699905$ and feedback substitution error is 0.0069504,
while the relative growth error is 44.84.

Its gain is close to the stationary boundary. A small change in feedback can produce a large
*relative* change in a near-zero growth rate. The eigenvalue condition recorded in the data
is a descriptive conditioning measure, not a rigorous nonlinear QSS bound.

For quantitative growth prediction, compare the matched full and reduced poles over the
intended parameter range and apply the desired error tolerance there. If only an asymptotic
scale separation is available, state that assumption explicitly; $q\ll1$ does not by itself
guarantee a prescribed relative growth error near onset. No universal cutoff is inferred.

**Table S7. P2/P3 turnover factors required on the evaluated grid.**

| family | growth-error tolerance | $\alpha$ required on grid |
|---|---|---|
| P2 | 0.05 | 398.107 |
| P2 | 0.1 | 199.526 |
| P3 | 0.05 | 794.328 |
| P3 | 0.1 | 398.107 |


The design was recorded in `figures/data/reviewer_extensions/protocol.json` before the
production calculations. All selected cases and counterexamples are retained. Source hashes
and case checksums are verified by `scripts/write_reviewer_reports.py`. The primary package
versions and numerical reproduction contract are given in `docs/REPRODUCIBILITY.md`.


## Reproduction

Run `python scripts/reviewer_linear.py`, followed by
`python scripts/write_reviewer_reports.py` and `python figures/plot_reviewer_extensions.py`.
Full inputs and pole diagnostics are in `qss_families.csv`, `qss_validity.csv` and
`qss_accuracy_thresholds.csv` under `figures/data/reviewer_extensions/`.

![**Figure S7. QSS validity. (A) P2/P3 relative growth error. (B) Error against the pole-to-relaxation ratio, coloured by eigenvalue condition. (C) Grid thresholds for all 52 families, sorted separately by tolerance. (D) Every family; repeated evaluations within a family are dependent.**](../figures/supplementary/FigS7_qss_validity.png){width=100%}
