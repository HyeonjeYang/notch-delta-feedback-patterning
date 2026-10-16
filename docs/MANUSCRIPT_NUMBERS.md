# Manuscript numbers (revised version)

Each value below is printed by the named script and stored in the named data file under
`figures/data/`. Baseline `p0 = model.params.P0` (literature-inspired set with `k_cis = 0.35`).
The continuation parameter is `lam_N`, with baseline 1.0, P2 = 2.0 and P3 = 5.623, all other
parameters at `p0`. The v1 values are listed where a number changed. The v1 sources are in
`archive/v1_physbio_submission/docs/MANUSCRIPT_NUMBERS.md`.

## Equilibrium and global sample (`scripts/global_sample.py` -> `global_sample_results.csv`)

| quantity | value | v1 |
|---|---|---|
| homogeneous equilibria per parameter set | exactly 1 (Proposition 1) | "multistable branches found: 0" |
| v1 "no positive fixed point" points that have an equilibrium | 733 / 733 | reported as nonexistent |
| max relative equilibrium residual, 4000 points | 3.1e-15 | |
| nonuniform-unstable fraction | **1453/4000 = 36.325%** | 1212/4000 = 30.3% |
| stable / stationary (1 RHP) / complex pair (2 RHP) / real pair (2 RHP) | 2547 / 1441 / 6 / 6 | |
| dominant unstable mode | mu = -3 at all 1453 points | |
| argmax over mu of G(mu) | mu = -3 at all 4000 points (Proposition 3) | |
| A_mu Hurwitz for every mu | all 4000 points (Proposition 2a) | |
| uniform mode unstable | 45 points, all complex pairs, all also K-unstable | |
| G(-3) > 1 agrees with instability | 3988/4000 (99.70%); 12 mismatches, all G < 1 with 2 RHP eigenvalues | 99.79% on 3267 points |
| same, restricted to the v1-valid 3267 points | 3260/3267 (7 mismatches) | 7 mismatches, attributed to "another mode" |
| example mismatch, point 1730 | G(-3) = 0.366; eigenvalues 0.159 +- 0.262i | |
| Nyquist winding = eigenvalue count at mu = -3 | 4000/4000 | |
| DNA reduction: unstable / extra / missed | 1479 / 26 / 0 | |
| DNA/explicit growth ratio where both unstable | median 4.90 (IQR 2.64-12.97); > 1 at 100% of points | |
| DN reduction: unstable / extra / missed | 1500 / 47 / 0 | |
| DN/explicit growth ratio | median 216 (IQR 41-1114) | |

The 12 two-eigenvalue points and their eigenvalues are listed in `two_eigenvalue_cases.csv`
(`scripts/nyquist_examples.py`).

## Mode dominance search (`scripts/mode_dominance_search.py` -> `mode_dominance_search.csv`)

| box | samples | nonuniform-unstable | another mode dominant | another mode unstable, K stable |
|---|---|---|---|---|
| v1 box | 20 000 | 7292 | 0 | 0 |
| 10x wider box | 20 000 | 5755 | 0 | 0 |

## lam_N continuation (`scripts/linear_analysis.py`, `scripts/verify_critical_point.py`)

| quantity | value |
|---|---|
| critical lam_N (brentq on Re s = 0, bracket [1.2, 1.4], xtol 1e-13) | 1.304815343 |
| critical lam_N (separate brentq on G(-3) = 1) | 1.304815343, difference 1.1e-15 |
| G(-3) at the pole crossing | 1.000000000000 |
| G(+6) at the pole crossing | -0.0628 |
| max Re of the mu = 6 pole over lam_N in [0.5, 6.5] | -0.100000 |
| G(-3) at p0 / P2 / P3 | 0.4558 / 1.8088 / 3.3242 |
| Schur vs direct characteristic polynomial, max coefficient difference (9 cases) | 1.14e-13 |
| Schur vs direct, max root mismatch | 7.1e-14 |
| det J_mu vs det(J_loc) det(I - mu K), max relative difference | 9.5e-16 |
| dispersion formula vs 576 x 576 adjacency spectrum | 2.842e-14 |

## Kinetic control (`scripts/kinetic_control.py`)

| quantity | value |
|---|---|
| d lambda_pattern / d ln p at p0, top six | n_S 0.232, d_A -0.179, lam_N 0.177, K_d -0.118, lam_S -0.060, d_S +0.044 |
| d ln G(-3) / d ln p at p0, top six | n_S 6.91, d_A -5.32, lam_N 5.32, K_d -3.52, lam_S -1.80, d_S +1.80 |
| parameters whose 1-D continuation over [0.1, 10] x p0 crosses lambda_pattern = 0 | K_d, d_A, d_S, lam_S, lam_N, n_S |
| 2-D maps: points / unstable / G > 1 agreement / missing equilibria | 1875 / 1087 / 100% / 0 (v1: 115 missing) |

## Proneural-timescale family (`scripts/timescale_family.py` -> `timescale_family.csv`)

| quantity | value |
|---|---|
| critical lam_N for alpha in {1e-2, 1e-1, 1, 1e1, 1e2, 1e4} | 1.304815342591 in every case |
| spread of G(-3) at P2 across alpha | 3.8e-15 |
| P2 growth rate across alpha in [1e-2, 1e4] | 8.1e-4 to 1.594 |
| P2 growth rate: explicit / slow-S estimate d_S(G-1) / DNA / DN | 0.103067 / 0.0809 / 1.597479 / 11.277 |
| P3 growth rate: explicit / DNA | 0.216300 / 2.799016 |

## Nonlinear growth validation (`scripts/growth_validation.py` -> `growth_validation.csv`)

The initial state is `x* + eps v cos(phi)`, where `v` is the exact eigenvector and `eps` gives
a relative Delta amplitude of 1e-6. The window runs from t = 0 until the K-mode amplitude of D reaches 1e-3 D*.

| run | predicted | fitted | window | relative error | slope on [0, 8] |
|---|---|---|---|---|---|
| explicit, P2 | 0.103067 | 0.103077 | [0, 67.25] | 1.0e-4 | 0.1031 |
| explicit, P3 | 0.216300 | 0.216305 | [0, 32.0] | 2.4e-5 | 0.2162 |
| DNA reduction, P2 | 1.597479 | 1.597504 | [0, 4.5] | 1.6e-5 | (saturates) |
| explicit, P2, v1 initial state | 0.103067 | 0.113653 | [0, 64.25] | 0.10 | **0.443** |
| explicit, P3, v1 initial state | 0.216300 | 0.257989 | [0, 30.0] | 0.19 | **0.568** |

| quantity | value |
|---|---|
| time to relative K-mode amplitude 0.1 at P2: explicit / DNA | 110.75 / 7.25 |
| final P2 Delta field, max abs(D_explicit - D_DNA) | 7.4e-11 |

## Removed from the manuscript

Numbers that appeared in v1 but no longer appear in the revision. They are still
reproducible from the archive.

* The G_eff heuristic (93.88% best-threshold accuracy).
* The nonlinear saturation statistics (q_amp, Moran's I, solver comparison).
* The exploratory stochastic check (`stochastic_check/`).
