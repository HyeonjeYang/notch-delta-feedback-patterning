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
| max relative equilibrium residual, 4000 points | 2.8e-15 | |
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
| 10x wider box | 20 000 | 5756 | 0 | 0 |

## lam_N continuation (`scripts/linear_analysis.py`, `scripts/verify_critical_point.py`)

| quantity | value |
|---|---|
| critical lam_N (brentq on Re s = 0, bracket [1.2, 1.4], xtol 1e-13) | 1.304815343 |
| critical lam_N (separate brentq on G(-3) = 1) | 1.304815343, difference 6.7e-15 |
| G(-3) at the pole crossing | 1.000000000000 |
| G(+6) at the pole crossing | -0.0628 |
| max Re of the mu = 6 pole over lam_N in [0.5, 6.5] | -0.100000 |
| G(-3) at p0 / P2 / P3 | 0.4558 / 1.8088 / 3.3242 |
| Schur vs direct characteristic polynomial, max coefficient difference (9 cases) | 1.71e-13 |
| Schur vs direct, max root mismatch | 9.95e-14 |
| det J_mu vs det(J_loc) det(I - mu K), max relative difference | 1.18e-15 |
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
| spread of G(-3) at P2 across alpha | 2.7e-15 |
| P2 growth rate across alpha in [1e-2, 1e4] | 8.1e-4 to 1.594 |
| P2 growth rate: explicit / slow-S estimate d_S(G-1) / DNA / DN | 0.103067 / 0.0809 / 1.597479 / 11.277 |
| P3 growth rate: explicit / DNA | 0.216300 / 2.799016 |

## Nonlinear growth validation (`scripts/growth_validation.py` -> `growth_validation.csv`)

The initial state is `x* + eps v cos(phi)`, where `v` is the exact eigenvector and `eps` gives
a relative Delta amplitude of 1e-6. The window runs from t = 0 until the K-mode amplitude of D reaches 1e-3 D*.

| run | predicted | fitted | window | relative error | slope on [0, 8] |
|---|---|---|---|---|---|
| explicit, P2 | 0.103067 | 0.103077 | [0, 67.25] | 9.9e-5 | 0.1030 |
| explicit, P3 | 0.216300 | 0.216305 | [0, 32.0] | 2.1e-5 | 0.2162 |
| DNA reduction, P2 | 1.597479 | 1.597498 | [0, 4.5] | 1.2e-5 | (saturates) |
| explicit, P2, v1 initial state | 0.103067 | 0.113653 | [0, 64.25] | 0.10 | **0.443** |
| explicit, P3, v1 initial state | 0.216300 | 0.257290 | [0, 30.0] | 0.19 | **0.562** |

| quantity | value |
|---|---|
| time to relative K-mode amplitude 0.1 at P2: explicit / DNA | 110.75 / 7.25 |
| final P2 Delta field, max abs(D_explicit - D_DNA) | 8.6e-11 |

## Removed from the manuscript

Numbers that appeared in v1 but no longer appear in the revision. They are still
reproducible from the archive.

* The G_eff heuristic (93.88% best-threshold accuracy).
* The nonlinear saturation statistics (q_amp, Moran's I, solver comparison).
* The exploratory stochastic check (`stochastic_check/`).

## Revision audit on 17 October 2026

`python scripts/audit_revision.py` recomputes every 4000-point full and reduced mode scan,
static gain and Nyquist count. `revision_audit.json` records its summary and
`audit_mode_scan.csv` retains all 54 distinct graph eigenvalues for every point.

* 1451 points also have unstable subdominant nonuniform modes. At P2, 21 distinct
  nonuniform graph eigenvalues are unstable. The earlier sentence that only the uniform
  mode could also be unstable was false; K-mode dominance is retained (Fig. S3).
* The complete global results are reproduced within 4e-12 absolute error in the stored
  statistics. The 36.325% fraction, 12 gain mismatches, and 26/47 reduced-only points persist.
* Drift Hill functions now match the stated equations and analytic derivatives without
  the historical S floor of 1e-8 or A/K cap of 1e5. No global-sample equilibrium reached those
  cutoffs. In the wider random search, removing the floor changes the unstable count by one:
  5756 instead of 5755, for 13048 total across both boxes. Dominance conclusions are unchanged.
* The category string `stationary` in the CSV is a legacy label for an odd RHP count; a
  parameter sample does not determine its bifurcation history. The six real-pair gain
  mismatches are not evidence by themselves of a Hopf crossing.
* The amplitude rule for growth fitting is prespecified, while each endpoint is determined
  by its trajectory. Equality of nonlinear equilibrium sets is not equality of their stability
  or basins; the common P2 final pattern is limited to the specified simulations.
* The graph determinant identity requires nonsingular J_loc. It does not verify the extra
  hypotheses of the [Arcak input-output theorems](https://arxiv.org/pdf/1201.1477) or the
  [Ferreira and Arcak patterned-state analysis](https://epubs.siam.org/doi/10.1137/130910142).
* The Markdown manuscript and response are updated. Existing DOCX files predate this audit;
  rebuild and visually verify them before submission. Figure 1A remains intentionally blank
  for the author's circuit schematic.

## Reproducibility follow-up

The numerical environment and full sampling protocols are documented in
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). The requirements were corrected to the observed
and tested environment. Both studies regenerate their frozen inputs; the additional search
now retains all 40000 parameter vectors in `mode_search_samples.csv.gz`.

`numerical_robustness.json` records the targeted tolerance checks: 1453 unstable points at
absolute thresholds 0 through 1e-6; Nyquist counts invariant over 2000/4000/8000 grid sizes
at 28 selected points; maximum tighter-run growth-rate error about 1.02e-4 relative to the
linear prediction. The full and reduced fitted slopes need not reproduce in their last
printed digits across environments. Those numerical checks alone do not test independent-seed
robustness; the additional study below now does so.

## Independent sampling robustness

Source: `scripts/robustness_study.py` and the frozen protocol/results under
`figures/data/robustness/`. `scripts/write_robustness_report.py` generates the full seed-level
table and Supplementary Note 2 in `manuscript/robustness_supplement.md`.

| quantity | result |
|---|---|
| new points / designs | 115000 / 30 |
| primary original-box study | 10 independent seeds x 4000 points |
| primary nonuniform-unstable count | 14949 / 40000 |
| mean primary fraction / approximate 95% MC interval | 37.37% / 37.07–37.67% |
| primary individual fraction range | 36.650–37.925% |
| primary pooled DNA/full growth-ratio median | 5.411923 |
| primary DNA replicate-median range | 5.250–5.609 |
| primary DN replicate-median range | 201.9–240.3 |
| primary DC misses / DNA-only / DN-only | 150 / 216 / 445 |
| original-box size checks | 5 seeds each at n = 1000, 2000, 8000; separate non-nested LHS |
| size-check mean fractions at n = 1000, 2000, 8000 | 37.10%, 37.69%, 37.39% |
| wider-box study | 5 independent seeds x 4000 points; rate/constant bounds [0.01,100] x p0 |
| wider-box mean fraction / approximate 95% MC interval | 28.89% / 28.58–29.20% |
| wider-box DNA / DN replicate-median ranges | 2.99–3.64 / 1493–2019 |
| all new unstable evaluations / DC misses | 41307 / 392 |
| reduction-missed full instabilities / reduced-growth reversals | 0 / 0 in analysed samples |
| faster competing nonuniform modes | 0 in analysed samples |
| targeted Nyquist checks / disagreements | 842 / 0; all gain misses plus closest-to-unit-gain points |
| 80-digit checked cases | 11; stability signs and K-mode RHP counts agree |
| near-marginal full-model point | wider seed 2026101803, ID 2229, growth 4.772659e-9 |
| max equilibrium production/loss backward error | 8.2e-15 |

The original 36.325% and median ratio 4.90 remain reference-sample results. The primary
replicated estimate supersedes neither those stored inputs nor their reproducible counts;
it estimates the chosen parameter measure with less sampling uncertainty. Interval units are
independent LHS replicates. Do not pool sizes sharing seeds into one iid confidence interval.
Changing bounds changes the integration measure, so neither fraction nor median is universal.

## Initial conditions and mechanistic predictions

Sources: `scripts/initial_condition_robustness.py`, `scripts/write_initial_condition_report.py`,
the frozen protocol, per-case arrays and `runs.csv`, `pairs.csv`, `summary.json` under
`figures/data/initial_conditions/`. Supplementary Note 3 gives all six seed-group summaries.

| quantity | result |
|---|---|
| production trajectories / main pairs / independent-S controls | 112 / 48 / 16 |
| operating points / seeds / amplitudes | P2 and P3 / 2026101901–2026101908 / 1e-6, 1e-4, 1e-2 |
| fixed horizons | P2: 600; P3: 300 |
| full later at CV(D)=0.1 | 48 / 48 paired cases |
| full/DNA onset-time ratios, P2 | 10.95994–15.10064 |
| full/DNA onset-time ratios, P3 | 12.33987–12.89529 |
| full later at secondary CV(D)=0.5 | 48 / 48 pairs; ratios P2: 11.86819–14.84389, P3: 12.08007–12.77015 |
| translation-aligned endpoint correlations, P2 | 0.58435–0.83779 |
| translation-aligned endpoint correlations, P3 | 0.86596–0.95273 |
| aligned relative endpoint L2 distances, all main pairs | 0.25884–0.77774 |
| endpoint exact-K power fraction, 96 main-pair endpoints | 0.002258–0.423762 |
| endpoint exact-K power fraction, all 112 including free-S controls | 0.002258–0.801365 |
| stationary candidates (scaled drift < 1e-7) | 112 / 112; not a nonlinear stability proof |
| min recorded concentration | 2.19147e-10 |
| largest valid-window relative L2 Delta linearisation error | 0.00180027 (0.180027%) |
| full onset free-S / QSS-S, P2 | 1.10087–1.22166 (median 1.17239) |
| full onset free-S / QSS-S, P3 | 1.09969–1.19326 (median 1.15587) |
| tighter full-horizon BDF cases / independent early DOP853 cases | 4 / 4 |
| largest onset difference in solver checks | 0.005884 time units |

`scripts/mechanistic_predictions.py` writes `figures/data/mechanistic_predictions.csv`.
At compensated alpha=0.25, 1, 4, P2 K-mode growth rates are 0.02203829, 0.10306672 and
0.38307960; critical lam_N remains 1.304815343. The clamped-A proneural relaxation
times are 40, 10 and 2.5; half-power angular frequencies are 0.025, 0.1 and 0.4.
All are dimensionless, model-based values; no tissue-specific prediction accuracy is implied.

## Time-discretisation audit

Sources: `scripts/timestep_audit.py`, `scripts/write_timestep_report.py`, the frozen
protocol and the per-case records under `figures/data/timestep_audit/`.

| quantity | result |
|---|---|
| representative lattice initial states / new full-horizon trajectories | 14 / 70 |
| internal step caps at original BDF tolerance | 0.5, 0.25, 0.125 |
| tighter BDF / independent Radau | rtol=1e-11, atol=1e-14; caps 0.125 / 0.25 |
| largest relative CV=0.1 or 0.5 event-time change | 0.0001631001 (0.01631001%) |
| largest raw endpoint relative L2 difference, all species | 1.113310e-10 |
| largest absolute endpoint K-power fraction change | 4.022338e-13 |
| largest new endpoint scaled drift | 5.481252e-9 |
| minimum accepted adaptive concentration | 2.456092e-10 |
| mode/model/fixed-step spectral comparisons | 1296 |
| coarse RK4 spectral count failures | 216, all at P0 with h=0.04 or 0.08 |
| spectral count failures at h=0.02 / 0.01 | 0 / 0 (648 finer-step comparisons) |
| nonlinear fixed RK4 lattice growth runs / independent Radau references | 16 / 4 |
| largest fixed h=0.01 Delta perturbation error vs Radau | below 3.8e-8 |
| largest fixed-step fitted-growth relative error vs eigenvalue | 1.004506e-4 |
| largest relative output-cadence slope change | 2.104783e-5 (0.002104783%) |
| reintegrated main pairs with full later onset | 30 / 30 |
| coarse fixed field accuracy failure | P3 DNA at h=0.08: 1.253026e-4 exceeds 1e-4 target |

These fixed-step failures are diagnostic controls, not the integration settings used in
the production adaptive runs. The original output interval does not specify the adaptive
internal dt. Supplementary Note 4 gives all nonlinear fixed-step growth, output-cadence
fits and assessment flags, retaining any coarse-step accuracy failures.

Post-audit window diagnostic: reusing the same high-accuracy amplitudes at theoretical
upper cutoffs 1e-3, 1e-4 and 1e-5 gives maximum relative slope discrepancies 1.004506e-4,
1.843778e-5 and 4.045175e-6. This extra diagnostic is labelled post-audit and does not
replace the frozen fitting/assessment protocol or any reference-manuscript fit.


## Additional matched controls: main Figure 6 and Notes 5-11

All files below are under `figures/data/reviewer_extensions/`. The frozen `protocol.json`
fixes designs before results; `linear_summary.json`, `geometry_summary.json` and
`pattern_summary.json` retain numerical-kernel hashes. The seven paper-facing notes are
generated by `scripts/write_reviewer_reports.py`, and all eight new figure pairs by
`figures/plot_reviewer_extensions.py`. See [the index](SUPPLEMENTARY_INDEX.md).

| Main result | Source data | Computation |
|---|---|---|
| 52 QSS families, 3172 evaluations; P2/P3 10% alpha thresholds 199.53/398.11 | qss_validity.csv, qss_accuracy_thresholds.csv, qss_families.csv | reviewer_linear.py |
| 27 of 677 q <= 0.01 evaluations still exceed 10% relative growth error | qss_validity.csv | reviewer_linear.py |
| P2 clock 0.06451837, P3 error -16.5105%; held-out LS errors -7.48% to -6.77% | clock_recalibration.csv, clock_mode_predictions.csv, clock_frequency_response.csv | reviewer_linear.py |
| All 12 DC misses enter via oscillatory linear boundaries, 6 later real-pair collisions | oscillatory_boundaries.csv, real_pair_collisions.csv, bifurcation_poles.csv.gz | reviewer_linear.py |
| 23/24/25 critical lambda_N 1.308113653/1.304815343/1.307902938 | geometry_stationary_thresholds.csv | reviewer_geometry.py |
| All 18 periodic pairs retain later full onset; 52 geometry/control trajectories | geometry_pairs.csv, geometry_runs.csv, geometry/*.npz | reviewer_geometry.py |
| 12 complete spectra stable; max leading Re -0.0976573; max recovery ratio 1.831e-9 | pattern_stability.csv, pattern_stability/*.npz | reviewer_pattern_stability.py |
| Margin >= 0.05: 4048 points, median ratio 3.45319, median absolute gap 0.82313 | growth_boundary_strata.csv, growth_error_points.csv.gz | reviewer_linear.py |
| Reweighted fractions 39.1287% / 38.3219%, ESS 22497.85 / 16699.60 | prior_sensitivity.csv, prior_replicates.csv | reviewer_linear.py |
| Exact production-lag similarity checked on 216000 blocks; max residual 1.111e-16 | lag_similarity_checks.csv | reviewer_linear.py and reviewer_models.py |
| Clamped-A factors 1.001/1.2/2: normalised differences 0.0631%/11.9854%/31.9231% | mechanism_step_responses.csv | reviewer_linear.py |
| Separate Hill thresholds absorbed by S-unit rescaling; 100 state checks | separate_Hill_threshold_checks.csv | reviewer_linear.py |

All spectrum and percentage calculations are continuous-time and have no integration dt.
Oscillatory linear boundaries are Hopf candidates; no cycle or nonlinear Hopf claim is made.
Full versus reduced sampling findings are empirical; alternative lag linear equivalence
and Hill rescaling are algebraic. Complete patterned spectra establish local stability only
for the selected endpoints. Prior reweighting has unchanged support and is illustrative.
