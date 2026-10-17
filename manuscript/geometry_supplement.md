# Supplementary Note 7. Finite lattice, open boundaries and static heterogeneity

## Periodic size and commensurability

The original 24-cell side is divisible by three. To test dependence on this choice, we use
23, 24 and 25 cells per side, P2/P3, three PCG64 seeds (2026102001–2026102003), and full/DNA
models: 36 runs, 18 matched pairs. Each pair has identical initial $D,N,A$ fields with
relative amplitude $10^{-4}$, and full-model $S=S_q(A)$. Random fields are spatially centred
and normalised as in Note 3. Periodic graph sums retain six neighbours.

**Table S8. Actual minimum graph eigenvalue and stationary onset.**

| periodic side length | $\mu_{\min}$ | multiplicity | $\lambda_N^c$ |
|---|---|---|---|
| 23 | -2.97648 | 6 | 1.30811 |
| 24 | -3 | 2 | 1.30482 |
| 25 | -2.97797 | 6 | 1.3079 |

The maximum relative threshold change from 24 is approximately 0.253%. Thus exact
$\mu=-3$ commensurability changes onset slightly but is not required for instability or
the timing separation at these operating points. Across all 18 pairs, the full/DNA time
ratio to Delta CV 0.1 is
12.7457–14.7206.

Endpoint structures vary. To compare spectra across sizes, we report both exact-minimum-shell
power and power in the fixed band $\mu\leq-2.8$. Neither statistic implies a perfectly
three-sublattice final pattern. P0 is a stable periodic control in both models.

## Open edges and quenched heterogeneity

Six further runs use a 24-by-24 **open** graph at P0/P2/P3 in both models. Neighbours outside
the array are dropped; degrees range from two to six. Raw sums are retained without degree
normalisation. The periodic homogeneous reference is therefore **not an equilibrium** of
this open graph. P0 approaches a nonuniform boundary-forced state (Delta CV about 0.0698)
even though its periodic control returns to homogeneity. This is a useful negative control
against calling every nonzero CV a spontaneous instability. We report three-cell-rim and
interior CV separately. A degree-normalised open graph would be a different coupling model.

Eight runs add time-independent, mean-one lognormal Notch-production multipliers with
target CV 0.05 or 0.1 at P2/P3, using a separate PCG64(seed+1) stream; the realised multipliers
are stored and have sample mean exactly one. These are **quenched parameter differences**,
not dynamical noise. Full-model CV crossings remain later in all four pairs, but deterministic
forcing advances pattern contrast and reduces the time ratio. Neither open nor heterogeneous
crossing times are interpreted as homogeneous-instability thresholds.

All 52 trajectories use BDF with an analytic sparse Jacobian, rtol $10^{-10}$,
atol $10^{-13}$, maximum internal step 0.5, and recording interval 0.5. Horizons are
600 at P0/P2 and 300 at P3. All recorded concentrations are positive; the maximum final
scaled drift is $7.55\times10^{-10}$. These finite controls do not constitute a survey of
arbitrary contact graphs, tissue growth or stochastic dynamics.

## Reproduction

Run `python scripts/reviewer_geometry.py --reuse`, then the report and figure commands.
Omit `--reuse` for fresh trajectories. `geometry_runs.csv`, `geometry_pairs.csv`,
`geometry_stationary_thresholds.csv`, and all 52 initial/final field files are retained.

The design was recorded in `figures/data/reviewer_extensions/protocol.json` before the
production calculations. All selected cases and counterexamples are retained. Source hashes
and case checksums are verified by `scripts/write_reviewer_reports.py`. The primary package
versions and numerical reproduction contract are given in `docs/REPRODUCIBILITY.md`.

![**Figure S9. Geometry controls. (A) Stationary thresholds. (B) All 18 paired periodic onset-time ratios. (C) Fixed near-K-band endpoint power at both operating points. (D) Stable P0: periodic decay versus open-edge forcing; open-graph contrast is not a homogeneous-instability test.**](../figures/supplementary/FigS9_geometry_robustness.png){width=100%}
