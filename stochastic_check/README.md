# stochastic_check

Code only. One minimal, prespecified stochastic experiment on top of the
already validated deterministic P2 point (`lam_N=2.0`), plus a follow-up
control distinguishing initial-condition effects from ongoing-noise
effects. Not a parameter search, not a restart of any prior
stochastic-timing project.

**Results, figures, and data (`FINAL_REPORT.md`, `FINAL_CONTROL_REPORT.md`,
`Figure_*.{png,pdf}`, `*.csv`, `*.npz`) are not kept in this repository** --
they live in the main development repository,
`Neuroblast_differentiation_simulation/stochastic_check/` (same relative
layout, `final_control/` subdirectory included). This directory holds only
the four scripts that produced them.

**Question**: at P2, does weak-to-moderate additive white noise on Delta
appreciably change the rate or quality of emergence of the
already-deterministically-unstable three-sublattice pattern -- and is any
such effect noise-specific, or mainly an artifact of comparing a
symmetry-protected pure-eigenmode deterministic initial condition against
a multimode state?

## Scripts

| file | role |
|---|---|
| `run_stochastic_check.py` | main sweep: additive white noise on `D` only (`dD_i = F_D,i(x) dt + sigma*D_star*dW_i`, Euler-Maruyama), 4 prespecified `sigma` in `{0, 1e-4, 1e-3, 1e-2}` x 12 replicates, `dt=0.005`, `T=600`, plus one `dt=0.0025` convergence control. Uses `model.equations.drift` unmodified. |
| `final_control/analyze_structure_factor.py` | replays the already-completed sweep (identical seeds, verified byte-for-byte against the saved summary before use -- not a new run) to compute the mean-subtracted discrete structure factor `S(k)` of every final Delta field, and verifies its exact-K power reproduces the `P_mu=-3` projection metric to machine precision. |
| `final_control/run_random_ic_control.py` | 12 purely deterministic control runs (DOP853, `rtol=1e-9, atol=1e-11`, zero noise) at the same P2 point, replacing only the pure-eigenmode initial perturbation's spatial profile with a norm-matched random field, to separate initial-condition effects from ongoing-noise effects. |
| `final_control/plot_final_control.py` | the one required final-control figure, reading only already-computed CSV/NPZ data. |

## Definitions (shared across all four scripts)

- `q_amp`, `std(D)`, Moran's I: `model.metrics.order_parameters`-equivalent
  (re-derived locally in `run_stochastic_check.py` for exact consistency).
- `P_mu=-3(t)`: `|Pi_{-3}(D(t)-mean(D(t)))|^2 / |D(t)-mean(D(t))|^2`, `Pi_{-3}`
  the QR-orthonormalized 2-D real eigenspace of the dominant `mu=-3` mode
  (`cos`/`sin` of the `(m,n)=(8,16)` plane wave).
- `t50`/`t90`: first *recorded* time (`RECORD_INTERVAL=1.0`) `q_amp`
  reaches 0.50/0.90; `NaN` if never reached by `T=600`. Single crossing,
  no persistence rule, prespecified before any run.
- sender/receiver ratio: `D > mean(D)` = sender-like, else receiver-like.
- `clip_fraction`: fraction of `(cell, step)` Delta updates that went
  negative before projection (main sweep only -- the deterministic control
  has no noise and no clipping step).

## Reproducing

```bash
python stochastic_check/run_stochastic_check.py                       # ~10 min, 48 runs + dt control
python stochastic_check/final_control/analyze_structure_factor.py     # ~5 min, replay for spectral analysis
python stochastic_check/final_control/run_random_ic_control.py        # ~2 min, 12 deterministic DOP853 runs
python stochastic_check/final_control/plot_final_control.py           # seconds
```

No parameter search anywhere in this chain; every sigma/seed/dt value is
prespecified in the scripts themselves, not chosen after inspecting
results.
