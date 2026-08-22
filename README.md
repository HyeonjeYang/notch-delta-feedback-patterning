# Notch–Delta Feedback Patterning

This repository studies spatial symmetry breaking in an explicit
Delta–Notch–activated-Notch–proneural (D–N–A–ASC) feedback architecture,
using linear stability analysis, spatial eigenmode decomposition, control
theory, and nonlinear deterministic simulation. It contains the frozen,
validated deterministic analysis behind a preprint, not the exploratory
history that produced it.

## Model

State variables: `D` = Delta, `N` = Notch, `A` = activated Notch /
NICD-like downstream activity, `S` = ASC / proneural activity. The
biological loop: `S_i -> D_i` (proneural drives Delta); Delta on a
neighbouring cell activates that neighbour's Notch/`A` (`neighbour D ->
Notch/A`, trans interaction); `A_i` represses `S_i` (`A_i -| S_i`),
closing the loop. A separate cis interaction couples `D` and `N` within
the same cell.

The exact implemented equations (`model/equations.py`,
`model/jacobian.py`), unchanged from the source repository:

```
D_dot = lam_D * h_D(S) - d_D*D - f_D*D*s_N - k_cis*N*D
N_dot = lam_N          - d_N*N - f_N*N*s_D - k_cis*N*D
A_dot = -d_A*A + f_N*N*s_D
S_dot = -d_S*S + lam_S*h_S(A)

h_D(S) = 1 / (1 + (Kd/S)^n_d)        (increasing in S)
h_S(A) = 1 / (1 + (A/Kd)^n_ASC)      (decreasing in A)

s_D = sum_{j~i} D_j,   s_N = sum_{j~i} N_j   (six-neighbour contact sums)
```

## Geometry

Cells are represented on a **periodic six-neighbour triangular contact
graph corresponding to a hexagonal cellular tiling** (`geometry/lattice.py`,
`geometry/graph_spectrum.py`). Its adjacency operator `W` is diagonalized
by its own spatial eigenmodes, `W v_mu = mu v_mu`. The uniform mode has
`mu = 6`; the dominant, three-sublattice pattern-forming mode has
`mu = -3` (the graph's minimum eigenvalue, attained only at the two
three-fold-periodic "K-like" points).

## Linear stability

A perturbation along one spatial eigenmode decouples the full
multicellular linearization into one 4x4 ODE per graph eigenvalue
(`model/jacobian.py`, `stability/mode_scan.py`):

```
d xi_mu / dt = J_mu xi_mu
det(sI - J_mu) = 0
```

The leading eigenvalue/pole of `J_mu` determines whether a spatial
perturbation at mode `mu` decays or grows.

## Control-theory interpretation

Partitioning `x = (D, N, A, S)` into `z = (D, N, A)` and the scalar `S`
gives the exact state-space form actually used in the analysis
(`control/state_space.py`):

```
z_dot = A_mu z + b S
S_dot = c^T z + a_S S

P_mu(s) = c^T (sI - A_mu)^{-1} b
L_mu(s) = P_mu(s) / (s - a_S)

det(sI - J_mu) = det(sI - A_mu) * (s - a_S) * [ 1 - L_mu(s) ]
```

verified against the direct 4x4 Jacobian eigenvalues to machine precision
(max coefficient difference `1.4e-13`, `control/schur_verification.py`).

**Key result**: the stationary pattern-forming bifurcation occurs when the
exact DC loop gain of the dominant three-sublattice mode reaches unity:

```
G_exact(-3) = L_{-3}(0) ~= 1   at   lam_N ~= 1.305
```

This was derived from the exact Schur-complement representation above, not
assumed, and independently confirmed by an unconstrained root-locus pole
crossing at the same `lam_N` (`control/root_locus.py`). Full derivation:
`docs/CONTROL_THEORY.md`.

## Main result

- Downstream Notch/proneural regulatory kinetics (`d_A, lam_N, Kd, n_ASC,
  d_ASC, lam_ASC`) — not the classical cis/trans Delta-Notch coupling
  (`k_cis, f_D, f_N`) — dominate the spatial stability margin in this
  realization.
- The dominant pattern-forming mode is the `mu = -3` three-sublattice mode
  on the six-neighbour contact graph, not a square-checkerboard mode.
- Geometry enters as a mode-dependent spatial gain: the same biochemical
  `S->D->A->S` loop is net-stabilizing at `mu=6` and net-destabilizing at
  `mu=-3`, purely through the sign of the graph eigenvalue.
- The bifurcation is an exact unity-DC-loop-gain condition,
  `G_exact(-3)=1`, not merely a numerically-observed sign change.
- Nonlinear deterministic growth from a tiny (1e-4 relative) perturbation
  matches the predicted unstable pole and saturates into a persistent,
  finite-amplitude three-sublattice pattern (validated at `lam_N=2.0` and
  `5.623`, zero numerical-integration clipping).

## Running the code

```bash
pip install -r requirements.txt
pytest                              # verifies the refactored code reproduces the frozen numbers
python figures/plot_fig1.py         # and plot_fig2.py .. plot_fig5.py
```

Each `figures/plot_fig*.py` reads only the small CSV/NPZ files already
included in `figures/data/` — no simulation or parameter search is
rerun. Figures are written to `figures/main/`.

```bash
python control/explain_geff.py     # optional: recomputes figures/data/geff_comparison.csv
                                    # from the included raw global sample (no new search)
```

To rerun a piece of the deterministic analysis itself (fast; no search):

```python
from model.params import P0, LAM_N_P2
from model.fixed_point import find_fixed_points
from stability.mode_scan import lambda_pattern
import dataclasses

p = dataclasses.replace(P0, lam_N=LAM_N_P2)
fp = find_fixed_points(p)[0]
print(lambda_pattern(fp, p))   # {'nonuniform_max_re_lambda': 0.103..., ...}
```

## Repository map

| directory | contents |
|---|---|
| `model/` | kinetic parameters, drift equations, homogeneous fixed-point solver, mode-resolved Jacobian, order-parameter metrics |
| `geometry/` | six-neighbour contact lattice, graph spectrum / dispersion relation, reciprocal-space growth map |
| `stability/` | spatial-mode scan (`lambda_pattern`), local sensitivity, 1-D parameter continuation |
| `control/` | Schur-complement state-space form, algebraic verification, root locus, exact/heuristic loop gain |
| `simulation/` | validated nonlinear deterministic ODE integration and controlled eigenmode perturbations |
| `figures/` | scripts for the five main figures, plus the small source data they read |
| `docs/` | condensed control-theory and geometry derivations, manuscript numbers |
| `tests/` | fidelity checks against the frozen source analysis |

## Reproducibility note

This repository contains the **frozen deterministic analysis** used for
the preprint: validated equations, fixed points, linear stability, control
theory, geometry, and nonlinear confirmation at two validated parameter
points (`lam_N = 2.0`, `5.623`). It intentionally excludes the exploratory
history that produced it — in particular, an earlier `dt=0.02`
Euler-Maruyama result at a different (`k_cis`-only) baseline was found to
be a numerical integration artifact and is **not** reproduced here. One
preliminary stochastic check has since been added (`stochastic_check/`,
scripts only — results live in the main development repository per
`stochastic_check/README.md`): at validated point P2 (`lam_N=2.0`),
additive white noise on Delta shortens formation time (`t50` -45% to -66%
across `sigma=1e-4..1e-2`) without changing final pattern amplitude, but
collapses final three-sublattice mode purity from 1.0 to ~0.13-0.22
(verdict S2); a deterministic random-IC control found ~60% of that purity
loss already occurs with zero noise, with ongoing noise contributing a
further, significant ~40% of the acceleration (verdict R3). This is one
check at one parameter point, not a systematic noise-timing study across
the bifurcation, which remains future work. The frozen
4000-point global parameter search that identified the dominant drivers
(`figures/data/global_parameter_samples.csv`) is included in full; its
root-finding step itself is not rerun here (no new search), but
`control/explain_geff.py` recomputes the `G_exact`-vs-`G_eff` comparison
from it and reproduces `figures/data/geff_comparison.csv` byte-for-byte.
