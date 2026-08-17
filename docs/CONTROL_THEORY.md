# Control-theoretic mechanism

Condensed from the source repository's `preprint_final_analysis/CONTROL_THEORY_DERIVATION.md`.
Full provenance, algebraic verification numbers, and the G_eff explanation
are reproduced here; see `../MANUSCRIPT_NUMBERS.md`-equivalent values in
`docs/MANUSCRIPT_NUMBERS.md` for the exact source of every number.

## State-space partition

At a homogeneous fixed point, a perturbation along graph eigenmode `mu`
obeys `d xi_mu/dt = J_mu xi_mu` (`model/jacobian.py`). Partition
`x = (D, N, A, S)` into `z = (D, N, A)` and the scalar `S`:

```
z_dot = A_mu z + b S
S_dot = c^T z + a_S S

A_mu = J_mu[0:3, 0:3]    b = J_mu[0:3, 3]    c^T = J_mu[3, 0:3]    a_S = J_mu[3, 3]
```

(`control/state_space.py::partition`).

## Exact Laplace-domain feedback representation

```
P_mu(s) = c^T (sI - A_mu)^{-1} b
L_mu(s) = P_mu(s) / (s - a_S)

det(sI - J_mu) = det(sI - A_mu) * (s - a_S) * [ 1 - L_mu(s) ]
```

The sign (`1 - L`, not `1 + L`) was derived from the block-determinant
expansion, not assumed. Verified numerically (`control/schur_verification.py`)
against the direct 4x4 Jacobian eigenvalues for `mu in {6, -3, 1.5}` at
`lam_N in {1, 2, 5.623}`: **max coefficient difference 1.4e-13, max root
mismatch 4.7e-14** -- machine precision.

## Root locus and the exact bifurcation condition

Tracking all four poles of `J_mu` across a `lam_N` sweep
(`control/root_locus.py`): the `mu=-3` (dominant, three-sublattice) leading
pole crosses `Re(s)=0` at **`lam_N ~= 1.305`**, while `mu=6` (uniform) stays
stable (max `Re(s)=-0.1`) across the whole tested range `[0.5, 6.5]`.

Define `G_exact(mu) = L_mu(0)` (`control/state_space.py::G_exact`). At the
critical `lam_N`:

```
G_exact(-3) = 0.999936   -- an exact unity-DC-loop-gain condition, to 0.007%
G_exact(+6) = -0.063368  -- far from unity, consistent with uniform-mode stability
```

This precision was not built in; it is a direct consequence of the
Schur-complement identity being exact algebra, checked against an
independently-located pole crossing.

## G_eff: a reduced approximation, explained

The heuristic used in earlier analysis,

```
G_eff = J_D,S * J_A,D * J_S,A / (|J_D,D| * |J_A,A| * |J_S,S|)
```

(`control/loop_gain.py::G_eff`) is the dominant-path (`S->D->A->S`)
approximation of `G_exact`: it drops the N-mediated coupling terms and
approximates `[A_mu^{-1}]_{A,D}` by `J_AD(mu)/(J_DD*J_AA)`. Over the frozen
global parameter sample (n=3267 valid points,
`figures/data/global_parameter_samples.csv`, included in full --
`control/explain_geff.py` reproduces the comparison below byte-for-byte
from this file, no new search):

| criterion | classification accuracy |
|---|---|
| `G_exact(-3) > 1` | **99.79%** |
| `G_eff` (best threshold) | 93.88% |

`G_eff` should be treated as an explanatory heuristic, not cited as exact.

## Spatial feedback sign (interpretation of the dominant path)

```
S -> D            positive   (h_D increasing in S)
D_neighbour -> A  sign(mu)   (J_AD(mu) = f_N*N* * mu)
A -> S            negative   (h_S decreasing in A)
```

The three-link product flips sign with `mu`: net-stabilizing at the
uniform mode (`mu=6`), net-destabilizing at the dominant pattern mode
(`mu=-3`) -- purely a consequence of the graph eigenvalue, not any change
in biochemical rate constants. This is the precise sense in which "tissue
geometry acts as a mode-dependent spatial gain in the biochemical feedback
loop." Stated as an interpretation of the dominant path; the full,
un-reduced mechanism is `G_exact(mu)` itself.
