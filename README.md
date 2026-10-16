# Notch-Delta feedback patterning

Analysis code, data and figures for a study of spatial symmetry breaking in an explicit
Delta-Notch-activated-Notch-proneural (D-N-A-S) lateral-inhibition model on a periodic
six-neighbour contact lattice.

The repository holds the **revised** analysis. The version that accompanied the first
manuscript submission is preserved unchanged in
[`archive/v1_physbio_submission/`](archive/v1_physbio_submission/) (tag `v1-physbio-submission`).
[`archive/README.md`](archive/README.md) lists the problems found in review and how each was
corrected.

## Main results

1. **One equilibrium.** For every positive parameter set the homogeneous equilibrium exists,
   is unique, and is found by scalar bracketing. In the 4000-point global sample,
   1453/4000 = 36.3% of points are unstable to a spatial mode. The v1 figure was 30.3%; the
   difference came from 733 solver failures that v1 counted as nonexistent equilibria.
2. **Static loop gain.** Split each mode-resolved Jacobian `J_mu` into the plant
   `(D, N, A)` and the feedback state `S`. The static loop gain `G(mu)` then satisfies
   `det J_mu = det(A_mu) a_S [1 - G(mu)]`. The plant is Hurwitz for every mode, so
   `G(mu) = 1` marks a zero eigenvalue and `G(mu) > 1` an odd number of positive real
   eigenvalues. This identity is standard (Schur complement, or equivalently graph
   input-output static gain). It is used here as a tool, not claimed as a result.
3. **Mode selection.** `G(mu)` is largest at the most negative graph eigenvalue for every
   positive parameter set. Stationary onset therefore always selects the three-sublattice
   (K-point) mode, and the uniform mode can never lose stability through a zero eigenvalue.
4. **Where `G` is not enough.** 12 of the 4000 points are unstable with `G(-3) < 1`. Each has
   exactly two right-half-plane eigenvalues (a complex pair or two real eigenvalues) in the
   `mu = -3` block. The Nyquist criterion counts them correctly at all 4000 points.
5. **What the explicit downstream layer adds.** Eliminating `S` (or `A` and `S`) at quasi-steady
   state keeps every equilibrium and the entire stationary boundary. The growth rate and the
   non-stationary instabilities change. Rescaling the proneural timescale leaves
   `lam_N^c = 1.304815343` unchanged while the P2 growth rate spans three decades. The
   matched reduction grows faster at every unstable point of the global sample (median
   4.9-fold; 15.5-fold at P2), yet both models reach the same final pattern.

Proofs are in [`docs/THEORY.md`](docs/THEORY.md). The source of every number is in
[`docs/MANUSCRIPT_NUMBERS.md`](docs/MANUSCRIPT_NUMBERS.md).

## Model

```text
dD/dt = lam_D h_D(S) - d_D D - f_D D s_N - k_cis N D
dN/dt = lam_N        - d_N N - f_N N s_D - k_cis N D
dA/dt =              - d_A A + f_N N s_D
dS/dt =              - d_S S + lam_S h_S(A)

h_D(S) = 1/(1 + (K_d/S)^n_D),   h_S(A) = 1/(1 + (A/K_d)^n_S),   s_X = sum over the six neighbours of X
```

## Running

```bash
pip install -r requirements.txt
pytest                                   # 340 checks, about 4 s
python scripts/reproduce_all.py          # every data file and figure, about 15 min
python scripts/reproduce_all.py --figures  # figures only, from figures/data/
```

## Layout

| path | contents |
|---|---|
| `model/` | parameters, drift, unique homogeneous equilibrium, mode-resolved Jacobian, matched QSS reductions |
| `geometry/` | six-neighbour contact lattice and its spectrum |
| `stability/` | mode scan, sensitivities, one-parameter continuations |
| `control/` | feedback partition and static gain, Nyquist count, stage gains, graph input-output form, root locus, Schur check |
| `simulation/` | lattice ODE integration (full and DNA reduction), eigenvector initial states, K-mode amplitude |
| `scripts/` | one script per data set (see `scripts/reproduce_all.py` for the order) |
| `figures/` | `plot_fig1.py`-`plot_fig5.py`, `plot_supplementary.py`, shared `style.py`; `data/`, `main/`, `supplementary/` |
| `docs/` | `THEORY.md` (propositions and proofs), `GEOMETRY.md`, `MANUSCRIPT_NUMBERS.md` |
| `manuscript/` | revised manuscript draft and response to reviewers |
| `tests/` | validation suite |
| `archive/` | frozen earlier versions (not imported or tested) |

Figure 1A is left blank on purpose. The circuit diagram is drawn by hand and overlaid.
