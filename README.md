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
   Ten independent 4000-point LHS replicates give 37.37% (approximate 95% Monte Carlo
   interval 37.07–37.67%) under the same parameter measure.
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
   matched reduction grows faster at every unstable point of the reference sample (median
   4.9-fold; 15.5-fold at P2); the pooled median is 5.41 in the independent 40000-point study.
   The specified pure-mode P2 runs reach the same final pattern; paired random starts
   give different endpoint fields.
6. **Sampling robustness.** 115000 new evaluations cover independent seeds, non-nested
   sample-size checks and a wider parameter box. Qualitative dynamic comparisons persist;
   changing the bounds changes the unstable fraction and growth-ratio magnitude. Every input,
   per-point result, seed-level summary and the protocol is retained under `figures/data/robustness/`.
7. **Initial-condition robustness.** 112 deterministic trajectories span two operating points,
   eight seeds and three amplitudes. In all 48 pairs the full model reaches CV(D)=0.1 later:
   10.96–15.10 times later at P2 and 12.34–12.90 at P3. Endpoint fields differ in every pair,
   including after periodic translation. Fixed-horizon residual checks do not prove attractor
   stability. Figure S5 and Supplementary Note 3 report the design and solver checks.
8. **Accuracy and mechanistic scope.** A reduced-model clock fitted at P2 underpredicts
   P3 growth by 16.5%. All twelve reference gain-below-one cases enter instability through
   oscillatory linear boundaries on the tested constant-gain paths; six later become real
   unstable pairs. A matched Delta-production lag is exactly equivalent for linear D,N,A
   responses, but differs by 31.9% of the total production response after a twofold clamped-A
   step. Linear growth identifies finite memory without uniquely locating it in proneural S.
   Main Figure 6 and [seven separate supplements](docs/SUPPLEMENTARY_INDEX.md) present these
   findings and QSS, geometry, pattern-stability and prior-sensitivity controls.

The Discussion and Table 3 distinguish the matched-control contribution from earlier
biological models and graph results. Proposed tests change proneural production and degradation
together to hold its static response fixed, or measure its finite response to clamped activated
Notch. The model has no tissue calibration; these predictions require time-resolved validation.

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
pip install -r requirements-lock.txt     # validated with CPython 3.13.5 on Windows x86_64
pytest                                   # numerical and sampling reproducibility checks
python scripts/reproduce_all.py          # reference analysis, robustness, audits and figures
python scripts/reproduce_all.py --figures  # figures only, from figures/data/
python scripts/reproduce_all.py --audit   # audit stored results and numerical settings
python scripts/reproduce_all.py --robustness # independent seeds, sizes, bounds and paper supplement
python scripts/reproduce_all.py --initial-conditions # paired random starts, predictions and Figure S5
python scripts/reproduce_all.py --dt-audit # actual step refinement, independent solvers and Figure S6
python scripts/reproduce_all.py --reviewer-extensions # seven matched controls, Notes 5–11, Figures 6/S7–S13
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
| `figures/` | main Figures 1–6 and supplementary S1–S13, plotting scripts and shared `style.py`; `data/`, `main/`, `supplementary/` |
| `docs/` | `THEORY.md` (propositions and proofs), `GEOMETRY.md`, `MANUSCRIPT_NUMBERS.md` |
| `manuscript/` | revised manuscript draft and response to reviewers (`.md` sources; `python manuscript/build.py` builds `.docx` and `.pdf`) |
| `tests/` | validation suite |
| `archive/` | frozen earlier versions (not imported or tested) |

Figure 1A is intentionally blank for the author's separately prepared circuit schematic.
All figures use a shared blue, orange, teal and ochre palette; the final-field map uses
a sequential cream-to-blue scale.

`python scripts/audit_revision.py` independently recalculates the 4000-point statistics and
writes the complete per-mode scan. Subdominant nonuniform modes are also unstable at 1451
points; dominance of the K mode does not mean it is the only unstable mode (Fig. S3).

The Markdown manuscript and response are the current sources. The checked-in DOCX files
predate this audit and must be rebuilt and visually checked before submission.

[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) specifies both random-sampling protocols,
the validated environment, input snapshots, numerical tolerances and run manifests.

[`manuscript/robustness_supplement.md`](manuscript/robustness_supplement.md) contains the
paper-facing protocol, results, Table S1 and Figure S4. Its tables are regenerated from data
by `scripts/write_robustness_report.py`. The main manuscript and response use the same results.

[`manuscript/initial_condition_supplement.md`](manuscript/initial_condition_supplement.md)
contains Supplementary Note 3, Tables S2–S3 and Figure S5. Its paired-field data and
protocol history are retained in `figures/data/initial_conditions/`. The Discussion's
compensated-turnover values are regenerated by `scripts/mechanistic_predictions.py`.

[`manuscript/timestep_supplement.md`](manuscript/timestep_supplement.md) contains
Supplementary Note 4, Tables S4–S6 and Figure S6. The frozen protocol and results are
under `figures/data/timestep_audit/`. This distinguishes output cadence from internal
adaptive steps, adds genuine fixed-step RK4 refinement for early growth, and tests
long-time endpoint agreement with step-capped BDF and independent Radau integration.
The 70 full-horizon checks change emergence times by at most 0.01631% and raw endpoint
relative L2 field differences remain below 1.12e-10. Coarse RK4 controls at P0 do fail: h=0.04 and
0.08 produce artificial instabilities. These failed controls are retained, so the combined
audit flag is false while all adaptive production-comparison targets pass; see
`figures/data/timestep_audit/assessment_interpretation.json`. No target was relaxed.
