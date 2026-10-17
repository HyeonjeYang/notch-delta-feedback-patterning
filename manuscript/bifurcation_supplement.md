# Supplementary Note 8. Oscillatory routes to gain-below-one instability

## Constant-gain continuation of every reference DC miss

The corrected reference sample has twelve unstable cases with $G(-3)<1$: six have an
unstable complex pair and six have two positive real poles at their original parameters.
All have dominant mode $\mu=-3$. To identify a route into these failures **within the tested
mode**, we follow $(\lambda_S,d_S)\mapsto\alpha(\lambda_S,d_S)$ at every one of the twelve
points. This holds the complete equilibrium and $G(\mu)$ fixed. We evaluate all 54 mode
blocks, including the uniform mode, over 401 logarithmically spaced $\alpha$ from
$10^{-6}$ to $10^2$, bracket the first stable-to-unstable crossing, and refine in $\log\alpha$.

All twelve first crossings occur in the K mode at nonzero frequency. The real-part residual
at a refined crossing is at most $1.57\times10^{-15}$; the derivative of the leading real
part with respect to $\log\alpha$ is positive in every case. Before and after crossing,
right-half-plane counts are computed over the 54 distinct blocks (without spatial
multiplicity). **We call these oscillatory linear boundaries or Hopf candidates.** We have
not computed a first Lyapunov coefficient or proved a nondegenerate nonlinear Hopf
bifurcation, the direction of a cycle branch, or cycle stability.

**Table S9. All first oscillatory boundaries on constant-gain paths.**

| reference ID | $\alpha$ | $\mu$ | $|\mathrm{Im}\,s|$ | $d\mathrm{Re}\,s/d\log\alpha$ | $G(-3)$ |
|---|---|---|---|---|---|
| 679 | 0.178925 | -3 | 0.106291 | 0.0894917 | 0.563539 |
| 1463 | 0.68719 | -3 | 0.0389695 | 0.0371507 | 0.581764 |
| 1543 | 0.168643 | -3 | 0.0905597 | 0.153262 | 0.818668 |
| 1669 | 0.288803 | -3 | 0.0178484 | 0.0210775 | 0.949905 |
| 1730 | 0.446628 | -3 | 0.235521 | 0.170988 | 0.365644 |
| 1829 | 0.081302 | -3 | 0.0369438 | 0.0597652 | 0.814166 |
| 2972 | 0.687272 | -3 | 0.0714944 | 0.117861 | 0.665258 |
| 3161 | 0.205375 | -3 | 0.0577945 | 0.0308902 | 0.449546 |
| 3165 | 0.0952798 | -3 | 0.0140526 | 0.0184115 | 0.917413 |
| 3746 | 0.34937 | -3 | 0.0390003 | 0.0870913 | 0.93095 |
| 3822 | 0.836747 | -3 | 0.0404116 | 0.0348826 | 0.780425 |
| 3919 | 0.473626 | -3 | 0.109471 | 0.0776732 | 0.794066 |

Six unstable complex pairs subsequently collide on the positive real axis before $\alpha=1$.
We locate this collision by the squared splitting of the leading pair. For example, sample
679 first crosses at $\alpha=0.1789253$ and frequency 0.1062911, then becomes a real pair at
$\alpha=0.7808506$, with mean real part 0.165840. Its $G(-3)=0.563539$ throughout.
Thus the real-pair misses need not enter instability through a zero eigenvalue on this path.
The Schur stationary condition and Nyquist unstable-pole count remain consistent.

**Table S10. Later unstable-pair collisions.**

| reference ID | $\alpha$ | collision real part | numerical pair splitting | $G(-3)$ |
|---|---|---|---|---|
| 679 | 0.780851 | 0.16584 | 1.77374e-08 | 0.563539 |
| 1543 | 0.31062 | 0.11712 | 7.3578e-08 | 0.818668 |
| 1669 | 0.86176 | 0.022821 | 1.6246e-08 | 0.949905 |
| 1829 | 0.150339 | 0.049732 | 1.99199e-08 | 0.814166 |
| 3165 | 0.205755 | 0.0199795 | 2.25178e-08 | 0.917413 |
| 3746 | 0.564796 | 0.045814 | 6.76503e-09 | 0.93095 |

At a collision, repeated eigenvalues are ill-conditioned. Stored tiny pair splittings are
numerical diagnostics, not a claim of exact root equality to that many digits. These paths
are constructed controls; no assertion is made that a biological tissue follows them or
that every possible route to the same endpoint must be oscillatory.

## Reproduction

`scripts/reviewer_linear.py` writes all K poles along the paths to
`bifurcation_poles.csv.gz`, the refined boundaries to `oscillatory_boundaries.csv` and
collisions to `real_pair_collisions.csv`. Use the report and figure commands in Note 5.

The design was recorded in `figures/data/reviewer_extensions/protocol.json` before the
production calculations. All selected cases and counterexamples are retained. Source hashes
and case checksums are verified by `scripts/write_reviewer_reports.py`. The primary package
versions and numerical reproduction contract are given in `docs/REPRODUCIBILITY.md`.

![**Figure S10. Constant-gain routes. (A,B) K poles for samples 1730 and 679; blue shows real parts, orange absolute imaginary parts. Vertical lines mark the first oscillatory boundary and, where present, later real-pair collision. (C) Crossing frequency for every miss. (D) Boundaries and collisions at unchanged G below one.**](../figures/supplementary/FigS10_oscillatory_routes.png){width=100%}
