# Theory: equilibrium, static loop gain, mode selection, and matched reductions

This file states the analytic results used in the revised manuscript and proves them. Every
statement is also checked numerically in `tests/test_core.py`.

## Notation

Each cell `i` carries Delta `D_i`, Notch receptor `N_i`, activated Notch / NICD-like activity
`A_i`, and proneural activity `S_i`. Cells are the vertices of a connected `z`-regular contact
graph with adjacency matrix `W`. On the periodic triangular lattice used here, `z = 6`. The neighbour sums are
`s_D,i = sum_j W_ij D_j` and `s_N,i = sum_j W_ij N_j`.

```
dD/dt = lam_D h_D(S) - d_D D - f_D D s_N - k N D
dN/dt = lam_N        - d_N N - f_N N s_D - k N D
dA/dt =              - d_A A + f_N N s_D
dS/dt =              - d_S S + lam_S h_S(A)

h_D(S) = 1 / (1 + (K_d/S)^n_D)    increasing
h_S(A) = 1 / (1 + (A/K_d)^n_S)    decreasing
```

`k` is `k_cis`, `lam_S` and `d_S` are `lam_ASC` and `d_ASC`, and `n_D` and `n_S` are `n_d`
and `n_ASC` in the code. All 13 parameters are positive.

`W` is symmetric with eigenvalues `mu` in `[-z, z]`. Linearising about a homogeneous
equilibrium `x* = (D, N, A, S)` and projecting onto an eigenvector of `W` with eigenvalue `mu`
gives one 4x4 block `J_mu` per graph eigenvalue (`model/jacobian.py`):

```
        [ a11       -D(k+f_D mu)   0      kappa_D ]
J_mu =  [ -N(k+f_N mu)   a22       0       0      ]
        [ f_N N mu     z f_N D    -d_A     0      ]
        [   0             0       -c_S   -d_S     ]

a11 = -(d_D + (z f_D + k) N),   a22 = -(d_N + (z f_N + k) D),
kappa_D = lam_D h_D'(S) > 0,    c_S = lam_S |h_S'(A)| > 0.
```

On the 24 x 24 periodic lattice the eigenvectors are plane waves indexed by integers `(m, n)`,
`0 <= m, n < 24`, with `mu(m, n) = 2 cos k_1 + 2 cos k_2 + 2 cos(k_1 - k_2)`, where
`k_1 = 2 pi m/24` and `k_2 = 2 pi n/24` (`docs/GEOMETRY.md`). The uniform mode is `(0, 0)`,
with `mu = 6`. The minimum, `mu = -3`, is attained at `(8, 16)` and `(16, 8)`. Its real
eigenvectors are the three-sublattice patterns.

## Proposition 1: unique homogeneous equilibrium

*For every positive parameter set there is exactly one homogeneous equilibrium. It is positive
and has `D` in `(0, lam_D/d_D)`.*

At a homogeneous state, `s_D = z D` and `s_N = z N`. The N, A and S equations give, in order,

```
N(D) = lam_N / (d_N + (z f_N + k) D)            strictly decreasing
A(D) = z f_N N(D) D / d_A                        strictly increasing
S(D) = (lam_S/d_S) h_S(A(D))                     strictly decreasing
```

so equilibria are the roots of

```
F(D) = lam_D h_D(S(D)) - D [d_D + (z f_D + k) N(D)].
```

The first term is decreasing, because `h_D` increases with `S` and `S` decreases with `D`. The
second term, `d_D D + (z f_D + k) lam_N D / (d_N + (z f_N + k) D)`, is strictly increasing.
So `F` is strictly decreasing. Moreover `F(0) = lam_D h_D(lam_S/d_S) > 0` and
`F(lam_D/d_D) < lam_D (h_D - 1) < 0`. Hence `F` has exactly one root, and `N`, `A`, `S` follow
from it. `model/fixed_point.py` finds the root by Brent bracketing on `[0, lam_D/d_D]`.

The v1 code used a multistart `scipy.optimize.root` search and reported "no positive fixed
point" when every start failed. In the 4000-point global sample that happened at 733 points.
All 733 have an equilibrium (residual normalised by `max(1, ||x*||_inf)`). Corrected statistics are listed in
`docs/MANUSCRIPT_NUMBERS.md`.

## Proposition 2: the open loop is stable for every mode, so the static gain counts eigenvalue parity

Partition `J_mu` into the plant `z = (D, N, A)` and the scalar feedback state `S`:
`A_mu = J_mu[0:3, 0:3]`, `b = J_mu[0:3, 3]`, `c^T = J_mu[3, 0:3]`, `a_S = -d_S`. Define

```
L_mu(s) = c^T (sI - A_mu)^{-1} b / (s - a_S),      G(mu) = L_mu(0) = c^T A_mu^{-1} b / a_S.
```

**(a) `A_mu` is Hurwitz for every `mu` in `[-z, z]`.** `A_mu` is block lower-triangular. Its
eigenvalues are `-d_A` and those of the D-N block `B_mu = [[a11, a12], [a21, a22]]`. The trace
of `B_mu` is negative. Its determinant is

```
Delta(mu) = a11 a22 - D N (k + f_D mu)(k + f_N mu).
```

For `|mu| <= z`, we have `|k + f mu| <= k + z f`. Therefore

```
Delta(mu) >= a11 a22 - D N (k + z f_D)(k + z f_N)
          = d_D d_N + d_D (z f_N + k) D + d_N (z f_D + k) N > 0.
```

So `B_mu` is Hurwitz, and `det A_mu = -d_A Delta(mu) < 0`. In this model, then, the
cis/trans D-N subsystem with Delta production held fixed cannot by itself destabilise any mode.
Every instability passes through the S-dependent loop.

**(b) Static identity.** The Schur complement gives
`det J_mu = det(A_mu) (a_S - c^T A_mu^{-1} b) = det(A_mu) a_S [1 - G(mu)]`. This identity is
standard. With `det(A_mu) a_S > 0` from (a):

* `G(mu) = 1` if and only if `J_mu` has a zero eigenvalue.
* `G(mu) > 1` if and only if `det J_mu < 0`, which holds if and only if `J_mu` has an odd
  number of positive real eigenvalues. Complex pairs contribute `|lambda|^2 > 0` to the
  determinant, and the 4x4 block has an even number of real eigenvalues.
* `G(mu) < 1` is consistent with 0, 2 or 4 right-half-plane eigenvalues. **It is not a
  stability certificate.**

**(c) Complete per-mode test.** The open loop has no right-half-plane poles (from (a), and
`a_S < 0`). By the Nyquist criterion, the number of right-half-plane eigenvalues of `J_mu`
equals the number of clockwise encirclements of `+1` by `L_mu(i omega)`, `omega` in `R`,
provided no closed-loop eigenvalue lies on the imaginary axis.
`G(mu)` is the `omega = 0` point of this curve.

In the global sample, the static test and the full eigenvalue count disagree at 12 points.
Each of the 12 has exactly two right-half-plane eigenvalues in the `mu = -3` block: a complex
pair at 6 points and two real eigenvalues at the other 6. All 12 have `G(-3) < 1`. Nyquist
winding equals the eigenvalue count at all 4000 points (`scripts/global_sample.py`).

## Proposition 3: the static gain is largest at the most negative graph eigenvalue

*For every positive parameter set and every connected `z`-regular graph, `G(mu)` attains its
maximum over the graph spectrum only at `mu_min`, and `G(z) < 0`.*

Write `G(mu) = kappa_D T_mu kappa_S` with `kappa_S = c_S/d_S` and
`T_mu = [A_mu^{-1}]_{A,D}` (`control/state_space.py::stage_gains`). Expanding the cofactor,

```
T_mu = f_N N n(mu) / (d_A Delta(mu)),      n(mu) = z k D - mu (d_N + k D).
```

So `G(mu) = K n(mu)/Delta(mu)` with `K = kappa_D kappa_S f_N N / d_A > 0` and `Delta > 0`.

* `n` is linear and decreasing in `mu`, with `n(mu_min) > 0` (because `mu_min < 0`) and
  `n(z) = -z d_N < 0`. Its root `mu_0` lies in `(mu_min, z)`. For `mu >= mu_0`,
  `G(mu) <= 0 < G(mu_min)`. In particular `G(z) < 0`, so the uniform mode never has a
  stationary instability. This agrees with Proposition 1: a zero eigenvalue at `mu = z` would be
  a fold of homogeneous equilibria.
* For `mu` in `(mu_min, mu_0)`, let `h(mu) = n(mu_min) Delta(mu) - n(mu) Delta(mu_min)`, so that
  `G(mu) < G(mu_min)` if and only if `h(mu) > 0`. The function `h` is quadratic in `mu` with
  leading coefficient `-n(mu_min) D N f_D f_N < 0`, so it is concave. Also `h(mu_min) = 0` and
  `h(mu_0) = n(mu_min) Delta(mu_0) > 0`. A concave function lies above its chord, so
  `h > 0` on `(mu_min, mu_0]`.

Consequence: along any parameter path that starts at a stable homogeneous state, the first
**stationary** instability occurs in the `mu_min` eigenspace. On the triangular contact lattice
this is the three-sublattice mode. On a bipartite graph, where `mu_min = -z`, it is the
checkerboard. Non-stationary instabilities are not covered by this argument. A search of 40 000
further random parameter sets found none in which another nonuniform mode dominated or
was unstable with K stable (`scripts/mode_dominance_search.py`). This search at sampled points
does not track bifurcations along parameter paths. Subdominant nonuniform modes are also
unstable at 1451 of the 4000 global-sample points. At P2, 21 distinct nonuniform graph
eigenvalues have positive growth rates (Fig. S3). The uniform mode has an unstable complex
pair at 45 of the 4000 points. All 45 also have an unstable three-sublattice mode.

In a single-channel lateral-inhibition model (Collier et al. 1996; Arcak 2013), `G(mu)` is
proportional to `-mu`, so this monotonicity holds trivially. Here it is not automatic. Cis
inactivation (`k`) and trans-depletion of ligand (`f_D`) make the neighbour coupling two-channel,
and `G(mu)` becomes a ratio of a linear and a quadratic function of `mu`. The proposition shows
that the extra channels never move the selected mode at onset.

## Relation to graph input-output criteria

For nonsingular `J_loc`, write `J_mu = J_loc + mu B C`, where `B` carries the neighbour sums `(s_D, s_N)` into the cell
and `C` reads out `(D, N)` (`control/graph_io.py`). Sylvester's identity gives

```
det J_mu = det(J_loc) det(I - mu K),     K = -C J_loc^{-1} B   (2x2 static I/O gain).
```

With one coupling channel, `K` is a scalar and `mu_min K = 1` is the static-gain condition used
in graph-based lateral-inhibition analyses (Arcak 2013; Rufino Ferreira & Arcak 2013). The zero
sets of `det(I - mu K)` and `1 - G(mu)` coincide, and the relative difference in
`figures/data/graph_io_check.csv` is at most 1e-15. Both criteria test only the
zero-frequency loop. The `S`-partition used here is convenient because it factors the gain into
named biochemical stages. It does not detect anything the graph I/O form misses.
This algebraic identity does not establish the monotonicity or globally stable single-cell
response assumptions in the cited graph-based theorems. Those theorems also address existence
and stability of patterned equilibria under their assumptions; they are not restricted to
checking a determinant. Their normalised contact matrix has eigenvalues `mu/z` on a regular
graph, with a corresponding factor `z` in the input gain.

## Proposition 4: matched reductions share every equilibrium and the stationary boundary

The **DNA** reduction slaves `S` to `A` (`S = S_q(A) = (lam_S/d_S) h_S(A)`). The **DN**
reduction also slaves `A` (`A = f_N N s_D / d_A`) (`model/reduced.py`). Both impose
`dS/dt = 0` (and `dA/dt = 0`) cell by cell. Their equilibria, homogeneous or patterned, are
therefore exactly those of the full model. Their mode-resolved Jacobians are Schur complements
of `J_mu`:

```
det J_mu = a_S det J_mu^DNA = det([[-d_A, 0], [-c_S, -d_S]]) det J_mu^DN,
```

and the prefactors are nonzero. So for every `mu` the full and reduced models acquire a zero
eigenvalue at exactly the same parameters. They share the stationary boundary, the critical
`lam_N`, and the selected mode at onset. Their spectra are otherwise different. They differ in
growth rates and in instabilities that do not pass through a zero eigenvalue.
The shared zero-eigenvalue set is an actual stability boundary only if an oscillatory
instability has not already intervened.

In the global sample the DNA reduction is unstable at every point where the full model is
unstable. It is also unstable at 26 further points, all with `G(-3) < 1`. Where both are unstable,
its growth rate is larger at every point, by a median factor of 4.9. The DN reduction has 47
extra instabilities and a median factor of 216. Equality of equilibrium sets does not imply
equality of their stability or basins. The common final P2 pattern is a numerical observation
for the specified initial states, not a general nonlinear equivalence theorem.

## Proposition 5: the proneural timescale sets the growth rate, not the onset

The family `(lam_S, d_S) -> alpha (lam_S, d_S)` leaves `S_q(A)`, every equilibrium, every
`G(mu)`, the DNA reduction, and the stationary boundary unchanged (`model/reduced.py`). Only the
relaxation rate of `S` changes. Two limits follow:

* `alpha -> infinity`: `S` is fast and the spectrum of `J_mu` tends to that of `J_mu^DNA`,
  together with one eigenvalue near `-alpha d_S`.
* `alpha -> 0`: `S` is slowest. The plant relaxes to `z = -A_mu^{-1} b S`, so
  `dS/dt = (a_S - c^T A_mu^{-1} b) S = a_S [1 - G(mu)] S`. The growth rate is therefore
  `lambda_mu ~ alpha d_S^(0) [G(mu) - 1]`, where `d_S^(0)` denotes the unscaled rate.
  This is the slow eigenvalue for small `alpha`; it is not a general finite-timescale
  stability criterion when `G < 1`.

At P2 (`lam_N = 2`, `alpha = 1`): `G(-3) = 1.809` and `lambda = 0.103`. The slow-S estimate
gives 0.081, the DNA reduction gives 1.597 (15.5 times faster), and the DN reduction gives 11.28
(109 times). Across `alpha` in `[1e-2, 1e4]`, `lam_N^c = 1.304815342591` in every case, while the
P2 growth rate ranges from 8.1e-4 to 1.59 (`scripts/timescale_family.py`).

## Independent sampling robustness

`manuscript/robustness_supplement.md` reports a protocol fixed before 115000 new evaluations.
Ten original-box LHS replicates of 4000 points give a mean nonuniform-instability fraction
37.37% (approximate 95% Monte Carlo interval across replicates 37.07–37.67%). The pooled
DNA/full median growth ratio is 5.41. Sample-size checks and a wider box retain the direction
of the dynamic bias while changing its magnitude. The wider-box mean unstable fraction is
28.89%, with DNA replicate medians 2.99–3.64.

No reduction-missed full instability, slower reduced growth, faster competing nonuniform
mode or failure of the static-gain ordering was found in those samples. Dynamic ordering
remains an empirical result over the stated measures, not an additional theorem. The one
near-marginal full growth rate is reported explicitly and confirmed positive at 80-digit
precision. Analytic propositions above retain their original mathematical scope.
