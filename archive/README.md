# Archive

Frozen earlier states of this repository. Nothing in here is imported, tested,
or used to build the current figures. Files are kept unchanged so that every
number in an earlier manuscript version can still be traced to its source.

| folder | state | git reference |
|---|---|---|
| `v1_physbio_submission/` | the full repository as it stood when the first manuscript version ("Mode-resolved feedback gain controls spatial symmetry breaking in an explicit Notch-Delta-proneural circuit") was submitted | commit `8503658`, tag `v1-physbio-submission` |

## Known problems in v1 (corrected in the current version)

These are the issues raised in review. Each one is fixed in the active tree, not
in this archive.

1. **Equilibrium solver failures counted as nonexistence.**
   `model/fixed_point.py::find_fixed_points` used a multistart `scipy.optimize.root`
   and treated a failure to converge as "no positive fixed point". In the
   4000-point global sample, 733 points were labelled this way. The
   homogeneous equilibrium in fact exists and is unique for every
   positive parameter set. All 733 have one, and the nonuniform-instability
   fraction is 1453/4000 = 36.325%, not 1212/4000 = 30.3%. The 2-D phase maps
   (`figures/data/two_dimensional_phase_maps.csv`, 115 grid points) had the same gaps.
2. **Explanation of the DC-gain mismatches.** The v1 text said that the
   mismatched samples were "plausibly" cases in which another spatial mode dominated.
   In all of them the dominant mode is mu = -3. Each mismatch has exactly two
   right-half-plane eigenvalues in the mu = -3 block, either as a complex pair or as two real
   eigenvalues. With an even number of such eigenvalues, det J_mu > 0, so the DC gain stays below 1.
3. **Growth-rate validation initial condition.**
   `simulation/perturbation.py::perturbed_field` multiplied each component of the
   eigenvector by that component's equilibrium value. The resulting initial state is no longer an
   eigenvector of J_mu, so several modes are excited at t = 0 and an early-time transient
   follows. The early-window slopes in `growth_rate_validation.csv`, for example 0.44 at P2 where 0.103 is
   predicted, show this transient.
4. **`stochastic_check/`.** This was exploratory and is not part of the revised
   manuscript. It now lives only here.
