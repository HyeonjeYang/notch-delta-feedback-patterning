# Geometry: six-neighbour contact lattice and mode selection

Condensed from the source repository's `preprint_final_analysis/GEOMETRY_DERIVATION.md`.

## Terminology

A **periodic six-neighbour triangular contact lattice, corresponding to a
hexagonal cellular tiling** -- not a "honeycomb graph" (that would be
degree-3 and bipartite; this contact graph is degree-6 and non-bipartite,
the triangular lattice's dual).

## Dispersion relation

`geometry/lattice.py::hex_adj`'s six neighbour offsets --
`(-1,0), (-1,1), (0,-1), (0,1), (1,0), (1,-1)` -- give, for a plane wave on
the periodic `rows x cols` grid (`k1=2*pi*m/rows, k2=2*pi*n/cols`):

```
mu(k1, k2) = 2*cos(k1) + 2*cos(k2) + 2*cos(k1 - k2)
```

Verified numerically against the explicit adjacency matrix's eigenvalues to
**2.8e-14** (`geometry/graph_spectrum.py::verify_dispersion`).

```
mu(Gamma, k=0)                    = 6.000000   (uniform mode, graph maximum)
mu(K-like, m=rows/3, n=2*rows/3)  = -3.000000  (graph minimum)
```

## Why K-like points give mu=-3, exactly

At `k1=2*pi/3, k2=4*pi/3`, all three cosine terms equal exactly `-1/2`:
`mu = 2*(-1/2)*3 = -3`. Exactly two `(m,n)` pairs on the 24x24 grid achieve
this: `(8,16)` and `(16,8)`, both multiples of `rows/3=8` -- i.e. **exactly
three-fold periodic**. This connects the graph-theoretic minimum eigenvalue
directly to a real-space three-sublattice pattern: the natural
generalization of a two-colour checkerboard to a non-bipartite lattice
where no two-colouring can make every edge alternate.

## Reciprocal-space growth map

At the validated P2 point (`lam_N=2.0`, `geometry/reciprocal_space.py`):

```
lambda_max(k=0)                  = -0.100000   (homogeneous perturbation stable)
lambda_max(K-like, (m,n)=(8,16)) =  0.103067   (strongest positive growth)
```

confirming the tissue behaves as a mode-selective spatial filter.
