# Manuscript numbers

Every value below was verified directly against the frozen source
analysis (`Neuroblast_differentiation_simulation/preprint_final_analysis/`)
before being copied here; none is restated from memory. Baseline `p0` =
`model.params.P0` (`LITERATURE_PARAMS` with `k_cis=0.35`); continuation
parameter `lam_N` (baseline 1.0, P2 2.0, P3 5.623), all else at `p0`.

| quantity | value |
|---|---|
| critical `lam_N` (mu=-3 leading pole crosses `Re(s)=0`) | 1.304815343 |
| `G_exact(-3)` at critical `lam_N` | 1.000000000000 |
| `G_exact(+6)` at critical `lam_N` | -0.0628 |
| `G_eff` at critical `lam_N` | ~3.84 |
| max `Re(s)` for mu=6 across `lam_N in [0.5,6.5]` | -0.100000 (always stable) |
| Schur vs. direct eigenvalue max coefficient diff (9 combinations) | 1.421e-13 |
| Schur vs. direct eigenvalue max root mismatch | 4.691e-14 |
| dispersion formula vs. adjacency-matrix spectrum, max deviation | 2.842e-14 |
| dominant unstable graph mode | (m,n)=(8,16) [and (16,8)], mu=-3 |
| `lambda_max(k=0)` at P2 | -0.100000 |
| `lambda_max` at K-like point, P2 | 0.103067 |
| P2 predicted `Re(lambda)` | 0.103067 |
| P2 observed growth (late window) | 0.0898-0.1029 |
| P3 predicted `Re(lambda)` | 0.216300 |
| P3 observed growth (late window, cleanest condition) | 0.233890 |
| P2 final `q_amp` (adaptive DOP853 / Euler dt=0.01 / dt=0.005) | 0.989569 (all three) |
| P3 final `q_amp` (adaptive / dt=0.01 / dt=0.005) | 0.999998 (all three) |
| Euler dt=0.01 clip fraction, P2 and P3 | 0 |
| final Moran's I, P2 and P3 | -0.5 (exact -- amplitude-blind by construction) |
| `G_exact(-3)>1` classification accuracy (frozen global sample, n=3267) | 99.79% |
| `G_eff` best-threshold classification accuracy | 93.88% (threshold 1.095) |
| nonuniform-unstable fraction, 4000-point global sample | 30.3% (1212/4000) |
| multistable branches found (entire source analysis) | 0 |

The three critical-`lam_N` rows above come from direct `brentq` root-finding
(bracket `[1.2, 1.4]`, `xtol=1e-13`), not grid interpolation; pole crossing
and unity-gain agree to `4.7e-15` in `lam_N` -- see `scripts/verify_critical_point.py`.

## Solver convergence (current, validated)

Adaptive DOP853 (`rtol=1e-9, atol=1e-11`) and Euler `dt=0.01`/`dt=0.005`
agree to 6 significant figures on final `q_amp` at both P2 and P3, with
zero Euler clipping. This is the current, valid convergence result for the
patterning regime studied in this repository. It supersedes an earlier,
invalidated `dt=0.02` finding at a *different* (`k_cis`-only) baseline --
see `README.md` "Reproducibility note". That invalidated result is not
reproduced here.
