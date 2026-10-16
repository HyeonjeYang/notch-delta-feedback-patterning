"""Recompute every statistic of the 4000-point global parameter sample.

Input: figures/data/global_sample_parameters.csv, the 13 kinetic parameters
of each point, unchanged from v1. They are a Latin hypercube sample
(scipy.stats.qmc.LatinHypercube, seed 20260817, n = 4000) over p0 x [0.1, 10]
for every rate and constant (log-uniform) and [1, 6] for the two Hill
exponents (uniform). ``latin_hypercube_parameters`` regenerates the sample, and
``main`` first checks that it matches the stored file exactly.
``v1_classification`` is carried along for traceability only.

Output: figures/data/global_sample_results.csv. One row per point:

* the unique homogeneous equilibrium (model.fixed_point) and its residual;
* the full mode scan: leading nonuniform growth rate, its mode, and the
  uniform-mode leading eigenvalue;
* the mu = -3 block: G(-3), the number of right-half-plane eigenvalues
  counted directly and by Nyquist winding, and whether any of them are
  complex;
* argmax over mu of G(mu);
* the same mode scan for the matched DNA and DN reductions (model/reduced.py).

Run time is about 1-2 minutes.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from tqdm import tqdm

from control.state_space import G_exact, nyquist_rhp_count, partition
from model.fixed_point import homogeneous_fixed_point, rhs_single_cell
from model.jacobian import jacobian_block
from model.params import NotchDeltaParams
from model.reduced import reduce_AS, reduce_S
from stability.mode_scan import DOMINANT_PATTERN_MU, UNIQUE_MU, lambda_pattern
from stability.sensitivity import PARAM_NAMES

DATA_DIR = Path(__file__).resolve().parents[1] / "figures" / "data"
LHS_SEED, LHS_N = 20260817, 4000
HILL = ("n_d", "n_ASC")


def latin_hypercube_parameters(seed: int = LHS_SEED, n: int = LHS_N) -> np.ndarray:
    """(n, 13) array in PARAM_NAMES order, exactly as drawn for v1."""
    from scipy.stats import qmc

    from model.params import P0
    unit = qmc.LatinHypercube(d=len(PARAM_NAMES), seed=seed).random(n=n)
    cols = []
    for i, name in enumerate(PARAM_NAMES):
        if name in HILL:
            cols.append(1.0 + unit[:, i] * 5.0)
        else:
            lo, hi = 0.1 * getattr(P0, name), 10.0 * getattr(P0, name)
            cols.append(np.exp(np.log(lo) + unit[:, i] * (np.log(hi) - np.log(lo))))
    return np.column_stack(cols)


def classify(lam_pattern: float, n_rhp: int, n_complex: int) -> str:
    if lam_pattern <= 0:
        return "stable"
    if n_rhp % 2 == 1:
        return "stationary"            # odd count: G(-3) > 1
    return "two_eig_complex" if n_complex else "two_eig_real"


def analyse(p: NotchDeltaParams) -> dict:
    fp = homogeneous_fixed_point(p)
    resid = float(np.max(np.abs(rhs_single_cell(fp, p))) / max(1.0, np.max(np.abs(fp))))
    full = lambda_pattern(fp, p)
    dna = lambda_pattern(fp, p, reduce=reduce_S)
    dn = lambda_pattern(fp, p, reduce=reduce_AS)
    J3 = jacobian_block(DOMINANT_PATTERN_MU, fp, p)
    ev3 = np.linalg.eigvals(J3)
    rhp = ev3.real > 0
    G_all = [G_exact(jacobian_block(mu, fp, p)) for mu in UNIQUE_MU]
    A_hurwitz = all(np.max(np.linalg.eigvals(partition(jacobian_block(mu, fp, p))[0]).real) < 0
                    for mu in UNIQUE_MU)
    return {
        "D_star": fp[0], "N_star": fp[1], "A_star": fp[2], "S_star": fp[3], "eq_residual": resid,
        "lambda_pattern": full["nonuniform_max_re_lambda"], "best_mu": full["nonuniform_best_mu"],
        "lambda_pattern_im": full["nonuniform_im_lambda"],
        "uniform_re": full["uniform_re_lambda"], "uniform_im": full["uniform_im_lambda"],
        "G_minus3": G_exact(J3), "G_argmax_mu": float(UNIQUE_MU[int(np.argmax(G_all))]),
        "A_mu_hurwitz_all_mu": A_hurwitz,
        "n_rhp_minus3": int(rhp.sum()), "n_rhp_minus3_complex": int((rhp & (np.abs(ev3.imag) > 1e-12)).sum()),
        "nyquist_rhp_minus3": nyquist_rhp_count(J3),
        "classification": classify(full["nonuniform_max_re_lambda"], int(rhp.sum()),
                                   int((rhp & (np.abs(ev3.imag) > 1e-12)).sum())),
        "lambda_DNA": dna["nonuniform_max_re_lambda"], "lambda_DNA_im": dna["nonuniform_im_lambda"],
        "best_mu_DNA": dna["nonuniform_best_mu"],
        "lambda_DN": dn["nonuniform_max_re_lambda"], "lambda_DN_im": dn["nonuniform_im_lambda"],
        "best_mu_DN": dn["nonuniform_best_mu"],
    }


def main() -> None:
    params = pd.read_csv(DATA_DIR / "global_sample_parameters.csv", float_precision="round_trip")
    regenerated = latin_hypercube_parameters()
    assert np.allclose(regenerated, params[PARAM_NAMES].values, rtol=1e-14, atol=0), "stored sample != seed"
    print(f"stored parameters match the Latin hypercube regenerated from seed {LHS_SEED}")
    rows = []
    for _, r in tqdm(params.iterrows(), total=len(params), desc="global sample"):
        p = NotchDeltaParams(*[float(r[name]) for name in PARAM_NAMES])
        rows.append({"sample_id": int(r["sample_id"]), "v1_classification": r["v1_classification"], **analyse(p)})
    out = pd.DataFrame(rows)
    out.to_csv(DATA_DIR / "global_sample_results.csv", index=False)

    n = len(out)
    unstable = out.lambda_pattern > 0
    pred = out.G_minus3 > 1
    print(f"points: {n};  max equilibrium residual {out.eq_residual.max():.1e}")
    print(f"v1 'no_positive_fixed_point' points now solved: {(out.v1_classification == 'no_positive_fixed_point').sum()}")
    print(f"nonuniform-unstable: {unstable.sum()}/{n} = {unstable.mean():.5f}")
    print(f"classification: {out.classification.value_counts().to_dict()}")
    print(f"dominant mode mu=-3 for every unstable point: {bool(np.allclose(out.best_mu[unstable], -3.0))}")
    print(f"argmax_mu G(mu) = -3 for every point: {bool(np.allclose(out.G_argmax_mu, -3.0))}")
    print(f"A_mu Hurwitz for every mu and point: {bool(out.A_mu_hurwitz_all_mu.all())}")
    print(f"max G(+6) sign check skipped; uniform-mode unstable: {(out.uniform_re > 0).sum()} "
          f"(oscillatory {((out.uniform_re > 0) & (out.uniform_im > 1e-9)).sum()}, "
          f"with K stable {((out.uniform_re > 0) & ~unstable).sum()})")
    print(f"G(-3)>1 vs instability: agree {(pred == unstable).sum()}/{n} ({(pred == unstable).mean():.4%}); "
          f"mismatches {(pred != unstable).sum()}, all with G<1 and 2 RHP eigenvalues: "
          f"{bool(((out.n_rhp_minus3 == 2) & ~pred)[pred != unstable].all())}")
    v1_valid = out.v1_classification != "no_positive_fixed_point"
    print(f"  on the v1-valid subset ({v1_valid.sum()}): mismatches {(pred != unstable)[v1_valid].sum()}")
    print(f"Nyquist winding == eigenvalue count at mu=-3: {(out.nyquist_rhp_minus3 == out.n_rhp_minus3).sum()}/{n}")
    for key in ("DNA", "DN"):
        u = out[f"lambda_{key}"] > 0
        both = unstable & u
        ratio = out[f"lambda_{key}"][both] / out.lambda_pattern[both]
        print(f"{key}: unstable {u.sum()}, disagree with full {(u != unstable).sum()} "
              f"(reduced-only {(u & ~unstable).sum()}, full-only {(~u & unstable).sum()}); "
              f"growth ratio reduced/full median {ratio.median():.2f} [IQR {ratio.quantile(.25):.2f}-{ratio.quantile(.75):.2f}]")
        stat_full = out.classification == "stationary"
        print(f"   stationary-full points also unstable in {key}: {(stat_full & u).sum()}/{stat_full.sum()}")
    print(f"wrote {DATA_DIR / 'global_sample_results.csv'}")


if __name__ == "__main__":
    main()
