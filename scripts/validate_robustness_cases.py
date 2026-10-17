"""Validate representative exceptions and marginal cases at 80-digit precision."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from model.fixed_point import homogeneous_fixed_point
from model.params import NotchDeltaParams
from model.reduced import reduce_AS, reduce_S
from stability.mode_scan import lambda_pattern
from stability.sensitivity import PARAM_NAMES
from scripts.robustness_study import kernel_digest, protocol

OUT = ROOT / "figures/data/robustness"


def load_samples():
    frames = []
    current_kernel = kernel_digest()
    for path in sorted(OUT.glob("*.csv.gz")):
        meta = json.loads(path.with_name(path.name.replace(".csv.gz", ".json")).read_text())
        if meta["kernel_sha256"] != current_kernel:
            raise RuntimeError(f"replicate was produced by a different computational kernel: {path.name}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != meta["csv_sha256"]:
            raise RuntimeError(f"corrupt replicate: {path.name}")
        df = pd.read_csv(path, float_precision="round_trip")
        for key in ("scenario", "seed", "n"):
            df[key] = meta[key]
        frames.append(df)
    if not frames:
        raise RuntimeError("Run scripts/robustness_study.py first")
    out = pd.concat(frames, ignore_index=True)
    if len(out) != protocol()["total_new_points"] or len(frames) != 30:
        raise RuntimeError("Robustness study is incomplete")
    return out


def high_precision(row):
    p_double = NotchDeltaParams(**{k: float(row[k]) for k in PARAM_NAMES})
    fp_double = homogeneous_fixed_point(p_double)
    mus = {"full": float(row.best_mu),
           "DNA": lambda_pattern(fp_double, p_double, reduce=reduce_S)["nonuniform_best_mu"],
           "DN": lambda_pattern(fp_double, p_double, reduce=reduce_AS)["nonuniform_best_mu"]}
    with mp.workdps(80):
        p = {k: mp.mpf(float(row[k])) for k in PARAM_NAMES}
        def chain(D):
            N = p["lam_N"] / (p["d_N"] + (6 * p["f_N"] + p["k_cis"]) * D)
            A = 6 * p["f_N"] * N * D / p["d_A"]
            S = p["lam_ASC"] / p["d_ASC"] / (1 + (A / p["Kd"]) ** p["n_ASC"])
            return N, A, S
        def F(D):
            N, A, S = chain(D)
            return p["lam_D"] / (1 + (p["Kd"] / S) ** p["n_d"]) - D * (p["d_D"] + (6 * p["f_D"] + p["k_cis"]) * N)
        lo, hi = mp.mpf(0), p["lam_D"] / p["d_D"]
        for _ in range(2000):
            mid = (lo + hi) / 2
            if F(mid) > 0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 8 * mp.eps * max(abs(mid), mp.mpf("1e-200")):
                break
        else:
            raise RuntimeError("high-precision scalar bracketing did not converge")
        D = (lo + hi) / 2
        N, A, S = chain(D)
        u, v = (p["Kd"] / S) ** p["n_d"], (A / p["Kd"]) ** p["n_ASC"]
        kd = p["lam_D"] * p["n_d"] * u / S / (1 + u) ** 2
        cs = -p["lam_ASC"] * p["n_ASC"] * v / A / (1 + v) ** 2
        def block(mu):
            mu = mp.mpf(str(mu))
            return mp.matrix([
                [-p["d_D"] - (6*p["f_D"]+p["k_cis"])*N, -D*(p["k_cis"]+p["f_D"]*mu), 0, kd],
                [-N*(p["k_cis"]+p["f_N"]*mu), -p["d_N"]-(6*p["f_N"]+p["k_cis"])*D, 0, 0],
                [p["f_N"]*N*mu, 6*p["f_N"]*D, -p["d_A"], 0],
                [0, 0, cs, -p["d_ASC"]]])
        def spectrum(J):
            return list(mp.eig(J, left=False, right=False))
        out = {}
        for name, mu in mus.items():
            J = block(mu)
            if name == "DNA":
                J = J[:3, :3] - J[:3, 3:4] * J[3:4, :3] / J[3, 3]
            elif name == "DN":
                J = J[:2, :2] - J[:2, 2:4] * J[2:4, 2:4] ** -1 * J[2:4, :2]
            eig = spectrum(J)
            lead = max(eig, key=mp.re)
            out[f"lambda_{name}_80digits"] = mp.nstr(mp.re(lead), 30)
            out[f"imag_{name}_80digits"] = mp.nstr(abs(mp.im(lead)), 30)
            out[f"lambda_{name}_80"] = float(mp.re(lead))
        Jk = block(-3)
        eigk = spectrum(Jk)
        out["n_rhp_K_80"] = sum(mp.re(e) > 0 for e in eigk)
        plant = Jk[:3, :3]
        gain = (Jk[3:4, :3] * plant ** -1 * Jk[:3, 3:4])[0] / Jk[3, 3]
        out["G_minus3_80digits"] = mp.nstr(gain, 30)
        original = {"full": row.lambda_pattern, "DNA": row.lambda_DNA, "DN": row.lambda_DN}
        out["signs_agree"] = all((out[f"lambda_{k}_80"] > 0) == (value > 0) for k, value in original.items())
        out["rhp_count_agrees"] = out["n_rhp_K_80"] == row.n_rhp_K
        out["max_abs_growth_difference"] = max(abs(out[f"lambda_{k}_80"] - value) for k, value in original.items())
        out["relative_D_equilibrium_difference"] = float(abs(D - mp.mpf(float(row.D_star))) / D)
        return out


def select_cases(df):
    selected = []
    mismatch = (df.lambda_pattern > 0) != (df.G_minus3 > 1)
    # Explanatory examples: largest growth within each already-defined pair class.
    for name, mask in (("static_miss_complex", mismatch & (df.n_complex_rhp_K > 0)),
                       ("static_miss_real", mismatch & (df.n_complex_rhp_K == 0))):
        for i in df[mask].sort_values(["lambda_pattern", "sample_id"], ascending=[False, True]).head(5).index:
            selected.append((name, i))
    marginal = (df[["lambda_pattern", "lambda_DNA", "lambda_DN"]].abs() <= 1e-8).any(axis=1)
    selected.extend(("near_marginal", i) for i in df[marginal].index)
    for name in ("DNA", "DN"):
        missed = (df.lambda_pattern > 1e-8) & (df[f"lambda_{name}"] < -1e-8)
        slower = (df.lambda_pattern > 1e-8) & (df[f"lambda_{name}"] > 1e-8) & (df.lambda_pattern - df[f"lambda_{name}"] > 1e-9)
        for label, mask, margin in ((f"{name}_missed", missed, np.minimum(df.lambda_pattern, -df[f"lambda_{name}"])),
                                    (f"{name}_slower", slower, df.lambda_pattern - df[f"lambda_{name}"])):
            # An empty frame must receive an empty Series: assigning the full
            # margin Series would expand it back to unrelated rows in pandas.
            for i in df[mask].assign(margin=margin.loc[mask]).sort_values(["margin", "sample_id"], ascending=[False, True]).head(5).index:
                selected.append((label, i))
    other = (df.lambda_pattern > 1e-8) & (df.dominance_excess > 1e-9)
    selected.extend(("other_mode_dominant", i) for i in df[other].nlargest(5, "dominance_excess").index)
    return selected


def main():
    df = load_samples()
    selected = select_cases(df)
    records = []
    cache = {}
    for kind, i in selected:
        row = df.loc[i]
        if i not in cache:
            cache[i] = high_precision(row)
        records.append({"case_class": kind, **{k: row[k] for k in ["scenario", "n", "seed", "sample_id", *PARAM_NAMES]},
                        "lambda_full_double": row.lambda_pattern, "lambda_DNA_double": row.lambda_DNA,
                        "lambda_DN_double": row.lambda_DN, **cache[i]})
        print(f"{kind}: {row.scenario}, n={row.n}, seed={row.seed}, id={row.sample_id}: "
              f"signs agree {cache[i]['signs_agree']}", flush=True)
    out = pd.DataFrame(records)
    out.to_csv(OUT / "high_precision_checks.csv", index=False)
    if not out.signs_agree.all() or not out.rhp_count_agrees.all():
        raise RuntimeError("High-precision classification differs; investigate before publication")
    print(f"80-digit checks passed for {len(cache)} distinct selected cases")


if __name__ == "__main__":
    main()
