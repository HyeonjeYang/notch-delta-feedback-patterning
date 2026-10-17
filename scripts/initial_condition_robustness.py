"""Paired deterministic lattice experiments with broad-spectrum random starts.

Write the protocol before computation with --protocol-only. Runs are resumable
with --reuse, which verifies the numerical kernel and each stored case hash.
"""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.sparse import bmat, eye, lil_matrix

from geometry.graph_spectrum import circulant_eigenvalue
from geometry.lattice import COLS, ROWS, NEIGHBOR_OFFSETS
from model.equations import hill_ASC
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import P0
from model.reduced import reduce_S
from simulation.nonlinear_ode import pack, rhs_lattice, rhs_lattice_DNA

OUT = ROOT / "figures/data/initial_conditions"
SEEDS = tuple(range(2026101901, 2026101909))
AMPLITUDES = (1e-6, 1e-4, 1e-2)
POINTS = {"P2": (2.0, 600.0), "P3": (5.623, 300.0)}
RTOL, ATOL, RECORD = 1e-9, 1e-11, .5
LINEAR_LIMIT = 1e-3


def protocol():
    return {
        "version": 2, "seeds": list(SEEDS), "amplitudes": list(AMPLITUDES),
        "points": {k: {"lam_N": v[0], "t_final": v[1]} for k, v in POINTS.items()},
        "lattice": [ROWS, COLS], "generator": "Generator(PCG64(seed))",
        "random_fields": "four independent U[-1,1] fields; subtract each spatial mean; divide each by its max absolute value",
        "main_pair": "same D,N,A fields x* (1+amplitude*xi); full S=S_q(A), DNA reconstructs that same S",
        "downstream_control": "full S=S* (1+amplitude*xi_S) at amplitude 1e-4, same D,N,A as main pair",
        "coupling": "same seed reuses fields across points and amplitudes; designs are paired, not independent samples",
        "total_primary_runs": 112,
        "solver": {"method": "BDF with graph Jacobian sparsity", "rtol": RTOL, "atol": ATOL, "record_interval": RECORD},
        "onset": "first upward continuous event std(D)/mean(D)=0.1; also record threshold 0.5; no crossing is censored at t_final",
        "linear_validation": "exact all-mode matrix propagator of each model, initial prefix with max |x/x*-1| <= 1e-3; no slope fit",
        "endpoint": "fixed horizon, not automatically a converged attractor",
        "endpoint_metrics": ["Delta CV", "K power / nonuniform Delta power", "max |dx/dt|/x*", "last-10-time-unit Delta change / D*"],
        "stationary_candidate_rule": "endpoint max |dx/dt|/x* < 1e-7, descriptive numerical criterion only",
        "pair_metrics": "onset-time ratio; endpoint Delta correlation maximised over periodic translations, no rotations; min relative L2 field distance under translations",
        "tight_checks": "first seed, amplitude 1e-6, both models at both points; BDF rtol=1e-11, atol=1e-14, recording=0.25, same fixed horizons",
        "solver_check": "same four cases, DOP853 rtol=1e-11, atol=1e-14 until baseline CV=0.5 event plus 5 time units; compare onset events",
        "protocol_history": "version 1 DOP853 pilot: one full P2 run completed; stiff DNA run interrupted for cost before a result was written. Solver changed to sparse BDF before production; design, thresholds and fixed horizons retained",
        "scope": "two operating points, bounded perturbations near the homogeneous state; no stochastic forcing, tissue fit, global basin survey or global stability proof",
    }


def kernel_hash():
    paths = [Path(__file__).resolve()]
    for directory in ("model", "simulation", "geometry"):
        paths.extend((ROOT / directory).glob("*.py"))
    data = "".join(p.relative_to(ROOT).as_posix() + "\0" + p.read_text(encoding="utf-8") + "\0"
                   for p in sorted(paths))
    return hashlib.sha256(data.encode()).hexdigest()


def random_fields(seed, rows=ROWS, cols=COLS):
    xi = np.random.Generator(np.random.PCG64(seed)).uniform(-1, 1, (4, rows, cols))
    xi -= xi.mean(axis=(1, 2), keepdims=True)
    return xi / np.max(np.abs(xi), axis=(1, 2), keepdims=True)


def initial_fields(fp, p, seed, amplitude, downstream="qss"):
    fields = fp[:, None, None] * (1 + amplitude * random_fields(seed))
    if downstream == "qss":
        fields[3] = p.lam_ASC / p.d_ASC * hill_ASC(fields[2], p)
    elif downstream != "free":
        raise ValueError(downstream)
    if not np.all(fields > 0):
        raise RuntimeError("nonpositive initial concentration")
    return fields


def spatial_metrics(fields):
    """Last two dimensions are the periodic lattice; allows a time batch."""
    mean = fields.mean(axis=(-2, -1))
    cv = fields.std(axis=(-2, -1)) / mean
    power = np.abs(np.fft.fft2(fields - mean[..., None, None])) ** 2
    denominator = power.sum(axis=(-2, -1))
    k = power[..., ROWS//3, 2*COLS//3] + power[..., 2*ROWS//3, COLS//3]
    purity = np.divide(k, denominator, out=np.zeros_like(k), where=denominator > 0)
    return cv, purity


def translation_comparison(a, b):
    """Direct circular cross correlation, then measure aligned field distance."""
    da, db = a-a.mean(), b-b.mean()
    cross = np.fft.ifft2(np.fft.fft2(da) * np.conj(np.fft.fft2(db))).real
    shift = np.unravel_index(np.argmax(cross), cross.shape)
    aligned = np.roll(b, shift, axis=(0, 1))
    denom = np.linalg.norm(da)*np.linalg.norm(db)
    corr = float(cross[shift]/denom) if denom > 0 else float("nan")
    distance = float(np.linalg.norm(a-aligned)/np.linalg.norm(a))
    return corr, distance, [int(s) for s in shift]


def linear_delta(initial, fp, p, model, times):
    """Independent spatial Fourier blocks, retaining all biochemical poles."""
    nsp = len(initial)
    j = np.stack([jacobian_block(circulant_eigenvalue(m,n), fp, p)
                  for m in range(ROWS) for n in range(COLS)])
    if model == "DNA":
        j = np.stack([reduce_S(b) for b in j])
    values, vectors = np.linalg.eig(j)
    f0 = np.fft.fft2(initial - fp[:nsp, None, None]).reshape(nsp, -1).T
    weights = np.linalg.solve(vectors, f0[..., None])[..., 0]
    fD = np.einsum("kj,tkj,kj->tk", vectors[:, 0, :],
                   np.exp(np.asarray(times)[:,None,None]*values[None,:,:]), weights)
    return np.fft.ifft2(fD.reshape(-1, ROWS, COLS)).real


def specifications():
    for point in POINTS:
        for amplitude in AMPLITUDES:
            for seed in SEEDS:
                for model in ("full", "DNA"):
                    yield point, amplitude, seed, model, "qss"
        for seed in SEEDS:
            yield point, 1e-4, seed, "full", "free"


def jacobian_sparsity(model):
    cells = ROWS*COLS
    W = lil_matrix((cells,cells),dtype=bool)
    for r in range(ROWS):
        for c in range(COLS):
            for dr,dc in NEIGHBOR_OFFSETS:
                W[r*COLS+c, ((r+dr)%ROWS)*COLS+(c+dc)%COLS] = True
    W = W.tocsr()
    I = eye(cells,format="csr",dtype=bool)
    blocks = [[I,W+I,I],[W+I,I,None],[W,I,I]] if model == "DNA" else [
        [I,W+I,None,I],[W+I,I,None,None],[W,I,I,None],[None,None,I,I]]
    return bmat(blocks,format="csr")


def run_case(point, amplitude, seed, model, downstream, tight=False, method="BDF", end_override=None):
    p = dataclasses.replace(P0, lam_N=POINTS[point][0])
    fp = homogeneous_fixed_point(p)
    ns = 4 if model == "full" else 3
    initial = initial_fields(fp, p, seed, amplitude, downstream)[:ns]
    rhs = rhs_lattice if model == "full" else rhs_lattice_DNA
    end = POINTS[point][1] if end_override is None else end_override
    record = .25 if tight else RECORD
    times = np.arange(0, end+record/2, record)

    def event10(t, y, p):
        d = y[:ROWS*COLS]
        return d.std()/d.mean() - .1

    def event50(t, y, p):
        d = y[:ROWS*COLS]
        return d.std()/d.mean() - .5

    event10.direction = event50.direction = 1
    start = time.perf_counter()
    options = {"jac_sparsity":jacobian_sparsity(model)} if method == "BDF" else {}
    sol = solve_ivp(rhs, (0, end), pack(*initial), args=(p,), method=method,
                    rtol=1e-11 if tight else RTOL, atol=1e-14 if tight else ATOL,
                    t_eval=times, events=(event10, event50), **options)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(f"{point}/{seed}/{model}: {sol.message}")
    fields = sol.y.T.reshape(len(times), ns, ROWS, COLS)
    cv, purity = spatial_metrics(fields[:,0])
    deviation = np.max(np.abs(fields/fp[:ns][None,:,None,None]-1),axis=(1,2,3))
    violations = np.flatnonzero(deviation > LINEAR_LIMIT)
    prefix = int(violations[0]) if len(violations) else len(times)
    linear_cv = np.full(len(times), np.nan)
    maxerr = float("nan")
    if prefix >= 2:
        pred = linear_delta(initial, fp, p, model, times[:prefix])
        actual = fields[:prefix,0]-fp[0]
        error = np.linalg.norm((actual-pred).reshape(prefix,-1),axis=1)
        norm = np.linalg.norm(actual.reshape(prefix,-1),axis=1)
        maxerr = float(np.max(error/np.maximum(norm,1e-30)))
        linear_cv[:prefix] = pred.std(axis=(1,2))/(fp[0]+pred.mean(axis=(1,2)))
    residual = float(np.max(np.abs(rhs(end,sol.y[:,-1],p).reshape(ns,ROWS,COLS))/fp[:ns,None,None]))
    late = int(round(10/record))
    change = float(np.max(np.abs(fields[-1,0]-fields[-1-late,0]))/fp[0])
    summary = dict(point=point, amplitude=amplitude, seed=seed, model=model, downstream=downstream,
                   tight=tight, t_final=end, onset_01=float(sol.t_events[0][0]) if len(sol.t_events[0]) else None,
                   onset_05=float(sol.t_events[1][0]) if len(sol.t_events[1]) else None,
                   final_cv=float(cv[-1]), final_k_fraction=float(purity[-1]),
                   final_scaled_rhs=residual, last10_D_change=change,
                   stationary_candidate=residual < 1e-7, min_concentration=float(sol.y.min()),
                   linear_n=prefix, linear_end=float(times[prefix-1]) if prefix else None,
                   linear_max_relative_D_error=maxerr if np.isfinite(maxerr) else None,
                   elapsed_seconds=time.perf_counter()-start, nfev=sol.nfev)
    return summary, dict(times=times, cv=cv, k_fraction=purity, linear_cv=linear_cv,
                         initial=initial, final=fields[-1], fp=fp)


def case_name(spec, tight=False):
    point, amp, seed, model, down = spec
    return f"{point}_a{amp:.0e}_{seed}_{model}_{down}" + ("_tight" if tight else "")


def compute(spec, digest, reuse, tight=False):
    name = case_name(spec,tight)
    path = OUT / "cases" / f"{name}.npz"
    meta = path.with_suffix(".json")
    if reuse and meta.exists() and path.exists():
        stored = json.loads(meta.read_text())
        if stored["kernel_sha256"] != digest or stored["npz_sha256"] != hashlib.sha256(path.read_bytes()).hexdigest():
            raise RuntimeError(f"stale or modified case: {name}")
        return stored["summary"]
    summary, arrays = run_case(*spec,tight=tight)
    np.savez_compressed(path, **arrays)
    meta.write_text(json.dumps(dict(kernel_sha256=digest, npz_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                   summary=summary), indent=2)+"\n",encoding="utf-8")
    print(f"{name}: onset={summary['onset_01']}, K={summary['final_k_fraction']:.3f}, "
          f"rhs={summary['final_scaled_rhs']:.2g} ({summary['elapsed_seconds']:.1f}s)",flush=True)
    return summary


def make_pairs(frame):
    pairs = []
    main = frame[frame.downstream == "qss"]
    for (point, amp, seed), group in main.groupby(["point","amplitude","seed"]):
        a = group[group.model == "full"].iloc[0]
        b = group[group.model == "DNA"].iloc[0]
        fa = np.load(OUT/"cases"/(case_name((point,amp,seed,"full","qss"))+".npz"))["final"][0]
        fb = np.load(OUT/"cases"/(case_name((point,amp,seed,"DNA","qss"))+".npz"))["final"][0]
        corr, distance, shift = translation_comparison(fa,fb)
        pairs.append(dict(point=point,amplitude=amp,seed=seed,full_onset_01=a.onset_01,DNA_onset_01=b.onset_01,
                          onset_ratio=a.onset_01/b.onset_01,full_onset_05=a.onset_05,DNA_onset_05=b.onset_05,
                          aligned_correlation=corr, aligned_relative_L2=distance,shift_r=shift[0],shift_c=shift[1],
                          full_k_fraction=a.final_k_fraction,DNA_k_fraction=b.final_k_fraction,
                          full_stationary=a.stationary_candidate,DNA_stationary=b.stationary_candidate))
    return pd.DataFrame(pairs)


def aggregate(frame, pairs, tight):
    groups = []
    for (point,amp), g in pairs.groupby(["point","amplitude"]):
        groups.append(dict(point=point,amplitude=amp,n=len(g),
            **{key: {"median":float(g[key].median()),"min":float(g[key].min()),"max":float(g[key].max())}
               for key in ("full_onset_01","DNA_onset_01","onset_ratio","aligned_correlation","aligned_relative_L2","full_k_fraction","DNA_k_fraction")}))
    checks = []
    for row in tight:
        base = frame[(frame.point==row["point"]) & (frame.amplitude==row["amplitude"]) &
                     (frame.seed==row["seed"]) & (frame.model==row["model"]) & (frame.downstream=="qss")].iloc[0]
        checks.append(dict(point=row["point"],model=row["model"],
                           onset_abs_difference=abs(row["onset_01"]-base.onset_01),
                           k_fraction_abs_difference=abs(row["final_k_fraction"]-base.final_k_fraction)))
    return dict(total_primary_runs=len(frame), total_tight_checks=len(tight),groups=groups,
                n_onset_01_censored=int(frame.onset_01.isna().sum()),
                n_onset_05_censored=int(frame.onset_05.isna().sum()),
                n_stationary_candidates=int(frame.stationary_candidate.sum()),
                n_full_later=int((pairs.full_onset_01>pairs.DNA_onset_01).sum()), n_main_pairs=len(pairs),
                min_concentration=float(frame.min_concentration.min()),
                max_linear_relative_D_error=float(frame.linear_max_relative_D_error.max()),
                tight_checks=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol-only",action="store_true")
    parser.add_argument("--reuse",action="store_true")
    parser.add_argument("--workers",type=int,default=4)
    args = parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    target = OUT/"protocol.json"
    current = protocol()
    if target.exists() and json.loads(target.read_text()) != current:
        raise RuntimeError("protocol differs from the frozen protocol; use a new study directory")
    if not target.exists():
        target.write_text(json.dumps(current,indent=2)+"\n",encoding="utf-8")
    if args.protocol_only:
        print(target)
        return
    (OUT/"cases").mkdir(exist_ok=True)
    digest = kernel_hash()
    if args.workers < 1:
        parser.error("--workers must be positive")
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(partial(compute,digest=digest,reuse=args.reuse), specifications()))
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT/"runs.csv",index=False)
    pairs = make_pairs(frame)
    pairs.to_csv(OUT/"pairs.csv",index=False)
    tight_specs = [(point,1e-6,SEEDS[0],model,"qss") for point in POINTS for model in ("full","DNA")]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        tight = list(pool.map(partial(compute,digest=digest,reuse=args.reuse,tight=True),tight_specs))
    summary = aggregate(frame,pairs,tight)
    solver_checks = []
    for spec in tight_specs:
        base = frame[(frame.point==spec[0]) & (frame.amplitude==spec[1]) & (frame.seed==spec[2]) &
                     (frame.model==spec[3]) & (frame.downstream=="qss")].iloc[0]
        end = float(np.ceil((base.onset_05+5)*4)/4)
        check, arrays = run_case(*spec,tight=True,method="DOP853",end_override=end)
        np.savez_compressed(OUT/"cases"/(case_name(spec)+"_DOP853.npz"),**arrays)
        solver_checks.append(dict(point=spec[0],model=spec[3],t_final=end,
            onset_01_difference=abs(check["onset_01"]-base.onset_01),
            onset_05_difference=abs(check["onset_05"]-base.onset_05),summary=check))
        print(f"cross-solver check {spec[0]} {spec[3]} complete",flush=True)
    summary["solver_checks"] = solver_checks
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(f"completed {len(frame)} primary runs and {len(tight)} tight checks",flush=True)


if __name__ == "__main__":
    main()
