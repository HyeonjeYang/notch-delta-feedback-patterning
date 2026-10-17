"""Finite-size, raw open-boundary and quenched-production controls."""
from __future__ import annotations
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import argparse,dataclasses,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from functools import partial
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from model.params import P0
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.reduced import reduce_S
from control.state_space import G_exact
from geometry.graph_spectrum import circulant_eigenvalue
from scripts.reviewer_models import adjacency,lattice_rhs,lattice_jacobian,random_initial
from scripts.reviewer_protocol import OUT,SEEDS,freeze,digest,sha

POINTS={"P0":(1.,600.),"P2":(2.,600.),"P3":(5.623,300.)}


def specifications():
    for point in ("P2","P3"):
        for n in (23,24,25):
            for seed in SEEDS:
                for model in ("full","DNA"):yield "periodic",point,n,seed,model,0.
    for point in POINTS:
        for model in ("full","DNA"):yield "open",point,24,SEEDS[0],model,0.
    for model in ("full","DNA"):yield "periodic", "P0",24,SEEDS[0],model,0.
    for point in ("P2","P3"):
        for cv in (.05,.1):
            for model in ("full","DNA"):yield "heterogeneous",point,24,SEEDS[0],model,cv


def name(spec):
    condition,point,n,seed,model,cv=spec
    return f"{condition}_{point}_{n}_{seed}_{model}_cv{cv}"


def run(spec,kernel,reuse=False):
    condition,point,n,seed,model,cv=spec
    path=OUT/"geometry"/(name(spec)+".npz");meta_path=path.with_suffix(".json")
    if reuse and path.exists() and meta_path.exists():
        meta=json.loads(meta_path.read_text())
        if meta["kernel_sha256"]!=kernel or meta["npz_sha256"]!=sha(path):raise RuntimeError("stale geometry case")
        return meta["summary"]
    p=dataclasses.replace(P0,lam_N=POINTS[point][0]);fp=homogeneous_fixed_point(p)
    W=adjacency((n,n),condition!="open");ns=4 if model=="full" else 3
    initial=random_initial(fp,p,(n,n),seed)[:ns].reshape(ns,-1)
    multiplier=None
    if cv:
        sig=np.sqrt(np.log1p(cv*cv))
        rng=np.random.Generator(np.random.PCG64(seed+1))
        multiplier=np.exp(sig*rng.standard_normal(n*n)-.5*sig*sig)
        multiplier/=multiplier.mean()
    rhs=lambda t,y:lattice_rhs(t,y,p,W,model,multiplier)
    jac=lambda t,y:lattice_jacobian(y,p,W,model)
    def event(t,y):
        d=y[:n*n];return d.std()/d.mean()-.1
    event.direction=1
    end=POINTS[point][1];times=np.arange(0,end+.1,.5)
    start=time.perf_counter()
    sol=solve_ivp(rhs,(0,end),initial.ravel(),method="BDF",jac=jac,rtol=1e-10,atol=1e-13,
                  max_step=.5,t_eval=times,events=event)
    if not sol.success:raise RuntimeError(sol.message)
    D=sol.y[:n*n].T.reshape(-1,n,n);final=sol.y[:,-1].reshape(ns,n,n)
    means=D.mean(axis=(1,2));cvs=D.std(axis=(1,2))/means
    last=final[0];dev=last-last.mean()
    power=abs(np.fft.fft2(dev))**2
    m,c=np.indices((n,n));mu=2*np.cos(2*np.pi*m/n)+2*np.cos(2*np.pi*c/n)+2*np.cos(2*np.pi*(m-c)/n)
    minmu=float(mu.min());shell=abs(mu-minmu)<1e-10
    near=mu<=-2.8
    denom=power.sum()
    exact_min_fraction=float(power[shell].sum()/denom) if denom>0 else np.nan
    near_fraction=float(power[near].sum()/denom) if denom>0 else np.nan
    boundary=np.zeros((n,n),dtype=bool);boundary[:3]=True;boundary[-3:]=True;boundary[:,:3]=True;boundary[:,-3:]=True
    def Cv(x):return float(x.std()/x.mean())
    s=dict(condition=condition,point=point,size=n,seed=seed,model=model,heterogeneity_cv=cv,
        actual_heterogeneity_cv=float(multiplier.std()) if multiplier is not None else 0.,
        onset_01=float(sol.t_events[0][0]) if len(sol.t_events[0]) else np.nan,
        onset_is_homogeneous_instability_readout=(condition=="periodic"),
        final_cv=float(cvs[-1]),bulk_cv=Cv(last[~boundary]),boundary_cv=Cv(last[boundary]),
        final_mean_D=float(last.mean()),final_min_D=float(last.min()),final_max_D=float(last.max()),
        final_scaled_rhs=float(np.max(abs(rhs(end,sol.y[:,-1]).reshape(ns,-1))/fp[:ns,None])),
        min_concentration=float(sol.y.min()),
        graph_min_mu=minmu if condition!="open" else np.nan,
        min_shell_multiplicity=int(shell.sum()) if condition!="open" else 0,
        min_shell_power=exact_min_fraction if condition!="open" else np.nan,
        near_K_band_power=near_fraction if condition!="open" else np.nan,
        graph_correlation=float(dev.ravel()@(W@dev.ravel())/(6*np.dot(dev.ravel(),dev.ravel()))) if denom>0 else np.nan,
        elapsed_seconds=time.perf_counter()-start,nfev=sol.nfev)
    np.savez_compressed(path,times=times,cv=cvs,initial=initial.reshape(ns,n,n),final=final,fp=fp,
                        notch_multiplier=multiplier if multiplier is not None else np.ones(n*n))
    meta_path.write_text(json.dumps(dict(kernel_sha256=kernel,npz_sha256=sha(path),summary=s),indent=2)+"\n",encoding="utf-8")
    print(name(spec),f"onset={s['onset_01']:.3f} CV={s['final_cv']:.3f} rhs={s['final_scaled_rhs']:.2g} ({s['elapsed_seconds']:.1f}s)",flush=True)
    return s


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--reuse",action="store_true")
    parser.add_argument("--workers",type=int,default=4);args=parser.parse_args()
    freeze();(OUT/"geometry").mkdir(exist_ok=True);kernel=digest(__file__)
    specs=list(specifications())
    if len(specs)!=52 or len(set(specs))!=52:raise RuntimeError("geometry design count changed")
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows=list(pool.map(partial(run,kernel=kernel,reuse=args.reuse),specs))
    frame=pd.DataFrame(rows);frame.to_csv(OUT/"geometry_runs.csv",index=False)
    critical=[]
    for n in (23,24,25):
        mus=np.array([circulant_eigenvalue(m,c,n,n) for m in range(n) for c in range(n)])
        mu=float(mus.min())
        def f(ln):
            p=dataclasses.replace(P0,lam_N=ln)
            return G_exact(jacobian_block(mu,homogeneous_fixed_point(p),p))-1
        critical.append(dict(size=n,min_mu=mu,min_mu_multiplicity=int(np.sum(abs(mus-mu)<1e-10)),
            stationary_lam_N=brentq(f,1.2,1.5,xtol=1e-13)))
    pd.DataFrame(critical).to_csv(OUT/"geometry_stationary_thresholds.csv",index=False)
    paired=frame.pivot(index=["condition","point","size","seed","heterogeneity_cv"],columns="model",values="onset_01").reset_index()
    paired["full_DNA_onset_ratio"]=paired["full"]/paired["DNA"]
    paired.to_csv(OUT/"geometry_pairs.csv",index=False)
    summary=dict(kernel_sha256=kernel,runs=len(frame),periodic_main_pairs=18,
        max_scaled_endpoint_drift=float(frame.final_scaled_rhs.max()),
        n_negative_recorded=int((frame.min_concentration<0).sum()),critical=critical,
        note="Open/heterogeneous onset contrast includes deterministic forcing and is not a homogeneous-instability threshold.")
    (OUT/"geometry_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":main()
