"""Complete patterned-state spectra and paired perturbation recovery checks."""
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
from scipy.linalg import eigvals
from scipy.sparse.linalg import spsolve,eigs
from model.params import P0
from model.fixed_point import homogeneous_fixed_point
from model.equations import hill_ASC
from scripts.initial_condition_robustness import case_name,kernel_hash
from scripts.reviewer_models import adjacency,lattice_rhs,lattice_jacobian
from scripts.reviewer_protocol import OUT,freeze,digest,sha


def selections():
    frame=pd.read_csv(ROOT/"figures/data/initial_conditions/runs.csv")
    rows=[];used=set()
    def add(row,reason):
        spec=(row.point,float(row.amplitude),int(row.seed),row.model,row.downstream)
        name=case_name(spec)
        if name in used:return False
        used.add(name);rows.append(dict(name=name,point=row.point,model=row.model,reason=reason));return True
    main=frame[frame.downstream=="qss"]
    for (point,model),g in main.groupby(["point","model"]):
        add(g[(g.seed==2026101901)&np.isclose(g.amplitude,1e-4)].iloc[0],"first seed, amplitude 1e-4")
    for (point,model),g in main.groupby(["point","model"]):
        for _,r in g.sort_values("final_scaled_rhs",ascending=False).iterrows():
            if add(r,"largest recorded drift among distinct remaining endpoints"):break
    for point,key in (("P2","final_k_fraction"),("P3","final_scaled_rhs")):
        g=frame[(frame.point==point)&(frame.downstream=="free")]
        add(g.sort_values(key,ascending=False).iloc[0],"free-S control: largest "+key)
    for model in ("full","DNA"):
        rows.append(dict(name="pure_P2_"+model,point="P2",model=model,reason="reconstructed pure-mode endpoint"))
    if len(rows)!=12 or len({r["name"] for r in rows})!=12:raise RuntimeError("endpoint selection count")
    return rows


def refine(y,p,W,model,fp):
    ns=len(fp);scale=np.repeat(fp,W.shape[0]);start=y.copy();history=[]
    for _ in range(9):
        f=lattice_rhs(0,y,p,W,model);res=float(np.max(abs(f)/scale));history.append(res)
        if res<=1e-11:break
        step=spsolve(lattice_jacobian(y,p,W,model),-f)
        for power in range(30):
            candidate=y+step*2.**(-power)
            if candidate.min()>0 and np.max(abs(lattice_rhs(0,candidate,p,W,model))/scale)<res:
                y=candidate;break
        else:raise RuntimeError("Newton line search failed")
    if history[-1]>1e-11:raise RuntimeError("Newton did not reach frozen tolerance")
    return y,history,float(np.linalg.norm((y-start)/scale)/np.linalg.norm(start/scale))


def run(index,row,kernel,reuse=False):
    path=OUT/"pattern_stability"/(row["name"]+".npz");meta_path=path.with_suffix(".json")
    pure=row["name"].startswith("pure_")
    source=ROOT/("figures/data/pattern_final_fields.npz" if pure else "figures/data/initial_conditions/cases/"+row["name"]+".npz")
    source_sha=sha(source)
    if not pure:
        old=json.loads(source.with_suffix(".json").read_text())
        if old["kernel_sha256"]!=kernel_hash() or old["npz_sha256"]!=source_sha:raise RuntimeError("stale original endpoint")
    if reuse and path.exists() and meta_path.exists():
        meta=json.loads(meta_path.read_text())
        if meta["kernel_sha256"]!=kernel or meta["source_sha256"]!=source_sha or meta["npz_sha256"]!=sha(path):raise RuntimeError("stale patterned case")
        return meta["summary"]
    start=time.perf_counter();model=row["model"];ns=4 if model=="full" else 3
    p=dataclasses.replace(P0,lam_N=2. if row["point"]=="P2" else 5.623)
    fp=homogeneous_fixed_point(p);W=adjacency((24,24))
    z=np.load(source)
    if pure:
        D=z["D_"+model+"_P2"].ravel()
        N=p.lam_N/(p.d_N+p.f_N*(W@D)+p.k_cis*D)
        A=p.f_N*N*(W@D)/p.d_A;S=p.lam_ASC/p.d_ASC*hill_ASC(A,p)
        original=np.array([D,N,A,S])[:ns].ravel()
    else:original=z["final"].ravel()
    y,history,change=refine(original,p,W,model,fp[:ns])
    J=lattice_jacobian(y,p,W,model)
    poles=eigvals(J.toarray(),overwrite_a=True,check_finite=False)
    leading=poles[np.argmax(poles.real)]
    # Shift off the dense leading root, independently recover a nearby eigenpair.
    v,w=eigs(J,k=3,sigma=complex(leading)+1e-6,which="LM",tol=1e-11,v0=np.ones(len(y)))
    best=int(np.argmin(abs(v-leading)));vec=w[:,best]
    sparse_error=float(abs(v[best]-leading))
    residual=float(np.linalg.norm(J@vec-v[best]*vec)/(np.linalg.norm(J.data)*np.linalg.norm(vec)))
    rng=np.random.Generator(np.random.PCG64(2026102201+index))
    xi=rng.uniform(-1,1,(ns,576));xi-=xi.mean(axis=1,keepdims=True);xi/=abs(xi).max(axis=1,keepdims=True)
    perturb=y*(1+1e-4*xi.ravel());times=np.arange(0,200.01,.5)
    rhs=lambda t,state:lattice_rhs(t,state,p,W,model)
    jac=lambda t,state:lattice_jacobian(state,p,W,model)
    sols=[]
    for initial in (y,perturb):
        sol=solve_ivp(rhs,(0,200),initial,method="BDF",jac=jac,rtol=1e-10,atol=1e-13,max_step=.5,t_eval=times)
        if not sol.success or sol.y.min()<=0:raise RuntimeError("invalid patterned continuation")
        sols.append(sol.y)
    scale=np.repeat(fp[:ns],576)
    distance=np.sqrt(np.mean(((sols[1]-sols[0])/scale[:,None])**2,axis=0))
    s=dict(**row,seed=2026102201+index,dimension=len(y),newton_residual=history[-1],newton_iterations=len(history)-1,
           refinement_relative_change=change,leading_real=float(leading.real),leading_imag=float(leading.imag),
           n_RHP=int(np.sum(poles.real>1e-9)),n_near_imaginary_axis=int(np.sum(abs(poles.real)<=1e-9)),
           sparse_dense_leading_error=sparse_error,sparse_eigenpair_residual=residual,
           recovery_final_over_initial=float(distance[-1]/distance[0]),
           recovery_max_over_initial=float(distance.max()/distance[0]),
           unperturbed_scaled_drift=float(np.max(abs(sols[0][:,-1]-y)/scale)),elapsed_seconds=time.perf_counter()-start)
    np.savez_compressed(path,original=original.reshape(ns,24,24),refined=y.reshape(ns,24,24),fp=fp,
                        poles=poles,times=times,paired_distance=distance,unperturbed_final=sols[0][:,-1],
                        perturbed_final=sols[1][:,-1],newton_history=history)
    meta_path.write_text(json.dumps(dict(kernel_sha256=kernel,source_sha256=source_sha,npz_sha256=sha(path),summary=s),indent=2)+"\n",encoding="utf-8")
    print(row["name"],f"leading={leading.real:.6g} RHP={s['n_RHP']} recovery={s['recovery_final_over_initial']:.3g} ({s['elapsed_seconds']:.1f}s)",flush=True)
    return s


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--reuse",action="store_true");parser.add_argument("--workers",type=int,default=4);args=parser.parse_args()
    freeze();(OUT/"pattern_stability").mkdir(exist_ok=True);kernel=digest(__file__)
    rows=selections();(OUT/"pattern_selection.json").write_text(json.dumps(rows,indent=2)+"\n")
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results=list(pool.map(partial(run,kernel=kernel,reuse=args.reuse),range(len(rows)),rows))
    frame=pd.DataFrame(results);frame.to_csv(OUT/"pattern_stability.csv",index=False)
    summary=dict(kernel_sha256=kernel,cases=len(frame),all_strictly_stable=bool((frame.leading_real< -1e-9).all()),
        leading_real_min=float(frame.leading_real.min()),leading_real_max=float(frame.leading_real.max()),
        total_RHP=int(frame.n_RHP.sum()),max_newton_change=float(frame.refinement_relative_change.max()),
        max_refined_residual=float(frame.newton_residual.max()),max_sparse_dense_error=float(frame.sparse_dense_leading_error.max()),
        max_eigenpair_residual=float(frame.sparse_eigenpair_residual.max()),recovery_ratio_max=float(frame.recovery_final_over_initial.max()),
        scope="12 selected 24x24 periodic states; local linear stability and one finite recovery check each, no global attractor theorem")
    (OUT/"pattern_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8");print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":main()
