"""Time-step, tolerance, solver and output-cadence checks without changing data.

The time-continuous Jacobian spectra are not obtained by time stepping. This
audit separately tests their reproduction by time-discrete dynamics and the
nonlinear lattice trajectories. Freeze the design with --protocol-only.
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
from scipy.linalg import expm

from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.params import P0
from model.reduced import reduce_S
from scripts.initial_condition_robustness import (OUT as IC_OUT, POINTS, SEEDS,
    case_name, initial_fields, jacobian_sparsity, kernel_hash as ic_kernel_hash,
    spatial_metrics)
from simulation.nonlinear_ode import pack, rhs_lattice, rhs_lattice_DNA
from simulation.perturbation import eigenmode_state, k_mode_amplitude, leading_real_eigvec
from stability.mode_scan import UNIQUE_MU

OUT = ROOT / "figures/data/timestep_audit"
FIXED_DT = (.08, .04, .02, .01)
CONFIGS = (
    ("BDF_h0.5", "BDF", .5, 1e-9, 1e-11),
    ("BDF_h0.25", "BDF", .25, 1e-9, 1e-11),
    ("BDF_h0.125", "BDF", .125, 1e-9, 1e-11),
    ("BDF_tight", "BDF", .125, 1e-11, 1e-14),
    ("Radau_tight", "Radau", .25, 1e-11, 1e-14),
)
GROWTH_END = {("P2","full"):70.4, ("P2","DNA"):4.8,
              ("P3","full"):32.8, ("P3","DNA"):2.72}


def lattice_cases():
    for point in POINTS:
        for amp in (1e-6, 1e-4, 1e-2):
            for model in ("full","DNA"):
                yield point,amp,SEEDS[0],model,"qss"
        yield point,1e-4,SEEDS[0],"full","free"


def protocol():
    return {
        "version":1,
        "scope":"time discretisation and arithmetic, not a spatial continuum limit or biological validation",
        "reference_inputs":"existing first-seed lattice fields; all three amplitudes and both points/models; two independent-S controls",
        "lattice_cases":[list(s) for s in lattice_cases()], "n_lattice_cases":14,
        "adaptive_configs":[dict(name=n,method=m,max_step=h,rtol=r,atol=a) for n,m,h,r,a in CONFIGS],
        "adaptive_total_runs":70,
        "adaptive_horizons":{k:v[1] for k,v in POINTS.items()},
        "baseline":"previous BDF runs with no max_step bound, rtol=1e-9, atol=1e-11",
        "adaptive_readout":"first upward CV(D)=0.1 and 0.5 events; final raw cellwise fields, CV, K fraction and scaled drift; actual accepted internal step sizes",
        "fixed_rk4_dt":list(FIXED_DT), "fixed_total_lattice_runs":16,
        "fixed_initialisation":"own exact K eigenvector; D amplitude 1e-6 D*, no concentration clipping",
        "fixed_compare_times":"every 0.08 time units shared by all fixed dt runs and the high-accuracy Radau reference",
        "fixed_growth_window":"same common times no later than log(1000)/leading_eigenvalue, independent of each computed trajectory",
        "fixed_reference":"Radau rtol=1e-12, atol=1e-15, max_step=0.01 through the early growth horizon; not a long-time endpoint test",
        "fixed_growth_horizons":{f"{p}_{m}":t for (p,m),t in GROWTH_END.items()},
        "matrix_propagator_test":"P0/P2/P3, full/DNA, every distinct graph eigenvalue: RK4 stability polynomial eigenvalues versus continuous spectrum; matrix-power propagator at t=2 versus expm(2J)",
        "output_cadence":[.5,.25,.125],
        "cadence_test":"resample the identical dense Radau solution; first-recorded-amplitude growth-window rule from the reference manuscript; internal steps unchanged",
        "assessment_thresholds":{
            "adaptive_event_relative_change":1e-3,
            "adaptive_final_relative_L2_all_species":1e-4,
            "adaptive_final_k_abs_change":1e-5,
            "adaptive_stationary_scaled_drift":1e-7,
            "fixed_growth_relative_error_vs_eigenvalue":5e-4,
            "fixed_relative_D_perturbation_error_vs_Radau":1e-4,
            "cadence_relative_slope_change":1e-4,
            "minimum_accepted_concentration":-1e-10},
        "assessment":"thresholds fixed before new runs; violations are retained and reported, never removed or relabelled as nonexistence",
        "limits":"14 representative initial states, not every seed or parameter sample; endpoint agreement does not prove attractor stability; finite-step convergence is not a mathematical exclusion of all numerical error",
    }


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def kernel_hash():
    paths=[Path(__file__).resolve(),ROOT/"scripts/initial_condition_robustness.py"]
    for folder in ("model","geometry","simulation","stability"):
        paths.extend((ROOT/folder).glob("*.py"))
    text="".join(p.relative_to(ROOT).as_posix()+"\0"+p.read_text(encoding="utf-8")+"\0" for p in sorted(paths))
    return hashlib.sha256(text.encode()).hexdigest()


def rk4_integrate(rhs, y0, dt, n_steps, record_every=1):
    """Classical fixed-step RK4, no adaptation or clipping, record a nested grid."""
    if dt<=0 or n_steps<1 or record_every<1 or n_steps%record_every:
        raise ValueError("invalid fixed time grid")
    y=np.array(y0,dtype=float,copy=True)
    states=[y.copy()]
    min_stage=float(y.min())
    for k in range(n_steps):
        t=k*dt
        a=rhs(t,y)
        yb=y+dt*a/2
        b=rhs(t+dt/2,yb)
        yc=y+dt*b/2
        c=rhs(t+dt/2,yc)
        yd=y+dt*c
        d=rhs(t+dt,yd)
        y=y+dt*(a+2*b+2*c+d)/6
        min_stage=min(min_stage,float(yb.min()),float(yc.min()),float(yd.min()),float(y.min()))
        if not np.all(np.isfinite(y)):
            raise RuntimeError("nonfinite RK4 trajectory")
        if (k+1)%record_every==0:
            states.append(y.copy())
    return np.arange(len(states))*dt*record_every,np.array(states),min_stage


def baseline(spec):
    path=IC_OUT/"cases"/(case_name(spec)+".npz")
    meta=json.loads(path.with_suffix(".json").read_text())
    if meta["kernel_sha256"] != ic_kernel_hash() or meta["npz_sha256"]!=sha(path):
        raise RuntimeError(f"stale baseline: {path}")
    with np.load(path) as z:
        arrays={k:z[k].copy() for k in ("initial","final","fp")}
    return arrays,meta["summary"],meta["npz_sha256"]


def adaptive_run(task,digest,reuse=False):
    spec,config=task
    point,amp,seed,model,downstream=spec
    label,method,h,rtol,atol=config
    stem=case_name(spec)+"_"+label
    path=OUT/"cases"/(stem+".npz")
    meta_path=path.with_suffix(".json")
    arrays,base,basehash=baseline(spec)
    if reuse and path.exists() and meta_path.exists():
        meta=json.loads(meta_path.read_text())
        if meta["kernel_sha256"]!=digest or meta["baseline_sha256"]!=basehash or meta["npz_sha256"]!=sha(path):
            raise RuntimeError(f"stale audit case: {stem}")
        return meta["summary"]
    p=dataclasses.replace(P0,lam_N=POINTS[point][0])
    rhs=rhs_lattice if model=="full" else rhs_lattice_DNA
    ns=len(arrays["initial"])
    end=POINTS[point][1]

    def event01(t,y):
        d=y[:576]
        return d.std()/d.mean()-.1

    def event05(t,y):
        d=y[:576]
        return d.std()/d.mean()-.5

    event01.direction=event05.direction=1
    start=time.perf_counter()
    sol=solve_ivp(lambda t,y:rhs(t,y,p),(0,end),pack(*arrays["initial"]),method=method,
        jac_sparsity=jacobian_sparsity(model),rtol=rtol,atol=atol,max_step=h,
        events=(event01,event05))
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(f"{stem}: {sol.message}")
    fields=sol.y[:,-1].reshape(ns,24,24)
    delta=fields-arrays["final"]
    dcv,kfrac=spatial_metrics(fields[0])
    events=[float(e[0]) if len(e) else None for e in sol.t_events]
    changes=[abs(events[i]-base[key])/base[key] if events[i] is not None else None
             for i,key in enumerate(("onset_01","onset_05"))]
    dt=np.diff(sol.t)
    scaled_rhs=float(np.max(np.abs(rhs(end,sol.y[:,-1],p).reshape(ns,24,24))/arrays["fp"][:ns,None,None]))
    row=dict(point=point,amplitude=amp,seed=seed,model=model,downstream=downstream,
        config=label,method=method,max_step=h,rtol=rtol,atol=atol,
        onset_01=events[0],onset_05=events[1],onset_01_relative_change=changes[0],onset_05_relative_change=changes[1],
        final_relative_L2_all=float(np.linalg.norm(delta)/np.linalg.norm(arrays["final"])),
        final_relative_L2_D=float(np.linalg.norm(delta[0])/np.linalg.norm(arrays["final"][0])),
        final_D_maxabs_over_Dstar=float(np.max(np.abs(delta[0]))/arrays["fp"][0]),
        final_k_fraction=float(kfrac),k_fraction_abs_change=abs(float(kfrac)-base["final_k_fraction"]),
        final_cv=float(dcv),cv_relative_change=abs(float(dcv)/base["final_cv"]-1),
        final_scaled_rhs=scaled_rhs,min_accepted_concentration=float(sol.y.min()),
        accepted_steps=len(dt),actual_step_min=float(dt.min()),actual_step_median=float(np.median(dt)),
        actual_step_max=float(dt.max()),fraction_steps_at_cap=float(np.mean(dt>=.999*h)),
        nfev=sol.nfev,elapsed_seconds=time.perf_counter()-start)
    np.savez_compressed(path,internal_times=sol.t,final=fields,initial=arrays["initial"])
    meta_path.write_text(json.dumps(dict(kernel_sha256=digest,baseline_sha256=basehash,npz_sha256=sha(path),summary=row),indent=2)+"\n",encoding="utf-8")
    print(f"{stem}: onset diff={max(c or 0 for c in changes):.2g}; final L2={row['final_relative_L2_all']:.2g}; steps={len(dt)} ({row['elapsed_seconds']:.1f}s)",flush=True)
    return row


def spectral_checks():
    rows=[]
    for point,lam_N in (("P0",1.),("P2",2.),("P3",5.623)):
        p=dataclasses.replace(P0,lam_N=lam_N)
        fp=homogeneous_fixed_point(p)
        for model in ("full","DNA"):
            for mu in UNIQUE_MU:
                J=jacobian_block(mu,fp,p)
                if model=="DNA":J=reduce_S(J)
                values=np.linalg.eigvals(J)
                target=expm(2*J)
                I=np.eye(len(J))
                for h in FIXED_DT:
                    M=h*J
                    R=I+M+M@M/2+np.linalg.matrix_power(M,3)/6+np.linalg.matrix_power(M,4)/24
                    discrete=np.log(np.abs(np.linalg.eigvals(R)))/h
                    numerical=np.linalg.matrix_power(R,int(round(2/h)))
                    rows.append(dict(point=point,model=model,mu=float(mu),dt=h,
                        continuous_max_real=float(values.real.max()),discrete_max_growth=float(discrete.max()),
                        continuous_rhp=int((values.real>0).sum()),discrete_unstable=int((discrete>0).sum()),
                        propagator_relative_error=float(np.linalg.norm(numerical-target)/np.linalg.norm(target))))
    frame=pd.DataFrame(rows)
    frame.to_csv(OUT/"spectral_checks.csv",index=False)
    return frame


def growth_run(spec):
    point,model=spec
    p=dataclasses.replace(P0,lam_N=POINTS[point][0])
    fp=homogeneous_fixed_point(p)
    ns=4 if model=="full" else 3
    J=jacobian_block(-3,fp,p)
    if model=="DNA":J=reduce_S(J)
    rate,vec=leading_real_eigvec(J)
    y0=pack(*eigenmode_state(fp[:ns],vec,1e-6))
    rhs0=rhs_lattice if model=="full" else rhs_lattice_DNA
    rhs=lambda t,y:rhs0(t,y,p)
    end=GROWTH_END[point,model]
    common=np.arange(int(round(end/.08))+1)*.08
    fit=common <= np.log(1000)/rate
    ref=solve_ivp(rhs,(0,end),y0,method="Radau",rtol=1e-12,atol=1e-15,
        max_step=.01,jac_sparsity=jacobian_sparsity(model),dense_output=True)
    if not ref.success:raise RuntimeError(ref.message)
    reference=ref.sol(common).T
    def amplitudes(states):
        return np.array([k_mode_amplitude(y[:576].reshape(24,24)) for y in states])
    reference_amp=amplitudes(reference)
    reference_slope=float(np.polyfit(common[fit],np.log(reference_amp[fit]),1)[0])
    rows=[]
    for h in FIXED_DT:
        every=int(round(.08/h))
        start=time.perf_counter()
        times,states,minimum=rk4_integrate(rhs,y0,h,int(round(end/h)),every)
        amps=amplitudes(states)
        slope=float(np.polyfit(times[fit],np.log(amps[fit]),1)[0])
        dn=(states[:,:576]-reference[:,:576])
        norm=np.linalg.norm(reference[:,:576]-fp[0],axis=1)
        error=float(np.max(np.linalg.norm(dn,axis=1)/norm))
        rows.append(dict(point=point,model=model,dt=h,predicted=rate,fitted_slope=slope,
            relative_slope_error=abs(slope/rate-1),reference_slope=reference_slope,
            slope_change_vs_reference=abs(slope/reference_slope-1),
            max_relative_D_perturbation_error=error,min_stage_concentration=minimum,
            fit_end=float(times[fit][-1]),elapsed_seconds=time.perf_counter()-start))
        np.savez_compressed(OUT/"cases"/f"RK4_{point}_{model}_h{h}.npz",times=times,amplitude=amps,
            reference_amplitude=reference_amp,final=states[-1].reshape(ns,24,24))
        print(f"RK4 {point} {model} h={h}: slope error={rows[-1]['relative_slope_error']:.2g}; field error={error:.2g}",flush=True)
    cadence=[]
    for record in (.5,.25,.125):
        times=np.arange(0,end+1e-10,record)
        amp=amplitudes(ref.sol(times).T)
        reached=np.flatnonzero(amp/fp[0]>=1e-3)
        if not len(reached):raise RuntimeError("cadence window did not reach its threshold")
        n=int(reached[0])+1
        slope=float(np.polyfit(times[:n],np.log(amp[:n]),1)[0])
        cadence.append(dict(point=point,model=model,record_interval=record,fitted_slope=slope,
            window_end=float(times[n-1]),predicted=rate,relative_slope_error=abs(slope/rate-1),
            reference_internal_steps=len(ref.t)-1))
    return rows,cadence


def assess(adaptive,growth,cadence,spectral):
    th=protocol()["assessment_thresholds"]
    events=adaptive[["onset_01_relative_change","onset_05_relative_change"]]
    cadence_changes=[]
    for _,g in cadence.groupby(["point","model"]):
        ref=float(g.loc[g.record_interval==.125,"fitted_slope"].iloc[0])
        cadence_changes.extend(abs(g.fitted_slope/ref-1).tolist())
    flags={
        "adaptive_events":bool(not events.isna().any().any() and events.to_numpy().max()<th["adaptive_event_relative_change"]),
        "adaptive_fields":bool(adaptive.final_relative_L2_all.max()<th["adaptive_final_relative_L2_all_species"]),
        "adaptive_k_fraction":bool(adaptive.k_fraction_abs_change.max()<th["adaptive_final_k_abs_change"]),
        "adaptive_stationarity":bool(adaptive.final_scaled_rhs.max()<th["adaptive_stationary_scaled_drift"]),
        "fixed_growth":bool(growth.relative_slope_error.max()<th["fixed_growth_relative_error_vs_eigenvalue"]),
        "fixed_fields":bool(growth.max_relative_D_perturbation_error.max()<th["fixed_relative_D_perturbation_error_vs_Radau"]),
        "output_cadence":bool(max(cadence_changes)<th["cadence_relative_slope_change"]),
        "positive_to_threshold":bool(min(adaptive.min_accepted_concentration.min(),growth.min_stage_concentration.min())>th["minimum_accepted_concentration"]),
        "spectral_classifications":bool((spectral.continuous_rhp==spectral.discrete_unstable).all()),
    }
    return dict(assessments=flags,all_assessments_pass=all(flags.values()),
        adaptive_runs=len(adaptive),fixed_lattice_runs=len(growth),spectral_comparisons=len(spectral),
        max_adaptive_event_relative_change=float(events.to_numpy().max()),
        max_final_relative_L2_all=float(adaptive.final_relative_L2_all.max()),
        max_final_relative_L2_D=float(adaptive.final_relative_L2_D.max()),
        max_k_fraction_abs_change=float(adaptive.k_fraction_abs_change.max()),
        max_final_scaled_rhs=float(adaptive.final_scaled_rhs.max()),
        min_accepted_or_RK_stage_concentration=float(min(adaptive.min_accepted_concentration.min(),growth.min_stage_concentration.min())),
        max_RK4_relative_slope_error=float(growth.relative_slope_error.max()),
        max_RK4_relative_D_perturbation_error=float(growth.max_relative_D_perturbation_error.max()),
        max_cadence_relative_slope_change=float(max(cadence_changes)),
        n_discrete_spectral_disagreements=int((spectral.continuous_rhp!=spectral.discrete_unstable).sum()),
        spectral_max_propagator_error_by_dt={str(h):float(g.propagator_relative_error.max()) for h,g in spectral.groupby("dt")},
        conclusions_scope=protocol()["limits"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol-only",action="store_true")
    parser.add_argument("--reuse",action="store_true")
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    if args.workers<1:parser.error("--workers must be positive")
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/"protocol.json"
    if path.exists() and json.loads(path.read_text())!=protocol():
        raise RuntimeError("audit protocol changed; preserve the original and use a new study")
    if not path.exists():path.write_text(json.dumps(protocol(),indent=2)+"\n",encoding="utf-8")
    if args.protocol_only:
        print(path)
        return
    (OUT/"cases").mkdir(exist_ok=True)
    digest=kernel_hash()
    spectral=spectral_checks()
    tasks=[(s,c) for s in lattice_cases() for c in CONFIGS]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        adaptive=pd.DataFrame(pool.map(partial(adaptive_run,digest=digest,reuse=args.reuse),tasks))
    adaptive.to_csv(OUT/"adaptive_runs.csv",index=False)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        runs=list(pool.map(growth_run,GROWTH_END))
    growth=pd.DataFrame([r for rows,_ in runs for r in rows])
    cadence=pd.DataFrame([r for _,rows in runs for r in rows])
    growth.to_csv(OUT/"fixed_growth.csv",index=False)
    cadence.to_csv(OUT/"output_cadence.csv",index=False)
    summary=assess(adaptive,growth,cadence,spectral)
    summary["kernel_sha256"]=digest
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    main()
