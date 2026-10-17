"""Frozen design and content hashes for the seven reviewer-question analyses."""
from __future__ import annotations
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"figures/data/reviewer_extensions"
SEEDS=(2026102001,2026102002,2026102003)


def protocol():
    return {
        "version":1,"date_client":"2026-10-17",
        "Q1_qss":"50 stratified full-unstable reference samples plus P2/P3, alpha logspace(-2,4,61); all graph modes; relative growth error, pole/dS and mode-specific feedback substitution error; empirical thresholds only",
        "Q2_recalibration":"single DNA clock fitted at P2; prediction on lam_N=linspace(1.31,6.5,161), all modes at P2, stable-P0 forced Delta-production response; separate training least-squares clock on lam_N in [1.5,3] and held-out [3.5,6.5]",
        "Q3_geometry":{"periodic_sizes":[23,24,25],"points":{"P2":2.,"P3":5.623},"seeds":list(SEEDS),"amplitude":1e-4,"horizons":{"P0":600.,"P2":600.,"P3":300.},"main_runs":36,
            "controls":"24x24 open P0/P2/P3, both models, first seed (6); periodic P0 both models first seed (2); periodic P2/P3 quenched Notch-production CV=0.05/0.1 both models first seed (8)","total_runs":52,
            "solver":"BDF analytic sparse Jacobian, rtol1e-10 atol1e-13, max_step0.5, recording0.5; identical initial D,N,A and reconstructed S within each model pair",
            "open_interpretation":"raw contact sums, varying degree; homogeneous reference is not an open-graph equilibrium; CV crossing is not interpreted as homogeneous-instability onset",
            "heterogeneity":"mean-one lognormal Notch production, fixed in time, separate PCG64(seed+1) stream; uniform reference is not an equilibrium"},
        "Q4_bifurcations":"all 12 reference DC misses; alpha logspace(-6,2,401), fp and G fixed, every distinct graph mode including uniform; bracket first loss, refine logalpha, report nonzero frequency/transversality as oscillatory linear boundary; no first Lyapunov coefficient or cycle-stability claim; follow unstable-pair discriminant to real-pair collisions",
        "Q5_pattern_stability":"12 distinct endpoint cases: first seed at1e-4 full/DNA P2/P3 (4), largest recorded drift in each point/model main group (4), P2 free-S largest K fraction and P3 free-S largest drift (2), two P2 pure-mode fields (2); sparse Newton refinement then complete dense spectrum; shifted sparse residual check; independent relative1e-4 positive perturbation and unperturbed continuation t200",
        "Q6_growth_prior":"40000 primary evaluations: absolute and bounded symmetric growth errors; cuts in raw growth and growth/maxabs(K-poles); reweight same stored inputs to triangular Hill prior and Gaussian-copula correlation rho0.75 of loglamS/logdS; report effective sample sizes and ten-replicate uncertainty; no new parameter draws or biologically measured prior claim",
        "Q7_mechanism":"matched finite Delta-production lag with Udot=dS*(lamD*hD(Sq(A))-U); exact full-block similarity T=diag(1,1,1,kappaD); verify reference4000 all graph modes and D,N,A forced responses; clamped-A steps factors1.001,1.2,2 over t0..80; two-Hill-threshold concentration-rescaling proof and 100 state checks",
        "outcome_policy":"retain counterexamples and failed checks; no global dynamic theorem, no tissue calibration or automatic promotion of small residuals to stability",
        "new_findings_policy":"decide main-figure promotion after results, document scope and prior standard identities",
    }


def freeze():
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/"protocol.json"
    p=protocol()
    if path.exists() and json.loads(path.read_text())!=p:raise RuntimeError("frozen reviewer protocol changed")
    if not path.exists():path.write_text(json.dumps(p,indent=2)+"\n",encoding="utf-8")
    return p


def digest(script):
    paths=[Path(script).resolve(),Path(__file__).resolve(),ROOT/"scripts/reviewer_models.py"]
    for folder in ("model","geometry","stability","control"):
        paths.extend((ROOT/folder).glob("*.py"))
    text="".join(p.relative_to(ROOT).as_posix()+"\0"+p.read_text(encoding="utf-8")+"\0" for p in sorted(set(paths)))
    return hashlib.sha256(text.encode()).hexdigest()


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
