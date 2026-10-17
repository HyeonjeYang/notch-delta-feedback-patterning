"""QSS validity, clock recalibration, bifurcation paths, priors and lag identity."""
from __future__ import annotations
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import dataclasses,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import ndtri,expit
from scipy.stats import t
from control.state_space import G_exact
from model.params import P0,NotchDeltaParams
from model.equations import hill_D,hill_ASC,drift
from model.fixed_point import homogeneous_fixed_point
from model.jacobian import jacobian_block
from model.reduced import reduce_S
from stability.mode_scan import UNIQUE_MU
from stability.sensitivity import PARAM_NAMES
from scripts.robustness_study import mode_blocks
from scripts.validate_robustness_cases import load_samples
from scripts.reviewer_models import production_lag_block,lag_similarity
from scripts.reviewer_protocol import OUT,freeze,digest,sha


def params(row):return NotchDeltaParams(**{k:float(row[k]) for k in PARAM_NAMES})


def qss_analysis(reference,inputs):
    unstable=reference[reference.lambda_pattern>0].sort_values(["lambda_pattern","sample_id"])
    selected=unstable.iloc[np.rint(np.linspace(0,len(unstable)-1,50)).astype(int)].sample_id.astype(int).tolist()
    families=[(f"sample_{sid}",params(inputs.loc[sid])) for sid in selected]
    families.extend((name,dataclasses.replace(P0,lam_N=ln)) for name,ln in (("P2",2.),("P3",5.623)))
    pd.DataFrame([dict(family=name,**dataclasses.asdict(p)) for name,p in families]).to_csv(OUT/"qss_families.csv",index=False)
    rows=[]
    for name,p in families:
        fp=homogeneous_fixed_point(p);blocks=mode_blocks(fp,p)
        reduced=np.stack([reduce_S(J) for J in blocks])
        red=float(np.linalg.eigvals(reduced[:-1]).real.max())
        for alpha in np.logspace(-2,4,61):
            B=blocks.copy();B[:,3,:]*=alpha
            ev=np.linalg.eigvals(B);r=ev.real.max(axis=1)
            mu_index=int(np.argmax(r[:-1]));J=B[mu_index]
            val,vec=np.linalg.eig(J);k=int(np.argmax(val.real));s=val[k]
            inv=np.linalg.inv(vec)
            cond=float(np.linalg.norm(vec[:,k])*np.linalg.norm(inv[k]))
            full=float(s.real);rate=alpha*p.d_ASC
            rows.append(dict(family=name,alpha=alpha,full_growth=full,DNA_growth=red,
                relative_growth_error=abs(red/full-1) if full>1e-8 else np.nan,
                scaled_absolute_error=abs(red-full)/p.d_ASC,
                pole_over_S_rate=abs(s)/rate,feedback_substitution_error=abs(s/(s+rate)),
                eigenvalue_condition=cond,best_mu=float(UNIQUE_MU[mu_index]),
                uniform_growth=float(r[-1]),gain_K=G_exact(blocks[0])))
    frame=pd.DataFrame(rows);frame.to_csv(OUT/"qss_validity.csv",index=False)
    accuracy=[]
    for name,g in frame.groupby("family",sort=False):
        for target in (.05,.1):
            # Require all later evaluated alphas to pass, not an isolated crossing.
            good=(g.relative_growth_error<target).to_numpy()
            tail=np.logical_and.accumulate(good[::-1])[::-1]
            ids=np.flatnonzero(tail)
            accuracy.append(dict(family=name,target=target,alpha_required=float(g.iloc[ids[0]].alpha) if len(ids) else np.nan))
    pd.DataFrame(accuracy).to_csv(OUT/"qss_accuracy_thresholds.csv",index=False)
    return dict(families=len(families),evaluations=len(frame),
        full_faster_than_DNA=int(((frame.full_growth>1e-8)&(frame.full_growth>frame.DNA_growth+1e-9)).sum()),
        other_nonuniform_dominant=int((frame.best_mu>-3+1e-9).sum()))


def recalibration():
    def blocks(ln):
        p=dataclasses.replace(P0,lam_N=ln);fp=homogeneous_fixed_point(p)
        J=jacobian_block(-3,fp,p)
        return J,reduce_S(J)
    J,R=blocks(2.)
    clock=float(np.linalg.eigvals(J).real.max()/np.linalg.eigvals(R).real.max())
    rows=[]
    for ln in np.unique(np.r_[np.linspace(1.31,6.5,161),2.,5.623]):
        J,R=blocks(ln)
        full=float(np.linalg.eigvals(J).real.max());red=float(np.linalg.eigvals(R).real.max())
        rows.append(dict(lam_N=ln,full_growth=full,DNA_growth=red,clock_fit_P2=clock,
            calibrated_growth=clock*red,calibrated_relative_error=(clock*red/full-1)))
    f=pd.DataFrame(rows)
    train=f[(f.lam_N>=1.5)&(f.lam_N<=3.)]
    ls=float(np.dot(train.full_growth,train.DNA_growth)/np.dot(train.DNA_growth,train.DNA_growth))
    f["training_LS_clock"]=ls;f["LS_prediction"]=ls*f.DNA_growth
    f["LS_relative_error"]=f.LS_prediction/f.full_growth-1
    f.to_csv(OUT/"clock_recalibration.csv",index=False)
    p=dataclasses.replace(P0,lam_N=2);fp=homogeneous_fixed_point(p)
    mode_rows=[]
    for mu in UNIQUE_MU:
        J=jacobian_block(mu,fp,p);R=reduce_S(J)
        a=float(np.linalg.eigvals(J).real.max());b=float(np.linalg.eigvals(R).real.max())
        mode_rows.append(dict(mu=mu,full_growth=a,DNA_growth=b,calibrated_growth=clock*b))
    pd.DataFrame(mode_rows).to_csv(OUT/"clock_mode_predictions.csv",index=False)
    fp=homogeneous_fixed_point(P0);J=jacobian_block(-3,fp,P0);R=reduce_S(J)
    B=np.array([1.,0.,0.,0.]);b=B[:3]
    fr=[]
    for w in np.r_[0.,np.logspace(-4,3,401)]:
        H=np.linalg.solve(1j*w*np.eye(4)-J,B)[0]
        Hred=np.linalg.solve(1j*w*np.eye(3)-clock*R,clock*b)[0]
        fr.append(dict(omega=w,full_re=H.real,full_im=H.imag,calibrated_re=Hred.real,calibrated_im=Hred.imag))
    pd.DataFrame(fr).to_csv(OUT/"clock_frequency_response.csv",index=False)
    hold=f[(f.lam_N>=3.5)&(f.lam_N<=6.5)]
    return dict(clock_at_P2=clock,P3_error=float(f.loc[np.isclose(f.lam_N,5.623),"calibrated_relative_error"].iloc[0]),
        training_clock=ls,heldout_LS_error_min=float(hold.LS_relative_error.min()),heldout_LS_error_max=float(hold.LS_relative_error.max()))


def bifurcations(reference,inputs):
    ids=reference.loc[(reference.lambda_pattern>0)&(reference.G_minus3<1),"sample_id"].astype(int)
    curves=[];boundaries=[];collisions=[]
    for sid in ids:
        p=params(inputs.loc[sid]);fp=homogeneous_fixed_point(p);base=mode_blocks(fp,p)
        def eigs(loga):
            B=base.copy();B[:,3,:]*=np.exp(loga)
            return np.linalg.eigvals(B)
        def margin(loga):return float(eigs(loga).real.max())
        grid=np.linspace(np.log(1e-6),np.log(1e2),401)
        margins=[]
        for la in grid:
            ev=eigs(la);margins.append(ev.real.max())
            for e in ev[0]:curves.append(dict(sample_id=sid,alpha=np.exp(la),re=e.real,im=e.imag,gain_K=G_exact(base[0])))
        hit=np.flatnonzero((np.array(margins[:-1])<0)&(np.array(margins[1:])>=0))
        if not len(hit):raise RuntimeError(f"no first stable-to-unstable crossing in frozen path: {sid}")
        i=int(hit[0]);la=brentq(margin,grid[i],grid[i+1],xtol=1e-13)
        ev=eigs(la);mi,ei=np.unravel_index(np.argmax(ev.real),ev.shape);pole=ev[mi,ei]
        slope=(margin(la+1e-5)-margin(la-1e-5))/2e-5
        boundaries.append(dict(sample_id=sid,alpha=np.exp(la),mu=float(UNIQUE_MU[mi]),
            frequency=abs(pole.imag),real_residual=abs(pole.real),transverse_logalpha_slope=slope,
            G_K=G_exact(base[0]),n_rhp_before=int((eigs(la-1e-4).real>0).sum()),
            n_rhp_after=int((eigs(la+1e-4).real>0).sum())))
        # The largest-real-part pair remains the unstable pair after this boundary.
        def pair_discriminant(z):
            v=eigs(z)[0];pair=v[np.argsort(v.real)[-2:]]
            return float(((pair[0]-pair[1])**2).real)
        g=np.linspace(la+1e-4,0,250);dis=np.array([pair_discriminant(z) for z in g])
        h=np.flatnonzero((dis[:-1]<0)&(dis[1:]>0))
        if len(h):
            k=int(h[0]);root=brentq(pair_discriminant,g[k],g[k+1],xtol=1e-12)
            vals=eigs(root)[0];pair=vals[np.argsort(vals.real)[-2:]]
            collisions.append(dict(sample_id=sid,alpha=np.exp(root),pair_real_mean=float(pair.real.mean()),
                pair_split=float(abs(pair[0]-pair[1])),G_K=G_exact(base[0])))
    pd.DataFrame(curves).to_csv(OUT/"bifurcation_poles.csv.gz",index=False,compression="gzip")
    pd.DataFrame(boundaries).to_csv(OUT/"oscillatory_boundaries.csv",index=False)
    pd.DataFrame(collisions).to_csv(OUT/"real_pair_collisions.csv",index=False)
    return dict(paths=len(boundaries),nonzero_frequency_crossings=int(sum(x["frequency"]>1e-8 for x in boundaries)),
        all_first_at_K=all(abs(x["mu"]+3)<1e-9 for x in boundaries),real_pair_collisions=len(collisions),
        max_crossing_real_residual=max(x["real_residual"] for x in boundaries))


def weighted_median(values,weights):
    ix=np.argsort(values);v=np.asarray(values)[ix];w=np.asarray(weights)[ix]
    return float(v[np.searchsorted(np.cumsum(w),.5*w.sum())])


def growth_and_priors():
    data=load_samples();d=data[data.scenario=="primary"].copy()
    d["r_char"]=[max(abs(np.linalg.eigvals(jacobian_block(-3,np.array([r.D_star,r.N_star,r.A_star,r.S_star]),params(r)))))
        for _,r in d.iterrows()]
    u=d[d.lambda_pattern>0].copy()
    u["ratio"]=u.lambda_DNA/u.lambda_pattern
    u["absolute_gap"]=u.lambda_DNA-u.lambda_pattern
    u["symmetric_error"]=2*u.absolute_gap/(u.lambda_DNA+u.lambda_pattern)
    u["relative_margin"]=u.lambda_pattern/u.r_char
    rows=[]
    for kind,cuts in (("relative_margin",[0.,1e-4,1e-3,1e-2,.05]),("lambda_pattern",[1e-3,1e-2,.1])):
        for cut in cuts:
            g=u[u[kind]>=cut]
            rows.append(dict(cut_kind=kind,cut=cut,n=len(g),ratio_median=float(g.ratio.median()),
                gap_median=float(g.absolute_gap.median()),symmetric_error_median=float(g.symmetric_error.median()),
                ratio_q25=float(g.ratio.quantile(.25)),ratio_q75=float(g.ratio.quantile(.75))))
    pd.DataFrame(rows).to_csv(OUT/"growth_boundary_strata.csv",index=False)
    u[["seed","sample_id","lambda_pattern","lambda_DNA","r_char","ratio","absolute_gap","symmetric_error","relative_margin"]].to_csv(OUT/"growth_error_points.csv.gz",index=False,compression="gzip")
    def tri(x):return np.where(x<=3,(x-1)/5,(6-x)/7.5)
    hill_weight=25*tri(d.n_d.to_numpy())*tri(d.n_ASC.to_numpy())
    us=(np.log(d.lam_ASC/P0.lam_ASC)+np.log(10))/(2*np.log(10))
    ud=(np.log(d.d_ASC/P0.d_ASC)+np.log(10))/(2*np.log(10))
    z1,z2=ndtri(us),ndtri(ud);rho=.75
    copula=np.exp((2*rho*z1*z2-rho*rho*(z1*z1+z2*z2))/(2*(1-rho*rho)))/np.sqrt(1-rho*rho)
    priors={"reference":np.ones(len(d)),"triangular_Hill":hill_weight,"correlated_S_turnover":np.asarray(copula)}
    summaries=[];reps=[]
    unstable=d.lambda_pattern.to_numpy()>0
    for name,w in priors.items():
        for seed,g in d.groupby("seed"):
            ix=d.index.get_indexer(g.index);ww=w[ix]
            reps.append(dict(prior=name,seed=int(seed),fraction=float(np.dot(ww,g.lambda_pattern>0)/ww.sum()),ESS=float(ww.sum()**2/np.dot(ww,ww)),mean_weight=float(ww.mean())))
        f=np.array([r["fraction"] for r in reps if r["prior"]==name]);half=t.ppf(.975,9)*f.std(ddof=1)/np.sqrt(10)
        summaries.append(dict(prior=name,fraction=float(f.mean()),mc95_low=float(f.mean()-half),mc95_high=float(f.mean()+half),
            ESS=float(w.sum()**2/np.dot(w,w)),weight_mean=float(w.mean()),
            unstable_growth_ratio_median=weighted_median((d.lambda_DNA/d.lambda_pattern).to_numpy()[unstable],w[unstable])))
    pd.DataFrame(reps).to_csv(OUT/"prior_replicates.csv",index=False)
    pd.DataFrame(summaries).to_csv(OUT/"prior_sensitivity.csv",index=False)
    return dict(primary_points=len(d),growth_strata=rows,prior_summaries=summaries)


def mechanism(reference,inputs):
    rows=[]
    for sid in reference.sample_id.astype(int):
        p=params(inputs.loc[sid]);fp=homogeneous_fixed_point(p)
        T=lag_similarity(fp,p);inv=np.diag(1/np.diag(T))
        Bs=mode_blocks(fp,p)
        expected=np.einsum("ij,mjk,kl->mil",T,Bs,inv)
        alternative=np.stack([production_lag_block(mu,fp,p) for mu in UNIQUE_MU])
        error=float(np.max(abs(alternative-expected))/max(1,np.max(abs(expected))))
        # Compare characteristic coefficients; eigenvalue matching near repeated roots is ill-conditioned.
        rows.append(dict(sample_id=sid,similarity_relative_residual=error,kappa_D=T[3,3]))
    pd.DataFrame(rows).to_csv(OUT/"lag_similarity_checks.csv",index=False)
    p=dataclasses.replace(P0,lam_N=2.);fp=homogeneous_fixed_point(p)
    tgrid=np.linspace(0,80,401);steps=[]
    S0=fp[3];U0=float(p.lam_D*hill_D(S0,p))
    for factor in (1.001,1.2,2.):
        Afinal=fp[2]*factor;S1=float(p.lam_ASC/p.d_ASC*hill_ASC(Afinal,p));U1=float(p.lam_D*hill_D(S1,p))
        S=S1+(S0-S1)*np.exp(-p.d_ASC*tgrid)
        Ufull=p.lam_D*hill_D(S,p);Ualt=U1+(U0-U1)*np.exp(-p.d_ASC*tgrid)
        for tt,a,b in zip(tgrid,Ufull,Ualt):
            steps.append(dict(A_step_factor=factor,time=tt,full_production=a,lag_production=b,
                full_normalized=(a-U0)/(U1-U0),lag_normalized=(b-U0)/(U1-U0)))
    step=pd.DataFrame(steps);step.to_csv(OUT/"mechanism_step_responses.csv",index=False)
    forced=[]
    fp0=homogeneous_fixed_point(P0);T=lag_similarity(fp0,P0)
    B=np.array([1.,0.,0.,0.])
    for mu in (-3.,0.,6.):
        J=jacobian_block(mu,fp0,P0);L=production_lag_block(mu,fp0,P0)
        for w in np.logspace(-4,3,201):
            a=np.linalg.solve(1j*w*np.eye(4)-J,B)[:3]
            b=np.linalg.solve(1j*w*np.eye(4)-L,T@B)[:3]
            forced.append(dict(mu=mu,omega=w,relative_response_error=float(np.linalg.norm(a-b)/np.linalg.norm(a))))
    pd.DataFrame(forced).to_csv(OUT/"mechanism_frequency_checks.csv",index=False)
    # Different Hill thresholds can be represented by a change of the free S concentration scale.
    rng=np.random.Generator(np.random.PCG64(2026102101));rescale=[]
    for i in range(100):
        KA,KS=np.exp(rng.uniform(np.log(.03),np.log(3),2));factor=KA/KS
        X=np.exp(rng.uniform(np.log(.03),np.log(3),4));D,N,A,S=X
        p0=dataclasses.replace(P0,Kd=KA)
        fq=np.array(drift(D,N,A,S,p0,6*D,6*N))
        fq[0]+=p0.lam_D*(expit(p0.n_d*(np.log(S)-np.log(KS)))-hill_D(S,p0))
        pp=dataclasses.replace(p0,lam_ASC=factor*p0.lam_ASC)
        fs=np.array(drift(D,N,A,factor*S,pp,6*D,6*N))
        target=fq*np.array([1,1,1,factor])
        rescale.append(dict(case=i,K_A=KA,K_S=KS,relative_error=float(np.linalg.norm(fs-target)/max(1,np.linalg.norm(target)))))
    pd.DataFrame(rescale).to_csv(OUT/"separate_Hill_threshold_checks.csv",index=False)
    return dict(mode_similarity_checks=len(rows)*len(UNIQUE_MU),max_similarity_residual=max(r["similarity_relative_residual"] for r in rows),
        max_DNA_response_error=max(r["relative_response_error"] for r in forced),
        step_max_normalized_difference={str(f):float(np.max(abs(g.full_normalized-g.lag_normalized))) for f,g in step.groupby("A_step_factor")},
        max_threshold_rescaling_error=max(r["relative_error"] for r in rescale))


def main():
    freeze();inputs=pd.read_csv(ROOT/"figures/data/global_sample_parameters.csv",float_precision="round_trip").set_index("sample_id")
    reference=pd.read_csv(ROOT/"figures/data/global_sample_results.csv")
    summary={"kernel_sha256":digest(__file__),"input_sha256":sha(ROOT/"figures/data/global_sample_parameters.csv")}
    for name,fun,args in (("qss",qss_analysis,(reference,inputs)),("recalibration",recalibration,()),
            ("bifurcations",bifurcations,(reference,inputs)),("growth_priors",growth_and_priors,()),("mechanism",mechanism,(reference,inputs))):
        print(name,flush=True);summary[name]=fun(*args);print(json.dumps(summary[name],allow_nan=True),flush=True)
    (OUT/"linear_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":main()
