"""Main Figure 6 and one supplementary figure per additional reviewer question."""
from __future__ import annotations
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figures.style import apply_style,save,panel_label,BLUE,ORANGE,AQUA,YELLOW,INK_2
R=ROOT/"figures/data/reviewer_extensions"
COLORS={"full":BLUE,"DNA":ORANGE}


def read(name):return pd.read_csv(R/(name+".csv"))
def layout(nrows=2,ncols=2,size=(8.2,6.2)):
    f,a=plt.subplots(nrows,ncols,figsize=size,layout="constrained")
    for ax,letter in zip(np.asarray(a).flat,"ABCDEFGHI"):panel_label(ax,letter,x=-.14,y=1.02)
    return f,np.asarray(a).ravel()


def qss_curve(ax):
    q=read("qss_validity")
    for name,c in (("P2",BLUE),("P3",ORANGE)):
        g=q[q.family==name];ax.loglog(g.alpha,g.relative_growth_error,color=c,label=name)
    ax.axhline(.1,color=INK_2,ls="--",lw=1,label="10% error")
    ax.set(xlabel=r"compensated turnover $\alpha$",ylabel="relative growth error",title="When is QSS accurate?")
    ax.legend(fontsize=7)


def clock_curve(ax):
    g=read("clock_recalibration")
    ax.plot(g.lam_N,g.full_growth,color=BLUE,label="full")
    ax.plot(g.lam_N,g.calibrated_growth,color=ORANGE,ls="--",label="DNA, clock fitted at P2")
    for x,name in ((2.,"P2"),(5.623,"P3")):
        row=g[np.isclose(g.lam_N,x)].iloc[0];ax.plot(x,row.full_growth,"o",color=BLUE,ms=4)
        ax.annotate(name,(x,row.full_growth),xytext=(3,7),textcoords="offset points",fontsize=7)
    ax.text(.97,.12,"P3 prediction: −16.5%",transform=ax.transAxes,ha="right",fontsize=8,color=ORANGE)
    ax.set(xlabel=r"Notch production $\lambda_N$",ylabel=r"growth rate at $\mu=-3$",title="A fitted clock does not transfer exactly")
    ax.legend(loc="upper left",fontsize=7)


def pole_path(ax,sid=679):
    g=pd.read_csv(R/"bifurcation_poles.csv.gz");g=g[(g.sample_id==sid)&(g.alpha>=.02)&(g.alpha<=2)]
    # All four K-block poles; the fast stable poles are outside the display.
    ax.scatter(g.alpha,g.re,s=5,color=BLUE,alpha=.7,edgecolors="none",label="real part",rasterized=True)
    ax.scatter(g.alpha,abs(g.im),s=4,color=ORANGE,alpha=.7,edgecolors="none",label="|imaginary part|",rasterized=True)
    b=read("oscillatory_boundaries");b=b[b.sample_id==sid].iloc[0]
    ax.axvline(b.alpha,color=AQUA,ls="--",lw=1)
    c=read("real_pair_collisions");c=c[c.sample_id==sid]
    if len(c):ax.axvline(c.alpha.iloc[0],color=YELLOW,ls=":",lw=1.4)
    ax.axhline(0,color=INK_2,lw=.7);ax.set_xscale("log");ax.set_ylim(-.12,.35)
    ax.set(xlabel=r"compensated turnover $\alpha$",ylabel="K-mode poles",title="Oscillatory entry, then real unstable pair" if len(c) else "Oscillatory entry at unchanged static gain")
    ax.text(.03,.97,f"sample {sid}; G = {b.G_K:.3f} throughout",transform=ax.transAxes,va="top",fontsize=7)
    ax.legend(loc="lower right",fontsize=6.7)


def step_curve(ax,factor=2.):
    g=read("mechanism_step_responses");g=g[np.isclose(g.A_step_factor,factor)]
    ax.plot(g.time,g.full_normalized,color=BLUE,label="explicit proneural")
    ax.plot(g.time,g.lag_normalized,color=ORANGE,ls="--",label="production lag")
    diff=abs(g.full_normalized-g.lag_normalized).max()
    ax.set(xlim=(0,60),xlabel="time (dimensionless)",ylabel="fraction of final production change",title="Same linear response, distinct large step" if factor>1.01 else "Small step approaches linear equivalence")
    ax.text(.97,.12,f"clamped A: {factor:g}×\nmax difference {diff:.3%}" if factor<1.01 else f"clamped A: {factor:g}×\nmax difference {diff:.1%}",transform=ax.transAxes,ha="right",fontsize=8)
    ax.legend(fontsize=7.2,loc="center right",bbox_to_anchor=(1,.52))


def main_figure():
    f,a=layout();qss_curve(a[0]);clock_curve(a[1]);pole_path(a[2]);step_curve(a[3])
    save(f,"Fig6_dynamic_information")


def supplement_qss():
    q=read("qss_validity");f,a=layout();qss_curve(a[0])
    positive=q.dropna(subset=["relative_growth_error"])
    sc=a[1].scatter(positive.pole_over_S_rate,positive.relative_growth_error,c=np.log10(positive.eigenvalue_condition),cmap="viridis",s=7,alpha=.6,edgecolors="none",rasterized=True)
    a[1].set(xscale="log",yscale="log",xlabel=r"$|s_{\rm lead}|/(\alpha d_S)$",ylabel="relative growth error",title="A slow pole alone is insufficient")
    a[1].axhline(.1,color=INK_2,ls="--",lw=1);a[1].axvline(.01,color=INK_2,ls=":",lw=1)
    f.colorbar(sc,ax=a[1],label="log10 eigenvalue condition",shrink=.8)
    t=read("qss_accuracy_thresholds")
    for target,c in ((.05,BLUE),(.1,ORANGE)):
        g=t[t.target==target].sort_values("alpha_required");a[2].plot(np.arange(1,len(g)+1),g.alpha_required,color=c,label=f"{target:.0%} error")
    a[2].set(yscale="log",xlabel="family rank (separately sorted)",ylabel="first sustained passing α on grid",title="Required turnover is family dependent");a[2].legend()
    for _,g in q.groupby("family"):
        a[3].loglog(g.alpha,g.relative_growth_error,color=BLUE,alpha=.3,lw=.8)
    a[3].axhline(.1,color=ORANGE,ls="--",lw=1);a[3].set(xlabel=r"$\alpha$",ylabel="relative growth error",title="All 52 matched families")
    save(f,"FigS7_qss_validity",True)


def supplement_clock():
    f,a=layout();clock_curve(a[0]);g=read("clock_recalibration")
    a[1].plot(g.lam_N,100*g.calibrated_relative_error,color=ORANGE,label="P2 fit")
    a[1].plot(g.lam_N,100*g.LS_relative_error,color=AQUA,label="training-range LS fit")
    a[1].axvspan(1.5,3,color=AQUA,alpha=.08);a[1].axvspan(3.5,6.5,color=BLUE,alpha=.06)
    a[1].axhline(0,color=INK_2,lw=.7);a[1].set(xlabel=r"$\lambda_N$",ylabel="growth prediction error (%)",title="Training and held-out conditions");a[1].legend()
    m=read("clock_mode_predictions")
    a[2].plot(m.mu,m.full_growth,color=BLUE,label="full");a[2].plot(m.mu,m.calibrated_growth,color=ORANGE,ls="--",label="P2 clock × DNA")
    a[2].axhline(0,color=INK_2,lw=.7);a[2].set(xlabel=r"graph eigenvalue $\mu$",ylabel="leading real part",title="P2: one clock, all spatial modes");a[2].legend()
    h=read("clock_frequency_response");h=h[h.omega>0]
    for prefix,c,l in (("full",BLUE,"full"),("calibrated",ORANGE,"clock × DNA")):
        a[3].loglog(h.omega,abs(h[prefix+"_re"]+1j*h[prefix+"_im"]),color=c,label=l)
    a[3].set(xlabel=r"frequency $\omega$",ylabel="Delta response magnitude",title="Stable P0, Delta-production forcing");a[3].legend()
    save(f,"FigS8_clock_recalibration",True)


def supplement_geometry():
    f,a=layout();t=read("geometry_stationary_thresholds")
    a[0].plot(t['size'],t.stationary_lam_N,"o-",color=BLUE);a[0].set(xticks=[23,24,25],xlabel="periodic lattice side",ylabel=r"stationary threshold $\lambda_N^c$",title="Commensurability shifts onset slightly")
    g=read("geometry_pairs");p=g[(g.condition=="periodic")&(g.point!="P0")]
    for point,c,dx in (("P2",BLUE,-.06),("P3",ORANGE,.06)):
        h=p[p.point==point];a[1].scatter(h['size']+dx,h.full_DNA_onset_ratio,color=c,s=25,label=point,alpha=.8)
    a[1].set(xticks=[23,24,25],xlabel="periodic lattice side",ylabel="full / DNA time to CV = 0.1",title="All 18 pairs retain slower full onset");a[1].legend()
    runs=read("geometry_runs");h=runs[(runs.condition=="periodic")&(runs.point!="P0")]
    for model,dx in (("full",-.07),("DNA",.07)):
        z=h[h.model==model];a[2].scatter(z['size']+dx,z.near_K_band_power,color=COLORS[model],s=23,alpha=.7,label=model)
    a[2].set(xticks=[23,24,25],ylim=(0,1),xlabel="periodic lattice side",ylabel=r"endpoint power at $\mu\leq-2.8$",title="Endpoint structure remains variable");a[2].legend()
    base=R/"geometry"
    for condition,c in (("periodic",BLUE),("open",ORANGE)):
        z=np.load(base/f"{condition}_P0_24_2026102001_full_cv0.0.npz")
        a[3].semilogy(z['times'],np.maximum(z['cv'],1e-15),color=c,label=condition)
    a[3].set(xlabel="time (dimensionless)",ylabel="Delta CV at P0",title="Open edges provide deterministic forcing");a[3].legend()
    save(f,"FigS9_geometry_robustness",True)


def supplement_bifurcation():
    f,a=layout();pole_path(a[0],1730);pole_path(a[1],679)
    b=read("oscillatory_boundaries");x=np.arange(len(b));labels=b.sample_id.astype(str)
    a[2].scatter(x,b.frequency,color=BLUE,s=27);a[2].set(yscale="log",xticks=x,xticklabels=labels,xlabel="reference sample ID",ylabel="crossing frequency",title="All 12 first crossings are oscillatory")
    a[2].tick_params(axis="x",rotation=60,labelsize=6.5)
    c=read("real_pair_collisions")
    a[3].scatter(b.alpha,b.G_K,color=BLUE,label="oscillatory boundary")
    a[3].scatter(c.alpha,c.G_K,color=ORANGE,marker="s",label="later real-pair collision")
    for _,r in c.iterrows():
        initial=b[b.sample_id==r.sample_id].iloc[0];a[3].plot([initial.alpha,r.alpha],[r.G_K,r.G_K],color=INK_2,lw=.6,alpha=.5)
    a[3].set(xscale="log",xlabel=r"turnover $\alpha$",ylabel="unchanged static G(−3)",title="G stays below one along every path");a[3].legend(fontsize=6.5)
    save(f,"FigS10_oscillatory_routes",True)


def supplement_patterns():
    g=read("pattern_stability");f,a=layout()
    for model,c in COLORS.items():
        h=g[g.model==model];a[0].scatter(h.index,h.leading_real,color=c,label=model,s=30)
    a[0].set(ylim=(-.1005,-.0965),xlabel="selected endpoint index",ylabel="leading real eigenvalue",title="Complete spectra: strictly negative leaders");a[0].legend()
    for i,r in g.iterrows():
        z=np.load(R/"pattern_stability"/(r['name']+".npz"))
        a[1].semilogy(z['times'],z['paired_distance']/z['paired_distance'][0],color=COLORS[r.model],alpha=.65,lw=1)
        if r['name']=="pure_P2_full":
            ev=z['poles'];a[2].scatter(ev.real,ev.imag,color=BLUE,s=6,alpha=.5,rasterized=True)
            a[2].text(.03,.9,f"max |Im| = {abs(ev.imag).max():.2g}\n(round-off scale)",transform=a[2].transAxes,fontsize=7,color=INK_2)
    a[1].set(xlabel="time after relative 10⁻⁴ perturbation",ylabel="paired distance / initial distance",title="Recovery relative to unperturbed continuation")
    a[2].axvline(0,color=ORANGE,ls="--",lw=1);a[2].set(ylim=(-.05,.05),xlabel="real part",ylabel="imaginary part",title="Pure-mode full endpoint: all poles")
    a[3].semilogy(g.index,g.newton_residual,"o",color=BLUE,label="scaled Newton residual")
    a[3].semilogy(g.index,g.sparse_eigenpair_residual,"s",color=ORANGE,label="normalised sparse eigenpair residual")
    a[3].set(xlabel="selected endpoint index",ylabel="numerical residual",title="Independent sparse eigenpair check");a[3].legend(fontsize=6.5)
    save(f,"FigS11_patterned_state_stability",True)


def supplement_priors():
    f,a=layout();g=read("prior_sensitivity");x=np.arange(len(g))
    a[0].errorbar(x,100*g.fraction,yerr=[100*(g.fraction-g.mc95_low),100*(g.mc95_high-g.fraction)],fmt="o",color=BLUE,capsize=4)
    a[0].set(xticks=x,xticklabels=["reference","triangular\nHill","correlated\nS turnover"],ylabel="instability fraction (%)",title="Illustrative prior sensitivity")
    s=read("growth_boundary_strata");s=s[s.cut_kind=="relative_margin"]
    a[1].errorbar(np.arange(len(s)),s.ratio_median,yerr=[s.ratio_median-s.ratio_q25,s.ratio_q75-s.ratio_median],fmt="o-",color=ORANGE,capsize=3)
    a[1].set(xticks=np.arange(len(s)),xticklabels=[f"{x:g}" for x in s.cut],xlabel="minimum relative growth margin",ylabel="DNA / full growth (median, IQR)",title="Bias persists away from the boundary")
    p=pd.read_csv(R/"growth_error_points.csv.gz")
    for ax,field,label in ((a[2],"absolute_gap","absolute growth-rate difference"),(a[3],"symmetric_error",r"$2(\lambda_{DNA}-\lambda_F)/(\lambda_{DNA}+\lambda_F)$")):
        ax.scatter(p.relative_margin,p[field],color=BLUE,s=3,alpha=.3,edgecolors="none",rasterized=True)
        ax.set(xscale="log",xlabel="full growth / maximum K-pole magnitude",ylabel=label)
    a[2].set_yscale("log");a[2].set_title("Absolute differences avoid ratio inflation");a[3].set_title("Bounded symmetric growth error")
    save(f,"FigS12_growth_and_prior_sensitivity",True)


def supplement_mechanism():
    f,a=layout();step_curve(a[0],1.001);step_curve(a[1],1.2);step_curve(a[2],2.)
    g=read("mechanism_frequency_checks")
    for mu,c in ((-3,BLUE),(0,ORANGE),(6,AQUA)):
        h=g[np.isclose(g.mu,mu)];a[3].loglog(h.omega,h.relative_response_error.clip(lower=1e-18),color=c,label=f"μ = {mu}")
    a[3].set(xlabel=r"frequency $\omega$",ylabel="relative D,N,A response difference",title="Linear observable equivalence at stable P0");a[3].legend()
    a[3].text(.97,.17,"display floor 10⁻¹⁸",transform=a[3].transAxes,ha="right",fontsize=7,color=INK_2)
    save(f,"FigS13_mechanism_identifiability",True)


def main():
    apply_style()
    for function in (main_figure,supplement_qss,supplement_clock,supplement_geometry,supplement_bifurcation,supplement_patterns,supplement_priors,supplement_mechanism):function()


if __name__=="__main__":main()
