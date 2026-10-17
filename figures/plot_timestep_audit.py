"""Figure S6: actual step refinement and solver-dependent trajectory errors."""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,NullLocator
import numpy as np

from figures.style import (BLUE,ORANGE,AQUA,YELLOW,INK_2,MUTED,apply_style,panel_label,save)
from scripts.timestep_audit import CONFIGS
from scripts.write_timestep_report import verified_data


def main():
    a,g,c,s,summary=verified_data()
    apply_style()
    fig,axes=plt.subplots(2,2,figsize=(8.7,7.2),gridspec_kw={"wspace":.43,"hspace":.62})
    ax=axes[0,0]
    wrong=s.assign(wrong=s.continuous_rhp!=s.discrete_unstable).groupby("dt").wrong.sum().sort_index()
    ax.plot(wrong.index,wrong.values,color=ORANGE,marker="s")
    ax.set_xscale("log")
    ax.set_xticks([.01,.02,.04,.08],["0.01","0.02","0.04","0.08"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(-6,122)
    ax.set_yticks([0,27,54,81,108])
    ax.text(.055,.72,"P0: false instability\nat h=0.04 and 0.08",transform=ax.transAxes,
            color=INK_2,fontsize=7.5)
    ax.set_xlabel("Fixed RK4 step h (dimensionless)")
    ax.set_ylabel("Mode blocks with wrong RHP count")
    ax.set_title("Coarse RK4 creates numerical instability",loc="left")
    ax=axes[0,1]
    for (point,model),color,marker in zip((("P2","full"),("P2","DNA"),("P3","full"),("P3","DNA")),
                                         (BLUE,ORANGE,AQUA,YELLOW),("o","s","^","D")):
        group=g[(g.point==point)&(g.model==model)].sort_values("dt")
        ax.loglog(group.dt,group.max_relative_D_perturbation_error,color=color,marker=marker,label=f"{point} {model}")
    ax.axhline(1e-4,color=MUTED,lw=.9,ls=":")
    ax.set_xticks([.01,.02,.04,.08],["0.01","0.02","0.04","0.08"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Fixed RK4 step h (dimensionless)")
    ax.set_ylabel("Max. relative Delta perturbation error")
    ax.set_title("Nonlinear early growth vs Radau",loc="left")
    ax.legend(loc="upper left",bbox_to_anchor=(0,.92),fontsize=7)
    names=[x[0] for x in CONFIGS]
    labels=["BDF\n"+r"$h\leq0.5$","BDF\n"+r"$h\leq0.25$","BDF\n"+r"$h\leq0.125$","BDF\ntighter","Radau\ntighter"]
    for ax,field,threshold,title,ylabel in (
        (axes[1,0],"events",.1,"Pattern emergence times","Relative event-time change (%)"),
        (axes[1,1],"final_relative_L2_all",1e-4,"Fixed-horizon endpoint fields","Relative L2 difference, all species")):
        for i,name in enumerate(names):
            group=a[a.config==name]
            values=100*group[["onset_01_relative_change","onset_05_relative_change"]].max(axis=1) if field=="events" else group[field]
            color=AQUA if name.startswith("Radau") else BLUE
            ax.scatter(i+np.linspace(-.11,.11,len(group)),values,s=18,color=color,alpha=.65)
        ax.axhline(threshold,color=MUTED,lw=.9,ls=":")
        if field=="events":
            ax.set_yscale("symlog",linthresh=1e-6,linscale=.5)
            ax.set_yticks([0,1e-6,1e-4,1e-2,.1])
        else:
            ax.set_yscale("log")
        ax.set_xticks(range(5),labels,fontsize=7)
        ax.set_xlim(-.45,4.45)
        ax.set_title(title,loc="left")
        ax.set_ylabel(ylabel)
    for ax,letter in zip(axes.flat,"ABCD"):
        panel_label(ax,letter,x=-.22)
    save(fig,"FigS6_timestep_audit",supplementary=True)


if __name__=="__main__":
    main()
