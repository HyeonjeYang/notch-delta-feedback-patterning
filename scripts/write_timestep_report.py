"""Generate Supplementary Note 4, retaining any failed accuracy assessments."""
from __future__ import annotations

import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

import numpy as np
import pandas as pd
from scripts.timestep_audit import OUT,CONFIGS,assess,kernel_hash,lattice_cases,sha


def verified_data():
    summary=json.loads((OUT/"summary.json").read_text())
    if summary["kernel_sha256"]!=kernel_hash():
        raise RuntimeError("time-step audit numerical sources have changed")
    a=pd.read_csv(OUT/"adaptive_runs.csv")
    if len(a)!=70 or a.duplicated(["point","amplitude","model","downstream","config"]).any():
        raise RuntimeError("incomplete adaptive audit")
    records=[]
    for spec in lattice_cases():
        from scripts.initial_condition_robustness import case_name
        for config in CONFIGS:
            path=OUT/"cases"/(case_name(spec)+"_"+config[0]+".npz")
            meta=json.loads(path.with_suffix(".json").read_text())
            if meta["kernel_sha256"]!=summary["kernel_sha256"] or meta["npz_sha256"]!=sha(path):
                raise RuntimeError(f"changed audit case: {path}")
            records.append(meta["summary"])
    pd.testing.assert_frame_equal(a,pd.DataFrame(records),check_dtype=False,rtol=1e-12,atol=1e-15)
    g=pd.read_csv(OUT/"fixed_growth.csv")
    c=pd.read_csv(OUT/"output_cadence.csv")
    s=pd.read_csv(OUT/"spectral_checks.csv")
    recomputed=assess(a,g,c,s)
    if recomputed["assessments"]!=summary["assessments"]:
        raise RuntimeError("audit assessment table disagrees with data")
    return a,g,c,s,summary


def main():
    a,g,c,s,summary=verified_data()
    text=["# Supplementary Note 4. Time-discretisation and solver audit\n",
        "The reported stationary threshold, mode spectra and instability fractions are computed "
        "from scalar equilibrium bracketing and the continuous-time Jacobian, not from stepping "
        "a trajectory forward. A temporal integration step therefore cannot create those "
        "particular classifications. Floating-point eigensolver errors remain a separate issue: "
        "selected cases were previously checked at 80-digit precision (Supplementary Note 2). "
        "This note tests whether time-discrete dynamics reproduce the spectra and whether "
        "the nonlinear growth, emergence times and endpoint fields depend on solver choices.\n",
        "**What dt means.** In the reference DOP853 runs and random-start BDF runs, the "
        "recording interval is an output sampling interval. It does not set the integrator's "
        "internal step. SciPy chooses adaptive internal steps using rtol and atol, while "
        "max_step imposes an upper bound. Event times are found from the solver's continuous "
        "interpolant. These distinctions follow the [SciPy solve_ivp documentation]"
        "(https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html); "
        "the installed version used for the audit is 1.16.3. Four checks are kept separate below.\n",
        "## 1. Genuine fixed-step refinement\n",
        "Classical RK4 is implemented independently of SciPy's adaptive integrators, without "
        "clipping. The full 24x24 lattice, rather than a three-cell symmetry reduction, is "
        "integrated with fixed dt=0.08, 0.04, 0.02 and 0.01 at P2/P3 for full/DNA models "
        "(16 trajectories). Each starts on its own exact K eigenvector at Delta amplitude "
        "1e-6 D*. Solutions are compared on the common 0.08 output grid to Radau at "
        "rtol=1e-12, atol=1e-15 and max_step=0.01. The fitting interval uses the same "
        "times for every integrator, ending at the last common time no greater than "
        "log(1000)/the predicted eigenvalue. Thus a shift of a sampled threshold cannot "
        "masquerade as a time-step effect. Finite-amplitude terms still cause a small "
        "growth-slope discrepancy from the exactly linear eigenvalue.\n",
        "**Table S4.** Fixed-step nonlinear lattice checks. Delta-field errors are relative "
        "to the growing Delta perturbation, not to the much larger background concentration.\n",
        "| point | model | dt | fitted growth | relative error vs eigenvalue | max relative Delta perturbation error vs Radau |",
        "|---|---|---|---|---|---|"]
    for row in g.itertuples():
        text.append(f"| {row.point} | {row.model} | {row.dt:g} | {row.fitted_slope:.9f} | {row.relative_slope_error:.3g} | {row.max_relative_D_perturbation_error:.3g} |")
    fine=g[g.dt==.01]
    window_rows=[]
    for (point,model),group in g.groupby(["point","model"]):
        rate=float(group.predicted.iloc[0])
        z=np.load(OUT/"cases"/f"RK4_{point}_{model}_h0.08.npz")
        times,amp=z["times"],z["reference_amplitude"]
        for upper in (1e-5,1e-4,1e-3):
            mask=times<=np.log(upper/1e-6)/rate
            slope=float(np.polyfit(times[mask],np.log(amp[mask]),1)[0])
            window_rows.append(dict(point=point,model=model,upper_relative_amplitude=upper,
                window_end=float(times[mask][-1]),fitted_slope=slope,predicted=rate,
                relative_slope_error=abs(slope/rate-1)))
    windows=pd.DataFrame(window_rows)
    windows.to_csv(OUT/"amplitude_window_diagnostic.csv",index=False)
    text.extend([f"\nAt dt=0.01, the largest Delta perturbation error is "
        f"{fine.max_relative_D_perturbation_error.max():.3g}; the largest fitted-growth error "
        f"against the continuous eigenvalue is {fine.relative_slope_error.max():.3g}. "
        "These are early-growth checks and do not validate fixed-step RK4 for the stiff "
        "late nonlinear regime.\n",
        "An additional post-audit diagnostic reuses the stored high-accuracy reference "
        "amplitudes, without new trajectories or changing the prespecified assessments. "
        "Reducing the upper amplitude of the common theoretical window from 1e-3 to "
        f"1e-4 and 1e-5 reduces the largest relative growth-slope discrepancy from "
        f"{windows[windows.upper_relative_amplitude==1e-3].relative_slope_error.max():.3g} "
        f"to {windows[windows.upper_relative_amplitude==1e-4].relative_slope_error.max():.3g} "
        f"and {windows[windows.upper_relative_amplitude==1e-5].relative_slope_error.max():.3g}. "
        "This supports finite-amplitude bias as the explanation for the residual slope "
        "difference once temporal error is small. The reference manuscript's reported "
        "fits and original window rule are retained; their last digits should not be "
        "interpreted as experimentally calibrated or exact accuracy.\n",
        "A separate linear test uses every one of the 54 distinct graph eigenvalues at "
        "P0/P2/P3 for full and DNA models. At each fixed dt the RK4 stability polynomial "
        "R(hJ)=I+hJ+(hJ)^2/2+(hJ)^3/6+(hJ)^4/24 is compared with the exact matrix "
        "exponential at t=2. The number of unstable discrete eigenvalues is compared with "
        "the continuous RHP count; graph modes are not inferred from a single spatial start. "
        f"There are {len(s)} mode/model/step comparisons and "
        f"{summary['n_discrete_spectral_disagreements']} count disagreements. This comparison "
        "is a check on the chosen steps at these operating points, not a statement that any "
        "explicit method at arbitrarily large dt preserves stability.\n",
        "The 216 disagreements are all at the continuously stable P0: h=0.04 and 0.08 "
        "each make all 54 mode blocks unstable in both models. These are genuine numerical "
        "instabilities of the coarse RK4 discretisation. All 648 comparisons at h=0.02 "
        "and 0.01 agree with the continuous RHP counts. The largest fine-step propagator "
        f"errors are {summary['spectral_max_propagator_error_by_dt']['0.02']:.3g} and "
        f"{summary['spectral_max_propagator_error_by_dt']['0.01']:.3g}, respectively. "
        "The failed coarse controls were not used in the adaptive production calculations.\n",
        "## 2. Internal step caps, tolerance refinement and an independent stiff solver\n",
        "Fourteen fixed initial states cover the first seed at every amplitude and both "
        "models/operating points (12 main cases), plus the two independent-S controls at "
        "amplitude 1e-4. For each case the existing BDF trajectory (rtol=1e-9, atol=1e-11, "
        "unbounded max_step) is compared with three step caps, 0.5, 0.25 and 0.125; "
        "BDF with cap 0.125 and tighter tolerances 1e-11/1e-14; and Radau with cap 0.25 "
        "and the same tighter tolerances. All 70 new trajectories run to the original fixed "
        "horizon, t=600 at P2 and t=300 at P3. Accepted internal step meshes are deposited. "
        "Cellwise endpoint differences use the original cell indices, without translation "
        "alignment, making this a stricter comparison than the full/DNA pattern comparison.\n",
        "**Table S5.** Largest deviations across the five new configurations for each fixed "
        "initial state. Event change takes the maximum over CV(D)=0.1 and 0.5.\n",
        "| point | model | initial S | amplitude | max relative event change | max relative endpoint L2, all species | max scaled endpoint drift |",
        "|---|---|---|---|---|---|---|"])
    for keys,group in a.groupby(["point","model","downstream","amplitude"]):
        point,model,down,amp=keys
        event=group[["onset_01_relative_change","onset_05_relative_change"]].to_numpy().max()
        text.append(f"| {point} | {model} | {down} | {amp:.0e} | {event:.3g} | {group.final_relative_L2_all.max():.3g} | {group.final_scaled_rhs.max():.3g} |")
    text.extend([f"\nAcross these cases the largest relative event-time change is "
        f"{100*summary['max_adaptive_event_relative_change']:.5f}%, the largest all-species "
        f"endpoint relative L2 difference is {summary['max_final_relative_L2_all']:.3g}, "
        f"and the largest absolute K-power-fraction change is {summary['max_k_fraction_abs_change']:.3g}. "
        f"The largest endpoint scaled drift is {summary['max_final_scaled_rhs']:.3g}. "
        "Endpoint agreement under refinement is numerical evidence of reproducibility for "
        "these starts; it does not prove attractor stability or equality of basins.\n",
        "## 3. Output sampling without changing the integration\n",
        "The identical dense Radau growth solution is resampled at recording intervals "
        "0.5, 0.25 and 0.125. The original first-recorded-amplitude rule is applied separately "
        "to each sampled curve. Its endpoint may change with cadence, although the "
        "underlying numerical trajectory and its internal steps are identical.\n",
        "**Table S6.** Growth fits under the original output-based window rule.\n",
        "| point | model | output interval | window end | fitted growth |",
        "|---|---|---|---|---|"])
    for row in c.itertuples():
        text.append(f"| {row.point} | {row.model} | {row.record_interval:g} | {row.window_end:g} | {row.fitted_slope:.9f} |")
    text.extend([f"\nThe largest relative slope change from the finest recording interval is "
        f"{summary['max_cadence_relative_slope_change']:.3g}. The output interval must therefore "
        "be reported for fitted slopes; it should not be described as the adaptive solver dt.\n",
        "## 4. Prespecified assessments and limits\n",
        "The protocol and numerical assessment thresholds were written before the new runs. "
        "Every assessment is reported below, including a failure if a coarse step exceeds "
        "the chosen accuracy target. A failed accuracy target need not mean an artificial "
        "instability; convergence at finer steps must be assessed separately.\n",
        "| assessment | within prespecified target |","|---|---|"])
    for key,value in summary["assessments"].items():
        text.append(f"| {key} | {'yes' if value else 'no: inspect the retained step-level results'} |")
    coarse_field=g[g.max_relative_D_perturbation_error>=1e-4]
    fine_spectrum=s[s.dt<=.02]
    paired=a[a.downstream=="qss"].pivot(index=["point","amplitude","config"],columns="model",values="onset_01")
    interpretation=dict(
        adaptive_production_checks_pass=all(summary["assessments"][k] for k in
            ("adaptive_events","adaptive_fields","adaptive_k_fraction","adaptive_stationarity")),
        fixed_fields_at_h_001_pass=bool((fine.max_relative_D_perturbation_error<1e-4).all()),
        finer_spectral_counts_agree=bool((fine_spectrum.continuous_rhp==fine_spectrum.discrete_unstable).all()),
        coarse_fixed_field_accuracy_failures=coarse_field[["point","model","dt","max_relative_D_perturbation_error"]].to_dict("records"),
        reintegrated_main_pairs=len(paired),full_later_in_reintegrated_pairs=int((paired["full"]>paired["DNA"]).sum()),
        original_all_assessments_pass=summary["all_assessments_pass"],
        explanation="Original coarse-step failures are retained; they do not describe the adaptive production settings. No target was relaxed.")
    (OUT/"assessment_interpretation.json").write_text(json.dumps(interpretation,indent=2)+"\n",encoding="utf-8")
    text.extend(["\nThe combined all_assessments_pass flag is false because the sweep includes "
        "failed coarse RK4 controls. All four adaptive production-comparison targets pass. "
        "The fixed-field target is exceeded only by P3 DNA at h=0.08 "
        f"(relative perturbation error {coarse_field.max_relative_D_perturbation_error.max():.3g}); "
        "it passes at every smaller tested step. The spectral failures are the coarse P0 "
        "cases described above. No failed control was discarded and no accuracy target "
        "was relaxed. This distinction is also recorded in assessment_interpretation.json. "
        f"Across the {len(paired)} reintegrated matched model pairs, the full model remains "
        f"later at CV=0.1 in all {int((paired['full']>paired['DNA']).sum())} comparisons.\n"])
    text.extend([f"\nThe minimum accepted adaptive concentration or RK4 stage concentration is "
        f"{summary['min_accepted_or_RK_stage_concentration']:.3g}. No positive concentration "
        "floor or post-step clipping is used. The fixed-step implementation is independently "
        "tested for fourth-order convergence against the nonlinear exact solution "
        "y'=y^2, y(0)=1. Tests supplement, and do not replace, the lattice comparisons.\n",
        "This audit covers fourteen representative lattice initial states, not every seed or "
        "parameter vector. The 24x24 contact graph is a discrete biological modelling choice; "
        "this is not a spatial mesh-convergence study or a thermodynamic-limit result. "
        "Agreement among time integrators does not establish tissue validity, remove model "
        "misspecification or mathematically exclude all floating-point error. The deposited "
        "data support only the stated resolution and solver comparisons.\n",
        "![**Figure S6. Time-step and solver checks.** (A) Number of mode blocks with an "
        "incorrect discrete RHP count, pooled over P0/P2/P3 and both models, versus fixed "
        "RK4 step. All failures occur at P0. (B) Maximum relative Delta perturbation error versus the high-accuracy Radau "
        "reference in the nonlinear early-growth checks. (C) Relative event-time changes "
        "from the original BDF run, taking the larger of the CV=0.1 and 0.5 differences; "
        "each point is one of fourteen matched initial states. (D) Raw all-species endpoint "
        "relative L2 difference from the original run, without spatial alignment, at the "
        "fixed horizons. Panel C uses a symmetric logarithmic scale with a linear region "
        "at zero. Horizontal dotted lines show prespecified accuracy targets, not "
        "biological confidence intervals.](../figures/supplementary/FigS6_timestep_audit.png){width=100%}\n"])
    path=ROOT/"manuscript/timestep_supplement.md"
    path.write_text("\n".join(text)+"\n",encoding="utf-8")
    print(f"wrote {path}")


if __name__=="__main__":
    main()
