# Audit and re-run procedure

This file gives the step-by-step procedure, from static checks to full recomputation, for the
results of the manuscript *Actuator-Aware Observer Inputs for Cascade ADRC of Quadrotors Under
Rotor Lag, Measurement Delay and Sensor Noise*. Steps 1–3 use only the stored results and take
minutes. Steps 4–6 recompute results and take hours (MATLAB) to days (PX4 SITL).

Requirements are listed in `README.md`. The analysis scripts (steps 2–4) need the Statistics and
Machine Learning Toolbox (`tinv`).

## 1. Integrity of the stored results

```matlab
cd analysis
M = results_manifest;   % recomputes SHA-256 of every result file listed in the manifest
```

Every entry of `results/analysis/results_manifest.csv` maps a reported table or figure to the
script that produced it, the result file and the design/gain source. A row with `MISSING` or a
changed hash means the file differs from the published version.

## 2. Controller equations against the implementation

```matlab
cd analysis
pid_equation_check       % PID baseline equations (manuscript Sec. II) vs. chan_output/chan_update_base
```

Writes `results/analysis/pid_equation_consistency.csv`; the expected maximum error is 0 for the
control, integrator and rate-filter states of every PID class and design.

## 3. Regeneration of tables and figures from stored results

```matlab
cd analysis
clean_regeneration_test  % copies scripts + stored results to an empty folder, regenerates, compares
regenerate_tables        % sitl_vs_matlab, stack_hover_table, stack_tracking_table (in place)
make_figures_setpoint    % budget-performance curve, noise activity, mismatch counts
```

`clean_regeneration_test` writes `results/analysis/clean_regeneration_test.txt` (expected:
all compared files identical). It starts a separate `matlab -batch` session with the MATLAB found
on the system path.

## 4. Static and linear checks of the designs

```matlab
cd MATLAB_simulation
run_design               % designs and audits (sensitivity, response, covariance ceilings)
```

and in `analysis/`: `settling_check` (response requirement of every design). The reported
gains are the best feasible candidates found by the stated finite search.

## 5. Full MATLAB recomputation

```matlab
cd MATLAB_simulation
run_pipeline             % = run_design, run_compute, finish_checks
```

`run_compute` runs the MEX equivalence checks, the paired Monte Carlo trials, the linear grids and
the hover runs; `finish_checks` runs the Simulink replication (`Simulink/`), output validation,
figures and export. Seeds are fixed (see `sim_config.m` and the manuscript appendix), so the
recomputed MATLAB metrics should agree with the stored ones to round-off. The redesigned
controllers (strict PID, commanded-input, tilt-input and estimator-aware designs, PID set-point
weights) are recomputed by the corresponding scripts in `analysis/`, as listed in the manifest.

## 6. PX4 software-in-the-loop

See "PX4 SITL" in `README.md`: build with `PX4_SITL/do_build.m`, then run the campaigns
(`campaign_main`, `campaign_activity`, `campaign_estimator`, `campaign_redesign`) with the fixed
run plans `results/analysis/run_plan*.csv`. The lockstep loop is not deterministic: repeated
campaigns reproduce the reported SITL statistics only approximately (manuscript Sec. V-J), and the
number of skipped lockstep steps depends on host load. Refitting the estimator-noise models needs
the raw logs of the work folder, which are not part of this package.
