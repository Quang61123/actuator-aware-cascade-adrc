# Cascade LADRC with actuator-aware observer inputs for quadrotors

Code, controller gains, random seeds and per-seed results for the manuscript
*Actuator-Aware Observer Inputs for Cascade ADRC of Quadrotors Under Rotor Lag, Measurement Delay and Sensor Noise*
(X.-Q. Dao, N.-B. Le, L.-M. Nguyen). Repository: https://github.com/Quang61123/actuator-aware-cascade-adrc

MATLAB/Simulink code for the design, simulation and PX4 software-in-the-loop (SITL) evaluation of
cascade linear active disturbance rejection control (LADRC) of a quadrotor, with actuator-aware
(AA) extended-state-observer inputs, budget-constrained covariance design, PID baselines and
estimator-aware (-E) redesigns for a PX4 + EKF2 stack.

## Requirements

- MATLAB with Control System Toolbox, Statistics and Machine Learning Toolbox, Optimization Toolbox
  and Parallel Computing Toolbox (optional; loops fall back to serial). Tested with MATLAB R2026b on Windows.
- Simulink (Simulink model only).
- MATLAB Coder and a C/C++ compiler to rebuild the MEX runners (`build_mex`, `build_extra_mex`).
  The supplied `.mexw64` files are Windows builds; on other platforms rebuild them or use the
  MATLAB fallback (slower).
- PX4 SITL only: UAV Toolbox Support Package for PX4 Autopilots, Embedded Coder, WSL with a PX4
  build (`PX4-Ubuntu-22.04`, PX4 v1.16.1) and Python 3 for the simulator relay.

## Folder layout

| Folder | Content |
|---|---|
| `MATLAB_simulation/` | Plant, sensor and actuator models, controller (`ctrl_law`), channel models and covariance tools, primary design pipeline, Monte Carlo trials, checks and figures. Configuration in `sim_config.m`. |
| `Simulink/` | Programmatic Simulink model (`create_simulink_model` -> `control_system.slx`) using the same controller and plant functions, and its validation. |
| `PX4_SITL/` | Code-generated PX4 controller (`create_sitl_controller`, `sitl_ctrl_step`), lockstep SITL runner (`sitl_run`), campaigns, log extraction and reports. `wsl/` holds the airframe, the launch script and the relay. |
| `estimator_noise/` | Estimator-error extraction from SITL logs, delay estimation, coloured-noise model fitting and estimator-aware predictions. |
| `analysis/` | Response-requirement checks, redesigns (strict PID, commanded-input and tilt-input ablations, estimator-stack designs, set-point weights), SITL/MATLAB summaries, table regeneration and the result manifest. |
| `results/` | Stored results (gains, CSV summaries, MAT files, per-run SITL metrics). |
| `figures/` | Created by the figure scripts. |

## Typical use

Primary design, simulations, checks and figures:

```matlab
cd MATLAB_simulation
run_pipeline            % = run_design, run_compute, finish_checks
```

Individual steps: `run_design` (design and audits), `run_compute` (MEX checks, Monte Carlo trials,
linear grids, hover runs), `finish_checks` (Simulink, validation, figures, export). Monte Carlo
results are checkpointed; a checkpoint is reused only if protocol and gains are unchanged.

Simulink model:

```matlab
cd Simulink
create_simulink_model
out = run_simulink(6, 'AA-Budget', 40);
```

Analysis scripts (MATLAB only, from stored results):

```matlab
cd analysis
regenerate_tables        % sitl_vs_matlab, stack_hover_table, stack_tracking_table
make_figures_setpoint    % budget-performance curve, noise activity, mismatch counts
results_manifest         % result item -> script -> file -> SHA-256
clean_regeneration_test  % regenerates the tables in an empty folder and compares
```

## PX4 SITL

`px4sitl_work_dir` returns the local work folder for logs, ulogs, build files and launch scripts
(default `%USERPROFILE%\px4sitl_work`; override with the environment variable `PX4SITL_WORK`).
Copy `PX4_SITL/wsl/*` into that folder. Build the controller with `PX4_SITL/do_build.m`, then run a
campaign (`campaign_main`, `campaign_activity`, `campaign_estimator`, `campaign_redesign`);
parallel workers are selected with `SITL_WORKER` / `SITL_NWORKERS`. Run plans in
`results/analysis/run_plan*.csv` are fixed before running and every listed run is reported.

Large raw logs (ulogs and controller-rate logs in the work folder) are not part of this package;
scripts that refit noise models from those logs (`estimator_noise/fit_noise_per_design`,
`analysis/fit_noise_models`, `validate_noise_models`, ...) need them. Their file prefixes follow
`run_sitl_estimator` (`est_`, `est_E_`, `est_cal_`, `yaw_ceiling_`) and `run_sitl_redesign`.

## Controller kinds and ids

`ctrl_law` kinds: 1/2 PID, 4 RO-LADRC, 7 AA RO-LADRC, 8 RO-HalfGain, 9 AA-HalfGain,
10 A4-CommandedTilt (commanded-tilt horizontal ESO input). PX4 controller ids are listed in
`create_sitl_controller.m` and `sitl_run.m`.

## Notes

- Seeds are paired across controllers and conditions; runs are not independent samples.
- Stored design structures and result tables keep the identifiers they were generated with
  (design names, run tags).

- Audit and full re-run procedure: `REPRODUCE.md`.

## Superseded rows in first-campaign result files

Some first-campaign files keep rows for designs that were redesigned later (response requirement, extended ratio range, PID set-point rule). The manuscript uses the redesigned results; `results/analysis/results_manifest.csv` maps every reported table and figure to its file.

| File | Superseded rows | Use instead |
|---|---|---|
| `results/summary.csv`, `results/trials.csv` | RO-HalfGain (E1, S5); PID-Budget-Strict (E0); PID-Budget S5 with the searched set-point weights | `analysis/matlab_ROx*.csv`; `analysis/matlab_PIDS.csv`; `analysis/setpoint_summary.csv`, `analysis/setpoint_trials.csv` |
| `results/paired.csv` | A1-Commanded (original grid) | `analysis/matlab_paired.csv` (A1x) |
| `results/analysis/tilt_summary.csv` | A4-CommandedTilt (retuned) rows | `analysis/tilt_scan_summary.csv` (retuned, full requirements) |

The PID-Budget, AA-Budget and AA-HalfGain rows of the first-campaign files remain valid (Table V).

## Citation

If you use this code, please cite the manuscript above (full reference to be added on publication).

## License

MIT License, see `LICENSE`.
