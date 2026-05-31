# Autonomous Vertical Landing of Reusable Rockets: An Integrated GNC Architecture with Model Predictive Control and a First-Principles 6-DOF Nonlinear Plant

Companion code for the undergraduate thesis by Matheus Martins Marinato, UFSM, 2026, published at https://repositorio.ufsm.br/handle/1/25118.

MATLAB and Simulink implementation of an integrated guidance, navigation,
and control (GNC) architecture for the terminal vertical landing phase
of a Falcon-9-class reusable first stage, evaluated in a six-degree-of-freedom
nonlinear plant model derived from first principles.

## Repository contents

All files live in a single flat directory. Keep them together for the
scripts to run.

### Simulation models

- `Integrated_GNC_and_Rocket_Model_Sim.slx` — top-level closed-loop
  Simulink model integrating plant, navigation, guidance, and control.
- `rocket_6dof_sim.slx` — old version, no navigation pipeline, used
  for comparison.
- `planta_foguete.slx` — plant subsystem with Thurst=0, used for the Drop Test performed by the Guidance.
- `rocket_6dof.m` — nonlinear 6-DOF dynamics function shared between
  plant and EKF prediction step.

### Plant submodels

- `aero_model.m` — aerodynamic drag force.
- `gravity_body.m` — gravity in the body frame.
- `tvc_model.m` — thrust-vector-control force and moment.
- `cg_model.m` — center-of-gravity migration with propellant burn.
- `inertia_model.m` — inertia-tensor scheduling.
- `kinematics_euler.m` — Euler-angle kinematics.
- `rotation_matrix.m` — body-to-inertial rotation.

### Guidance

- `run_scurve.m` — quintic minimum-jerk trajectory generator with the
  Predictive True-State Scan ignition algorithm.
- `gerar_traj_rapido.m` — fast trajectory regenerator used inside the
  tuning pipeline.

### Navigation

- `init_ekf.m` — initialization of the bias-augmented Extended Kalman
  Filter and the tracker filter (covariances, gains, time constants).

### Control

- `init_mpc.m` — initialization of the LTV-MPC controller (weights,
  horizons, box constraints, QP setup).

### Tuning pipeline

- `otimizar_ganhos_bayesopt.m` — two-stage tuning pipeline combining
  Latin-Hypercube exploration and Bayesian optimization.

### Validation campaign

- `envelope_mais_monte_carlo_final.m` — operational-envelope sweep and
  dense Monte Carlo campaign with stochastic wind disturbances.
- `metrica_qualidade_pouso.m` — touchdown-quality metrics for a single
  simulation.
- `simulate_openloop.m` — open-loop simulation utility.

### Vehicle and simulation parameters

- `falcon9_params.m` — vehicle parameters and **initial conditions**
  (this is where ICs are set for single runs).
- `init_sim.m` — simulator initialization (solver settings, sample times).

### Verification suite (6-DOF model)

- `test_1dof.m` — one-degree-of-freedom verification battery.
- `test_3dof.m` — three-degree-of-freedom verification battery.
- `test_6dof.m` — six-degree-of-freedom verification battery.

Each test script contains usage instructions in its header.

### Plotting

- `plot_6dof.m` — plotting utilities for closed-loop responses.
- `plotar_graficos_pro_relatoriofinal.m` — figure generation for the
  thesis report.

### Reference data

- `mc_denso.mat` — dense Monte Carlo results (1000 runs).
- `crashes_dense.csv` — crash diagnostics from the dense campaign.

## Requirements

- MATLAB R2023b or later.
- Simulink.
- Optimization Toolbox (for `quadprog`).
- Statistics and Machine Learning Toolbox (for `bayesopt`).

## How to run

### Single closed-loop simulation

1. Open `falcon9_params.m` and set the **initial conditions** for the
   run.
2. From the MATLAB command window, in this order:
```matlab
   init_sim
   init_mpc
   init_ekf
   run_scurve
```
3. Open and run the Simulink model `Integrated_GNC_and_Rocket_Model_Sim.slx`.
4. After the run, evaluate the touchdown:
```matlab
   metrica_qualidade_pouso
```

### Operational-envelope sweep and dense Monte Carlo campaign

```matlab
init_ekf
envelope_mais_monte_carlo_final
```

The campaign takes approximately eight hours on a desktop for 1000 runs. Results are saved to `mc_denso.mat` and
`mc_envelope.mat`, with crash diagnostics in `crashes_dense.csv` and
`crashes_envelope.csv`.

### Plant-model verification tests

Run the scripts `test_1dof.m`, `test_3dof.m`, and `test_6dof.m` from
the command window. Each script's header explains the specific test
battery it executes and the expected outputs.

## Notes

- All files must reside in the same working directory for the scripts
  to find each other.
- All scripts and parameters use SI units.
- The reference vehicle is the SpaceX Falcon 9 Block 5; parameters are
  taken from publicly available sources documented in the thesis.

## Citation

If you use this code, please cite:

> Marinato, M. M. *Autonomous Vertical Landing of Reusable Rockets:
> An Integrated GNC Architecture with Model Predictive Control and a
> First-Principles 6-DOF Nonlinear Plant.* Undergraduate Thesis,
> Federal University of Santa Maria, 2026.
> Available at the website: https://repositorio.ufsm.br/handle/1/25118

The archived snapshot corresponding to the thesis defense is available
on Zenodo with a permanent DOI (see release `v1.0-thesis`).

## License

MIT License. See `LICENSE` for details.
