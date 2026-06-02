# SFRIO-HyDCMG MATLAB Example

This repository provides a minimal MATLAB implementation of the
statistical-feasibility-based robust interval optimization (SFRIO) model
for hybrid renewable-energy data-center microgrid scheduling.

The code is prepared as a compact reproducibility package for the manuscript.
It includes the main optimization script, the uncertainty-set construction
routine, processed input profiles, and result-saving utilities.

## Requirements

- MATLAB
- YALMIP
- Gurobi with a valid MATLAB interface

The optimization model uses YALMIP functions such as `sdpvar`, `binvar`,
`sdpsettings`, `kkt`, `optimize`, and `value`.

## Quick Start

From MATLAB, run:

```matlab
cd path/to/SFRIO-HyDCMG-MATLAB
run_example
```

By default, `main_intervalrobust.m` runs one representative case with
`epsilon = 0.4` and `delta = 0.4`. To run other settings, edit:

```matlab
epsilon_list = [0.4];
delta_list = [0.4];
```

near the top of `main_intervalrobust.m`.

The empirical quantile constraints are calibrated directly with `epsilon`;
no additional relaxation or post-hoc adjustment is applied.

In this implementation, the chance constraints specify the target
probabilistic requirements for RE availability and data-center temperature.
They are not solved as a generic chance-constrained program. Instead, the
code uses empirical quantiles to build statistically feasible robust
surrogate constraints, following the statistical-feasibility-based robust
optimization formulation in the manuscript.

## Outputs

Results are written to the `results/` folder. The output file name follows:

```text
MultiDCMG3_DataInter_NoGuess_e*_d*.mat
```

Important saved variables include:

- `cost_mat`: total operation cost along the generated Pareto front.
- `wind_mat`: renewable-energy curtailment index along the Pareto front.
- `constraint_over_prob_list`: empirical violation probability summary.
- `temperature_over_prob_list`: empirical temperature-constraint violation probability.
- `re_over_prob_list`: empirical renewable-energy violation probability.
- `temperature_zone_over_prob_list`: zone-wise diagnostic temperature violation probability.
- `H_DC_solution_mat`, `P_DC_solution_mat`, `P_cool_solution_mat`, and
  `P_re_solution_mat`: selected decision-variable trajectories.
- `P_zone_solution_tensor` and `P_cool_zone_solution_tensor`: zone-wise
  diagnostic power and cooling allocations.

The legacy variable `prob_mat` is intentionally not saved.

## Repository Structure

```text
.
|-- main_intervalrobust.m              # Main SFRIO scheduling script
|-- run_example.m                      # Minimal entry point
|-- valuesave.m                        # Result evaluation and saving helper
|-- paradef.m                          # System parameters
|-- vardef.m                           # YALMIP decision variables
|-- curves_construction.m              # Delay-tolerant workload curves
|-- Construcrt_Box_uncertainty_set.m   # Data-driven uncertainty set construction
|-- Seperate_dataset.m                 # Sample split helper
|-- uncertainty_set_construct_wind_power.m
|-- data/                              # Processed input data
`-- results/                           # Generated results
```

## Notes

The data files in `data/` are processed numerical profiles used by the
example case. They are provided to allow direct execution of the optimization
script. The raw public data sources and preprocessing assumptions should be
cited according to the manuscript.

## License

This code is released under the MIT License.
