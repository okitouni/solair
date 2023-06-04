<!-- center the title below make text large and colored -->
<p align="center">
    <font size="6">Solair: Optimization of Air-Cooling for sCO2</font>
</p>

# Overview
The optimization process in the provided code is based on the TurboM and Turbo1 algorithms, which are Bayesian optimization methods. The main goal of the optimization is to find the optimal parameters for a heat exchanger design that minimizes the total cost while maintaining certain performance requirements.

The `optimization.py` file contains three main classes:

`CSP`: This class represents the objective function for the optimization problem. It takes an input vector x representing various design parameters and returns the total cost associated with those parameters.

`Simulator`: This class, defined in simulation.py, simulates a heat exchanger with given design parameters (tube properties, number of sub-heat exchangers, etc.). It calculates temperatures and pressures of air and CO2 at different points in the system.

`DynamicLength`: This class inherits from Simulator and extends its functionality to handle dynamic segment lengths in a tube.

The optimization process starts by initializing an instance of the CSP class, which defines bounds for each parameter in vector x. Then, depending on whether TurboM or Turbo1 is selected (based on command-line arguments), an instance of either TurboM or Turbo1 is created with appropriate settings such as lower bounds (lb), upper bounds (ub), maximum evaluations (max_evals), batch size, etc.

The optimization algorithm then iteratively evaluates candidate solutions (sets of design parameters) using Bayesian optimization techniques to find a solution that minimizes the objective function (total cost). At each iteration, it updates its internal model based on new evaluations and uses this information to guide future exploration/exploitation steps.

After completing all iterations, it reports the best-found solution along with its corresponding cost value. The results are saved into `.npz` files for further analysis.

An example of such analysis is presented in `analysis.ipynb`. This notebook loads the results from the `.npz` files and plots the cost values over iterations. It also plots the best-found solution and compares it to the baseline solution which was extracted from Ehsan et al.
## Requirements
You will probably need to run `pip install CoolProp` to install the right version of CoolProp. Otherwise, you'll also need some standard packages you can install with conda (or whatever else): gpytorch, matplotlib, etc.

## How to use
The high-level scripts you should be using are `run.py` and `optimization.py`.

Here's a breakdown of what everything does (roughly).
- `run.py`: Run a single simulation based on the CLI arguments in the file.
- `optimization.py`: Essentially runs `run.py` multiple times inside a TurBO optimization loop. Output is written into the directory `outputs`.
Here's an example of how to run the optimization script:
```
>>> python optimization -n 2 -o output
Using dtype = torch.float32 
Using device = cpu
Starting from fbest = 1.126e+06
4) New best: 9.57e+05
Done
Best Parameters: [1.35738670e-02 1.77686832e+00 1.36968161e+00 1.26345163e+00
1.33573795e+00 3.87276112e-03 7.82251524e-01 2.81413635e+01]
Best Cost: [957018.29638951]
Saved to: outputs/output.npz  outputs/output.log outputs/output_steps.log
```
You can run `python optimization --help` to get more details on the CLI arguments.
- `constants.py`: 
- `design.py`: 
- `simulator.py`: 
(Some of the simulation configuration is in `constants.py`) 



# Simulator Implementation
The simulator is designed to model and analyze the performance of a heat exchanger system, specifically focusing on air-cooled CO2 coolers. The implementation takes into account various design parameters, such as tube dimensions, fin properties, and flow conditions. It employs a combination of physics-based equations and empirical correlations to simulate the heat transfer process between CO2 and air streams.

## Assumptions
The following assumptions are made in the simulator:

Steady-state operation: The simulation assumes that the heat exchanger operates under steady-state conditions, with constant mass flow rates for both CO2 and air.

Constant properties: Thermophysical properties of fluids (CO2 and air) are assumed to be constant within each segment.

Uniform air distribution: Airflow across the tubes is assumed to be uniform.

Negligible pressure drop for air: The pressure drop in the air side is not considered in this simulation.

Segmented approach: The heat exchanger tubes are divided into multiple segments, with each segment being treated as an individual heat exchanger unit.

## Physics Equations
The simulator uses several physics-based equations to model different aspects of the heat exchanger:

Energy balance equation: Ensures that the energy transferred from CO2 to air through convection equals the energy change in CO2 due to temperature variation and pressure drop.

$$Q_{\rm co2} = Q_{\rm htc} $$

where $Q_{\rm co2}$ is the energy change in CO2, and $Q_{\rm htc}$ is the energy transferred through convection.

Log Mean Temperature Difference (LMTD): Calculates an average temperature difference between hot (CO2) and cold (air) streams across a segment.

$$\Delta T_m = {\rm lmtd}(T_{\rm air_{in}}, T_{\rm air_{out}}, T_{\rm co2_{in}}, T_{\rm co2_{out}})$$

Overall Heat Transfer Coefficient (OHTC): Computes an effective heat transfer coefficient based on tube geometry, fluid properties, and flow conditions.

   $$ {O_{\rm HTC}}_{(i,j)} = {\rm computeOHTC}(T_{\rm air_{in}}, T_{\rm co2_{out}}, P_{\rm air_{in}}, P_{\rm co2_{out}}, \dot{m}_{\rm air}, \dot{m}_{\rm co2}, {\rm tube state}) $$

Pressure drop calculation for CO2: Estimates pressure drop across a segment using Darcy-Weisbach equation or other appropriate correlations.

$$\Delta_p = {\rm drop_pressure}(P_{\rm co2_{in}}, T_{\rm co2_{in}}, \dot{m}_{\rm co2}, {\rm tube state})$$

## Solution Methodology
The simulator employs an iterative approach to solve for temperatures and pressures at different points within the system:

Initialization: Create a Simulator or DynamicLength instance with specified design parameters (tube object), verbosity level (verbose), maximum iterations (max_iterations), number of sub-heat exchangers (n_sub_shx), etc.

Tube-level solution: For each tube row in a sub-heat exchanger (SHX), call _solve_tube() method with initial conditions for CO2 pressure (p_co2_init), temperature (t_co2_init), and optionally air inlet temperature(s) (t_air_init). This method solves for temperatures and pressures at each segment along a tube by calling _solve_segment() method iteratively.

Segment-level solution: In _solve_segment(), it performs binary search on output CO₂ temperature until energy balance equation is satisfied within specified tolerance limits.

Sub-Heat Exchangers (SHXs): For intermediate SHXs after the first one (_intermediate_shx()), perform binary search on initial conditions until mean outlet temperature converges within tolerance limits.

Results collection: Store calculated temperatures and pressures for each segment in results dictionary during simulation process.

Execution: Call run() method on simulator instance to execute entire simulation process across all SHXs sequentially.

By combining these physics-based equations with iterative solution techniques like binary search and segmented approach for solving complex systems like multi-row heat exchangers with varying lengths or designs efficiently while maintaining accuracy in performance calculations.