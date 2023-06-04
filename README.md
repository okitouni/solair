# Solair
(WIP)
## Overview
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