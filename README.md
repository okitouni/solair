# Solair
(WIP)
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