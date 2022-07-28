from turbo import Turbo1
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from solair.simulation import Simulator, Simulator_Ehasan, DynamicLength
from solair.design import Tube
from solair.cost import calculate_total_cost_air_cooler
from solair.constants import constants
import torch

torch.manual_seed(0)
torch.cuda.manual_seed_all(0)
np.random.seed(0)

class Csp:
    def __init__(self):
        # tube_in_diameter, tube_outer_inner_diff, fin_in_diameter, fin_outer_inner_diff
        self.lb = np.array(
            [
                20e-3 - 10e-3,  # tube_in_diameter
                1.1, # tube_out_diameter = fraction of tube_in_diameter
                1.1,  # fin_in_diameter = fraction of tube_out_diameter
                1.1,  # fin_out_diameter = fraction of fin_in_diameter
            ]
        )
        self.ub = np.array(
            [
                20e-3 + 20e-3,  # tube_in_diameter
                2, # tube_out_diameter  - tube_in_diameter
                2, # fin_in_diameter = fraction of tube_out_diameter
                2, # fin_out_diameter = fraction of fin_in_diameter
            ]
        )

    def __call__(self, x):
        assert len(x) == len(self.ub)
        assert x.ndim == 1
        assert np.all(x <= self.ub) and np.all(x >= self.lb)
        tube_in_diameter = x[0]
        tube_out_diameter = tube_in_diameter * x[1]
        fin_in_diameter = tube_out_diameter  * x[2]
        fin_out_diameter = fin_in_diameter * x[3]

        tube = Tube(
            tube_in_diameter=tube_in_diameter,
            tube_out_diameter=tube_out_diameter,
            fin_in_diameter=fin_in_diameter,
            fin_out_diameter=fin_out_diameter,
        )
        sim = DynamicLength(tube, verbose=0, n_rows=4, n_sub_shx=1, fast=False,)
        sim.run()
        tube.n_segments = sim.n_segments
        # value = sim.results["t_co2"][-1][-1] # minimize the last temperature of the last tube
        cost = calculate_total_cost_air_cooler(
            constants.rho_steel,
            constants.cost_steel,
            constants.rho_alu,
            constants.cost_alu,
            tube,
        )

        return cost


if __name__ == "__main__":
    f = Csp()

    turbo1 = Turbo1(
        f=f,  # Handle to objective function
        lb=f.lb,  # Numpy array specifying lower bounds
        ub=f.ub,  # Numpy array specifying upper bounds
        n_init=20,  # Number of initial bounds from an Latin hypercube design
        max_evals=1000,  # Maximum number of evaluations
        batch_size=10,  # How large batch size TuRBO uses
        verbose=True,  # Print information from each batch
        use_ard=True,  # Set to true if you want to use ARD for the GP kernel
        max_cholesky_size=2000,  # When we switch from Cholesky to Lanczos
        n_training_steps=50,  # Number of steps of ADAM to learn the hypers
        min_cuda=1024,  # Run on the CPU for small datasets
        device="cpu",  # "cpu" or "cuda"
        dtype="float64",  # float64 or float32
    )

    turbo1.optimize()
    np.savez("fit_results.npz", X=turbo1.X, fX=turbo1.fX)