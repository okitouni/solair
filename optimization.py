# %%
from turbo import Turbo1, TurboM
import numpy as np
from solair.simulation import Simulator, Simulator_Ehasan, DynamicLength
from solair.design import Tube
from solair.cost import calculate_total_cost_air_cooler, calculate_sub_cost_air_cooler
from solair.constants import constants
import torch
import time

torch.manual_seed(0)
torch.cuda.manual_seed_all(0)
np.random.seed(0)


class Csp:
    def __init__(self, log=True):
        # tube_in_diameter, tube_outer_inner_diff, fin_in_diameter, fin_outer_inner_diff
        self.lb = np.array(
            [
                20e-3 - 10e-3,  # tube_in_diameter
                1.1,  # tube_out_diameter = multiple of tube_in_diameter
                1.1,  # fin_in_diameter = multiple of tube_out_diameter
                1.1,  # fin_out_diameter = multiple of fin_in_diameter
                1.1,  # tube transverse pitch = multiple of fin_out_diameter
                1e-3,  # fin_pitch
                0.1,  # fin_thickness = fraction of fin_pitch
            ]
        )
        self.ub = np.array(
            [
                20e-3 + 20e-3,  # tube_in_diameter
                2,  # tube_out_diameter = multiple of tube_in_diameter
                2,  # fin_in_diameter = multiple of tube_out_diameter
                2,  # fin_out_diameter = multiple of fin_in_diameter
                2,  # tube transverse pitch = multiple of fin_out_diameter
                4e-3,  # fin_pitch
                0.8,  # fin_thickness = fraction of fin_pitch
            ]
        )
        self.log = log
        if self.log:
            with open("output_steps.log", "w") as f:
                f.write(f"{'x':8s} {'tube_len [m]':>8s} {'costs':>8s}\n")

                 

    def __call__(self, x):
        assert len(x) == len(self.ub)
        assert x.ndim == 1
        assert np.all(x <= self.ub) and np.all(x >= self.lb)
        tube_in_diameter = x[0]
        tube_out_diameter = tube_in_diameter * x[1]
        fin_in_diameter = tube_out_diameter * x[2]
        fin_out_diameter = fin_in_diameter * x[3]
        tube_transverse_pitch = fin_out_diameter * x[4]
        fin_pitch = x[5]
        fin_thickness = x[6] * fin_pitch
        lifetime_years = 25
        LCOE_fanpower_cents = 0.05

        tube = Tube(
            tube_in_diameter=tube_in_diameter,
            tube_out_diameter=tube_out_diameter,
            fin_in_diameter=fin_in_diameter,
            fin_out_diameter=fin_out_diameter,
            fin_pitch=fin_pitch,
            fin_thickness=fin_thickness,
            tube_transverse_pitch=tube_transverse_pitch,
        )
        sim = DynamicLength(tube, verbose=0, n_rows=4, n_sub_shx=1, fast=False,)
        sim.run()
        tube.n_segments = sim.n_segments
        # value = sim.results["t_co2"][-1][-1] # minimize the last temperature of the last tube

        costs = calculate_sub_cost_air_cooler(
            constants.rho_steel,
            constants.cost_steel,
            constants.rho_alu,
            constants.cost_alu,
            lifetime_years,
            LCOE_fanpower_cents,
            tube,
        )
        cost = calculate_total_cost_air_cooler(*costs)
        if self.log:
            with open("output_steps.log", "a") as f:
                tube_len = tube.segment_length * tube.n_segments
                f.write(f"{x} {tube_len} {costs}\n")
        return cost


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-n", "--n_evals", help="number of evaluations for BO", type=int, default=100)
    parser.add_argument("-o", "--output", help="output file name for run results and logs", type=str, default="output")
    parser.add_argument("-m", "--turbo-m", help="use turboM instead of tubro1", action="store_true")
    parser.add_argument("-l", "--no-log", help="no log in every call", action="store_true")
    args = parser.parse_args()
    time_start = time.time()

    log = not args.no_log
    f = Csp(log=log)
    Turbo = TurboM if args.turbo_m else Turbo1
    kwargs = dict(
        f=f,  # Handle to objective function
        lb=f.lb,  # Numpy array specifying lower bounds
        ub=f.ub,  # Numpy array specifying upper bounds
        n_init=20,  # Number of initial bounds from an Latin hypercube design
        max_evals=args.n_evals,  # Maximum number of evaluations
        batch_size=10,  # How large batch size TuRBO uses
        verbose=True,  # Print information from each batch
        use_ard=True,  # Set to true if you want to use ARD for the GP kernel
        max_cholesky_size=2000,  # When we switch from Cholesky to Lanczos
        n_training_steps=50,  # Number of steps of ADAM to learn the hypers
        min_cuda=1024,  # Run on the CPU for small datasets
        device="cpu",  # "cpu" or "cuda"
        dtype="float32",  # float64 or float32 # Number of trust regions
    )
    if args.turbo_m:
        kwargs["n_trust_regions"] = 10
    turbo1 = Turbo(**kwargs)

    turbo1.optimize()
    best_index = turbo1.fX.argmin()
    print("Done")
    print("Best Parameters:", turbo1.X[best_index])
    print("Best Cost:", turbo1.fX[best_index])
    with open(args.output + ".log", "w") as f:
        date = time.strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"{date}\n")
        f.write("Args: " + str(args) + "\n")
        time_final = time.time()
        f.write("Total time: " + str(time_final - time_start) + "\n")
        f.write(str(turbo1.X[best_index]) + "\n")
        f.write(str(turbo1.fX[best_index]) + "\n")
    np.savez(f"{args.output}.npz", X=turbo1.X, fX=turbo1.fX)

# %%
