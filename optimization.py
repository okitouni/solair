# %%
from turbo import Turbo1, TurboM
import numpy as np
from solair.simulation import DynamicLength
from solair.design import Tube
from solair.cost import (
    calculate_total_cost_air_cooler,
    calculate_sub_cost_air_cooler,
)
from solair import constants
import torch
import time
import os

torch.manual_seed(0)
torch.cuda.manual_seed_all(0)
np.random.seed(0)


class CSP:
    def __init__(self, logfile="", t_air_inlet=20):
        self.lb = np.array(
            [
                20e-3 - 10e-3,  # tube_in_diameter
                1.1,  # tube_out_diameter = multiple of tube_in_diameter
                1.1,  # fin_in_diameter = multiple of tube_out_diameter
                1.1,  # fin_out_diameter = multiple of fin_in_diameter
                1.01,  # tube transverse pitch = multiple of fin_out_diameter
                1e-3,  # fin_pitch
                0.1,  # fin_thickness = fraction of fin_pitch,
                t_air_inlet,  # t_air_out
            ]
        )
        self.ub = np.array(
            [
                20e-3 + 20e-3,  # tube_in_diameter
                2,  # tube_out_diameter = multiple of tube_in_diameter
                2,  # fin_in_diameter = multiple of tube_out_diameter
                2.1,  # fin_out_diameter = multiple of fin_in_diameter
                2.1,  # tube transverse pitch = multiple of fin_out_diameter
                4e-3,  # fin_pitch
                0.8,  # fin_thickness = fraction of fin_pitch
                30,  # t_air_out
            ]
        )
        self.logfile = logfile
        if logfile:
            with open(self.logfile, "w") as f:
                f.write(
                    f"{'x array':8s} {'tube_len [m]':>8s} {'costs array':>8s}\n"
                )

        # t_air_in as a variable attribute
        self.t_air_inlet = t_air_inlet

    def _run_simulation(self, x):
        # error handling
        if len(x) != len(self.lb):
            raise ValueError(
                f"x has length {len(x)} but should have length {len(self.lb)}"
            )
        assert x.ndim == 1, "x should be a 1D array"
        for i in range(len(x)):
            if x[i] < self.lb[i] or x[i] > self.ub[i]:
                raise ValueError(
                    f"x[{i}] = {x[i]} is out of bounds, [{self.lb[i]}, {self.ub[i]}]"
                )

        tube_in_diameter = x[0]
        tube_out_diameter = tube_in_diameter * x[1]
        fin_in_diameter = tube_out_diameter * x[2]
        fin_out_diameter = fin_in_diameter * x[3]
        tube_transverse_pitch = fin_out_diameter * x[4]
        fin_pitch = x[5]
        fin_thickness = x[6] * fin_pitch
        t_air_out = x[7]

        constants_t = constants(self.t_air_inlet, t_air_out)
        tube = Tube(
            tube_in_diameter=tube_in_diameter,
            tube_out_diameter=tube_out_diameter,
            fin_in_diameter=fin_in_diameter,
            fin_out_diameter=fin_out_diameter,
            fin_pitch=fin_pitch,
            fin_thickness=fin_thickness,
            tube_transverse_pitch=tube_transverse_pitch,
            constants_t=constants_t,
        )
        sim = DynamicLength(tube, verbose=0, n_sub_shx=1)
        sim.run()
        tube.n_segments = sim.n_segments

        tube_cost, fan_cost = calculate_sub_cost_air_cooler(
            constants_t.rho_steel,
            constants_t.cost_steel,
            constants_t.rho_alu,
            constants_t.cost_alu,
            constants_t.lifetime_years,
            constants_t.LCOE_fanpower_cents,
            tube,
        )
        cost = calculate_total_cost_air_cooler(tube_cost, fan_cost, tube)
        return cost, tube_cost, fan_cost, tube

    def __call__(self, x):
        cost, tube_cost, fan_cost, tube = self._run_simulation(x)
        if self.logfile:
            with open(self.logfile, "a") as f:
                tube_len = tube.segment_length * tube.n_segments
                f.write(
                    f"{x} {tube_len}, costs: tube {tube_cost}, fan {fan_cost}, total {cost} \n"
                )
        return cost


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "-n",
        "--n_evals",
        help="number of evaluations for BO",
        type=int,
        default=100,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output file name for run results and logs",
        type=str,
        default="output",
    )
    parser.add_argument(
        "-m",
        "--turbo-m",
        help="use turboM instead of tubro1",
        action="store_true",
    )
    parser.add_argument(
        "-s",
        "--silent",
        help="No log in every call. Will still log results.",
        action="store_true",
    )
    args = parser.parse_args()
    time_start = time.time()

    # select cuda device if available
    device = "cuda" if torch.cuda.is_available() else "cpu"

    os.makedirs("outputs", exist_ok=True)

    filename = os.path.join("outputs", args.output)

    log_steps_file = f"{filename}_steps.log" if not args.silent else ""

    f = CSP(logfile=log_steps_file)
    Turbo = TurboM if args.turbo_m else Turbo1
    kwargs = dict(
        f=f,  # Handle to objective function
        lb=f.lb,  # Numpy array specifying lower bounds
        ub=f.ub,  # Numpy array specifying upper bounds
        n_init=min(
            args.n_evals, 20
        ),  # Number of initial bounds from an Latin hypercube design
        max_evals=args.n_evals + 1,  # Maximum number of evaluations
        batch_size=min(args.n_evals, 10),  # How large batch size TuRBO uses
        verbose=True,  # Print information from each batch
        use_ard=True,  # Set to true if you want to use ARD for the GP kernel
        max_cholesky_size=2000,  # When we switch from Cholesky to Lanczos
        n_training_steps=50,  # Number of steps of ADAM to learn the hypers
        min_cuda=1024,  # Run on the CPU for small datasets
        device=device,  # "cpu" or "cuda"
        dtype="float32",  # float64 or float32
    )
    if args.turbo_m:
        kwargs["n_trust_regions"] = 5  # Number of trust regions
    turbo1 = Turbo(**kwargs)

    turbo1.optimize()
    best_index = turbo1.fX.argmin()
    print("Done")
    print("Best Parameters:", turbo1.X[best_index])
    print("Best Cost:", turbo1.fX[best_index])
    # check if directories exist
    with open(f"{filename}.log", "w") as f:
        date = time.strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"{date}\n")
        f.write("Args: " + str(args) + "\n")
        time_final = time.time()
        f.write("Total time: " + str(time_final - time_start) + "\n")
        f.write(str(turbo1.X[best_index]) + "\n")
        f.write(str(turbo1.fX[best_index]) + "\n")
    np.savez(f"{filename}.npz", X=turbo1.X, fX=turbo1.fX)

    final_msg = f"Saved to: {filename}.npz  {filename}.log"
    if not args.silent:
        final_msg += f" {log_steps_file}"
    print(final_msg)
