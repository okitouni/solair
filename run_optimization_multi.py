from turbo import Turbo1, TurboM
import numpy as np
import torch
import time
import pandas as pd

from optimization import CSP

def run_optimization(turbo_m = True , 
                    log = True,
                    output_file_name = "output", 
                    n_evals = 1000, 
                    t_air_inlet = 20):

    """
    Run optimization for a CSP cooler with specified air inlet temperature.
    
    Parameters
    ----------
        turbo_m : bool, optional - Use TuRBO-M, by default True
        log : bool, optional - Use Create log file by default True
        output_file_name : str, optional - Name of the output file, by default "output"
        n_evals : int, optional - Number of evaluations, by default 1000
        t_air_inlet : float, optional - Inlet temperature of the air, by default 20
    """

    ############### MAIN OPTIMIZATION ################

    # start timer
    time_start = time.time()

    # select device
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # define objective function
    f = CSP(log=log, t_air_inlet=t_air_inlet)
    Turbo = TurboM if turbo_m else Turbo1
    kwargs = dict(
        f=f,  # Handle to objective function
        lb=f.lb,  # Numpy array specifying lower bounds
        ub=f.ub,  # Numpy array specifying upper bounds
        n_init=20,  # Number of initial bounds from an Latin hypercube design
        max_evals=n_evals,  # Maximum number of evaluations
        batch_size=5,  # How large batch size TuRBO uses
        verbose=True,  # Print information from each batch
        use_ard=True,  # Set to true if you want to use ARD for the GP kernel
        max_cholesky_size=2000,  # When we switch from Cholesky to Lanczos
        n_training_steps=50,  # Number of steps of ADAM to learn the hypers
        min_cuda=1024,  # Run on the CPU for small datasets
        device=device,  # "cpu" or "cuda"
        dtype="float32",  # float64 or float32 
    )
    if turbo_m:
        kwargs["n_trust_regions"] = 5 # Number of trust regions
    turbo1 = Turbo(**kwargs)

    # run optimization
    turbo1.optimize()
    # get best parameters
    best_index = turbo1.fX.argmin()
    print("Done")
    print("Best Parameters:", turbo1.X[best_index])
    print("Best Cost:", turbo1.fX[best_index])
    with open(output_file_name + ".log", "w") as f:
        date = time.strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"{date}\n")
        time_final = time.time()
        f.write("Total time: " + str(time_final - time_start) + "\n")
        f.write(str(turbo1.X[best_index]) + "\n")
        f.write(str(turbo1.fX[best_index]) + "\n")
    np.savez(f"{output_file_name}.npz", X=turbo1.X, fX=turbo1.fX)
    return turbo1

if __name__ == "__main__":

    turbo_m = False
    log = True
    n_evals = 100

    out_folder = "output/"

    # import data from csv file "tops_df.csv"
    tops_df = pd.read_csv("tops_df.csv")
    t_in_vec = tops_df["temp"].values

    for t_in in t_in_vec:
        output_file_name = f"{out_folder}output_{t_in}"
        turbo_out = run_optimization(turbo_m, log, output_file_name, n_evals, t_in)

        # add best parameters to tops_df
        best_index = turbo_out.fX.argmin()
        tops_df.loc[tops_df["temp"] == t_in, "best_params"] = str(turbo_out.X[best_index])
        tops_df.loc[tops_df["temp"] == t_in, "best_cost"] = turbo_out.fX[best_index]

    # save tops_df to csv file
    tops_df.to_csv(f"{out_folder}tops_df_optimized.csv", index=False)

