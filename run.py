from argparse import ArgumentParser
from time import time
from solair.simulation import Simulator
from solair.constants import constants
from solair.design import Tube

parser = ArgumentParser()
parser.add_argument("-n", "--n_segments", type=int, default=100)
parser.add_argument("-v", "--verbose", type=int, default=1)
parser.add_argument("-d", "--max_depth", type=int, default=100)


def main(args):
    start = time()
    sim = Simulator(
        Tube(),
        verbose=args.verbose,
        n_segments=args.n_segments,
        max_iterations=args.max_depth,
        fast=False,
    )
    sim.run()
    end = time()
    print("Time: {:.2f} seconds".format(end - start))


if __name__ == "__main__":
    print(f"Mass flow rate of air: {constants.m_air} kg/s")
    args = parser.parse_args()
    main(args)
