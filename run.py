from argparse import ArgumentParser
from lmtd_simulator import Simulator
from constants import constants


parser = ArgumentParser()
parser.add_argument("-n", "--n_segments", type=int, default=100)
parser.add_argument("-v", "--verbose", type=int, default=2)
parser.add_argument("-d", "--max_depth", type=int, default=20)



def main():
    args = parser.parse_args()
    simulator = Simulator()
    simulator.run(constants.t_co2_outlet, constants.t_air_inlet, max_segments=args.n_segments, verbose=args.verbose)
    return 1

if __name__ == "__main__":
    main()