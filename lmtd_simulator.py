from collections import defaultdict
import numpy as np
from typing import Tuple
from constants import constants
from utils import drop_pressure, lmtd, energy_co2
from design import compute_ohtc, Tube



class Simulator:
    def __init__(self, tube: Tube = None) -> None:
        """Initialize the simulator with a specific tube design.

        Args:
            tube (Tube, optional): Tube object associated with the specific design. Defaults to None.
        """        
        self.converged = False # Whether the simulation has converged
        self.temps = defaultdict(list) # Dictionary of temperatures for each segment
        self._verbose: int # Verbosity level
        self.tube = Tube() if tube is None else tube 

    def run(
        self,
        t_co2_in: float,
        t_air_in: float,
        p_co2_in: float =None,
        verbose: int =1,
        max_segments: int = 50,
        max_depth: int = 20,
    ) -> None: 
        """Run the simulation. The initial values here are actually outlet values for CO2 and Inlet for air.

        Args:
            t_co2_in (float): Initial temperature of the sCO2.
            t_air_in (float): Initial temperature of the air.
            p_co2_in (float, optional): Initial pressure of the CO2. Defaults to None.
            verbose (int, optional): Verbosity level. Defaults to 1.
            max_segments (int, optional): Number of segments along the tubes. Defaults to 50.
            max_depth (int, optional): Maximum number of iterations in the binary search at each segment. Defaults to 20.
        """    
        self._verbose = verbose
        if verbose > 0:
            print("Initial conditions:")
            print("t_co2_in:", t_co2_in, "t_air_in:", t_air_in)

        p_co2_in = constants.p_co2_outlet if p_co2_in is None else p_co2_in

        for i in range(max_segments):
            t_co2_in, t_air_out = self._temperature_search(
                p_co2_in=p_co2_in,
                t_co2_in=t_co2_in,
                t_air_in=t_air_in,
                max_depth=max_depth,
            )
            self.temps["t_co2"].append(t_co2_in)
            self.temps["t_air"].append(t_air_out)

            if verbose > 1:
                print(
                    f"Finished segment {i+1}",
                    f"t_co2_out: {t_co2_in}, t_air_out: {t_air_out}",
                )
            if self.converged:
                break

        if verbose > 0:
            print("Finished simulation")

        # testing using final temperature for second tube
        t_co2_in, t_air_out = self._temperature_search(
            p_co2_in=p_co2_in,
            t_co2_in=t_co2_in,
            t_air_in=self.temps["t_air"][0],
            max_depth=max_depth,
        )
        print(f"t_co2_out: {t_co2_in}, t_air_out: {t_air_out}")
        return

    def _temperature_search(
        self, p_co2_in: float, t_co2_in: float, t_air_in: float, max_depth: int = 20,
    ) -> Tuple[float]:
        """Search for the temperature of the sCO2 in a given segment.

        Args:
            p_co2_in (float): Initial pressure of the sCO2. 
            t_co2_in (float): Initial temperature of the sCO2.
            t_air_in (float): Initial temperature of the air.
            max_depth (int, optional): Max iteration depth in the search. Defaults to 20.

        Returns:
            Tuple[float]: The final temperature of the sCO2 and air across the segment.
        """
        # Currently, this is a binary search and only works when data is sorted
        # (i.e. q_htc - q_co2 is monotonic in sCO2 output temp).
        def binary_search_temps(left: float, right: float, max_depth: int) -> Tuple[float]:
            """Search for the temperature of the sCO2 in a given segment. This is a binary search which will
            terminate when the temperature is within a certain tolerance of the target temperature. It assumes
            that the difference in Q_htc and Q_co2 is monotonic in the guess temperature.

            Args:
                left (_type_): Left bound of the search.
                right (_type_): Right bound of the search.
                max_depth (_type_): Maximum number of iterations in the search.

            Raises:
                RecursionError: If the maximum number of iterations is reached.

            Returns:
                Tuple: The final temperature of the sCO2 and air across the segment.
            """            
            if max_depth == 0:
                raise RecursionError("Max depth reached")
            t_co2_out = (left + right) / 2
            p_co2_out = drop_pressure(p_co2_in)
            m_co2_per_tube = constants.m_co2 / constants.n_tubes_tot
            q_co2 = energy_co2(p_co2_in, p_co2_out, t_co2_in, t_co2_out, m_co2_per_tube)
            # TODO get m_air for segment
            m_air_segment = constants.m_air / constants.n_segments / constants.n_tubes_in_row # divided by number of tubes in a row as we are considering a single tube for all calculations            # compute air out temp from heat transfer into air
            t_air_out = q_co2 / (m_air_segment * constants.cp_air) + t_air_in

            if np.isclose(t_co2_out, constants.t_co2_inlet, rtol=0.001):
                self.converged = True
                if self._verbose > 0:
                    print("Converged")
                return t_co2_out, t_air_out

            if t_air_out > t_co2_in:
                # air out can't be hotter than sCO2 used to heat it
                # air is too hot -> decrease out temp
                # -> set current middle to be the new max
                # if self._verbose > 2:
                #     print("t_co2_out:", t_co2_out, "t_air_out:", t_air_out),
                right = t_co2_out
                return binary_search_temps(left, right, 1000)

            delta_t_m = lmtd(t_air_in, t_air_out, t_co2_in, t_co2_out)

            # Compute the OHTC from the tube design and theormodynamics probabilities
            ohtc = compute_ohtc(
                t_air = t_air_in, t_s=t_co2_in, p_air=constants.p_air_in, p_s=p_co2_in, m_air=m_air_segment, m_s=m_co2_per_tube,  tube=self.tube
            )
            print(ohtc)

            q_htc = ohtc * (delta_t_m)  # TODO OHTC is not constant
            if self._verbose > 2:
                print(
                    f"t_co2_out: {t_co2_out}, t_air_out: {t_air_out}, q_co2: {q_co2}, q_htc: {q_htc}"
                )
            if not np.isclose(q_htc, abs(q_co2), rtol=constants.tolerance):
                if q_htc < q_co2:
                    # too much co2 heat change -> lower delta T -> lower out temp
                    # -> set current middle to be the new max
                    right = t_co2_out
                else:
                    left = t_co2_out
                return binary_search_temps(left, right, max_depth - 1)
            else:
                return t_co2_out, t_air_out

        # iterate until temperature converges.
        # max_temp_allowed = np.inf # constants.t_co2_inlet is what you want to end with
        t_co2_out, t_air_out = binary_search_temps(
            t_co2_in, constants.t_co2_inlet, max_depth
        )
        return t_co2_out, t_air_out

