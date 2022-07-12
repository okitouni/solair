# Great job guys this looks good :)
from collections import defaultdict
import numpy as np
from typing import Tuple
from constants import constants
from utils import get_enthalpy, drop_pressure


def lmtd(t_air_in: float, t_air_out: float, t_co2_in: float, t_co2_out: float) -> float:
    """Compute the LMTD in a segment.

    Args:
        t_air_in (float): Segment input temperature of air.
        t_air_out (float): Segment output temperature of air.
        t_co2_in (float): Segment input temperature of sCO2.
        t_co2_out (float): Segment output temperature of sCO2.

    Raises:
        ValueError: Need t_co2_out > t_co2_in >= t_air_out > t_air_in for the physics to make sense. 

    Returns:
        float: The LMTD of the segment.
    """
    delta = (t_co2_in - t_air_out) - (t_co2_out - t_air_in)
    try:
        assert t_co2_out > t_co2_in and t_air_out > t_air_in and t_co2_in > t_air_out
    except AssertionError:
        raise ValueError(
            "Need t_co2_out > t_co2_in > t_air_out > t_air_in\n"
            "have: t_co2_out = {}, t_co2_in = {}, t_air_out = {}, t_air_in = {}".format(
                t_co2_out, t_co2_in, t_air_out, t_air_in
            )
        )
    return delta / np.log((t_co2_in - t_air_out) / (t_co2_out - t_air_in))


def energy_co2(
    p_in: float, p_out: float, t_in: float, t_out: float, m_co2: float
) -> float:
    """Compute the energy transferred out of the sCO2.

    Args:
        p_in (float): Initial pressure of the sCO2.
        p_out (float): Final pressure of the sCO2.
        t_in (float): Initial temperature of the sCO2.
        t_out (float): Final temperature of the sCO2.

    Returns:
        float: The energy transferred out of the sCO2.
    """
    h_in = get_enthalpy(p_in, t_in)
    h_out = get_enthalpy(p_out, t_out)
    return m_co2 * abs(h_out - h_in)


def energy_air(mdot: float, t0: float, t1: float) -> float:
    """Compute the energy transferred into the air.

    Args:
        mdot (float): Mass flow rate of the air.
        t0 (float): Temperature of the inlet air.
        t1 (float): Temperature of the outlet air.

    Returns:
        float: The energy transferred into the air.
    """
    return constants.cp_air * mdot * (t1 - t0)


class Simulator:
    def __init__(self) -> None:
        self.converged = False
        self.temps = defaultdict(list)
        self._verbose: int

    def run(
        self,
        t_co2_in,
        t_air_in,
        p_co2_in=None,
        verbose=1,
        max_segments: int = 50,
        max_depth: int = 20,
    ) -> None:
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
        def binary_search_temps(left, right, max_depth):
            if max_depth == 0:
                raise RecursionError("Max depth reached")
            t_co2_out = (left + right) / 2
            p_co2_out = drop_pressure(p_co2_in)
            m_co2_per_tube = constants.m_co2 / constants.n_tubes_tot
            q_co2 = energy_co2(p_co2_in, p_co2_out, t_co2_in, t_co2_out, m_co2_per_tube)
            # TODO get m_air for segment
            m_air_segment = constants.m_air / constants.n_segments
            # compute air out temp from heat transfer into air
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
            q_htc = constants.OHTC * (delta_t_m)  # TODO OHTC is not constant
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

