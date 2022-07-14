from constants import constants
import numpy as np
from utils import drop_pressure, energy_co2
from design import compute_ohtc
from collections import defaultdict


### TODO stop passing everything in as arguments 
class Simulator:
    def __init__(self, tube):
        self.tube = tube
        self.temps = defaultdict(list)
        self.converged = False
        self.n_segments = 100

    def _solve_segment(
        self,
        p_co2_in: float,
        t_co2_in: float,
        t_air_in: float,
        upstream: bool = False,
        max_temp: float = constants.t_co2.intlet,
        min_temp: float = constants.t_co2.outlet,
    ):
        # initial guess for t_c02_out
        if upstream:
            # going upstream co2 out gets hotter
            left_bound = t_co2_in
            right_bound = max_temp
        else:
            # going downstream co2 out gets colder
            left_bound = min_temp
            right_bound = t_co2_in
        # binary search for t_c02_out
        t_co2_out = self.binary_search(p_co2_in, t_air_in, left_bound, right_bound)
        # solve for t_air_out (AGAIN)
        p_co2_out = drop_pressure(t_c02_in=t_co2_in)  # TODO make it give you delta P
        q_co2 = abs(
            energy_co2(
                p_co2_in, p_co2_out, t_co2_in, t_co2_out, constants.m_co2_per_tube
            )
        )
        t_air_out = q_co2 / (constants.m_air_segment * constants.cp_air) + t_air_in

    def _solve_tube(
        self,
        p_co2_init: float,
        t_co2_init: float,
        t_air_init: float,
        upstream: bool = False,
    ):
        t_co2_out, t_air_out = t_co2_init, t_air_init
        for _ in range(self.n_segments):
            t_c02_out, t_air_out = self._solve_segment(
                p_co2_init, t_co2_init, t_co2_out, t_air_init, upstream
            )
            self.temps["t_co2"].append(t_c02_out)
            self.temps["t_air"].append(t_air_out)
        return t_c02_out, t_air_out

    def binary_search(
        self,
        p_co2_in: float,
        t_air_in: float,
        left_bound: float,
        right_bound: float,
        max_depth: int = 20,
        upstream: bool = False,
    ) -> float:
        """Binary search for the temperature of the final CO2. Depends on direction of the flow.

        Args:
            t_air_in (float): temperature of the initial air. This is always in the direction of air flow.
            left_bound (float): left bound of the search.
            right_bound (float): right bound of the search.
            max_depth (int, optional): maximum number of iterations. Defaults to 20.
            upstream (bool, optional): True if solver is going upstream of CO2 flow, False if going downstream. Defaults to False.

        Returns:
            float: temperature of the final CO2.
        """
        midpoint = (left_bound + right_bound) / 2
        if max_depth == 0:
            return midpoint

        out = self.energy_balance(
            p_co2_in=p_co2_in, t_co2_in=t_air_in, t_co2_out=midpoint, t_air_in=t_air_in
        )
        if out == 1:
            return midpoint  # TODO also return the temperature of the Air.
        elif out == 2:  # too little heat from co2
            if upstream:
                # increase co2 out temp
                return self.binary_search(
                    midpoint, right_bound, max_depth - 1, upstream
                )
            else:
                # decrease co2 out temp
                return self.binary_search(left_bound, midpoint, max_depth - 1, upstream)
        elif out == 3:  # too much heat from co2
            if not upstream:
                # increase co2 out temp
                return self.binary_search(
                    midpoint, right_bound, max_depth - 1, upstream
                )
            else:
                # decrease co2 out temp
                return self.binary_search(left_bound, midpoint, max_depth - 1, upstream)

    def energy_balance(
        self, p_co2_in: float, t_co2_in: float, t_co2_out: float, t_air_in: float
    ) -> int:
        """Check if the energy balance is satisfied. Compute the heat transfer coefficient from the geometry of the tube.
        Compare against the heat lost/gained by the CO2.

        Args:
            p_co2_in (float): pressure of the initial co2 (depends on the direction of the flow).
            t_co2_in (float): temperature of the initial co2 (depends on the direction of the flow).
            t_co2_out (float): temperature of the final co2 (depends on the direction of the flow). This is what you guess.
            t_air_in (float): temperature of the initial air. This is always in the direction of air flow.

        Returns:
            int: 1 if the energy balance is satisfied, 2 if co2 emitted too little heat, 3 if co2 emitted too much heat.
        """
        p_co2_out = drop_pressure(t_c02_in=t_co2_in)  # TODO make it give you delta P
        q_co2 = abs(
            energy_co2(
                p_co2_in, p_co2_out, t_co2_in, t_co2_out, constants.m_co2_per_tube
            )
        )
        t_air_out = q_co2 / (constants.m_air_segment * constants.cp_air) + t_air_in

        # compute the ohtc
        t_air = (
            t_air_in + t_air_out
        ) / 2  # TODO: check if this is correct (LMTD instead?)
        p_co2 = (p_co2_in + p_co2_out) / 2
        t_co2 = (t_co2_in + t_co2_out) / 2
        p_air = constants.p_air_in  # TODO: compute air pressure drop
        q_htc = compute_ohtc(
            t_air,
            t_co2,
            p_air,
            p_co2,
            constants.m_air_segment,
            constants.m_co2_per_tube,
            self.tube,
        )
        if np.isclose(q_htc, q_co2, rtol=constants.tolerance):
            return 1
        else:
            if q_htc > q_co2:
                return 2
            else:
                return 3

    # def get_t_air_out(t_co2_in: float, t_air_in: float, p_co2_in: float, t_co2_out: float) -> float:
    #     p_co2_out = drop_pressure(t_c02_in=t_co2_in) # TODO make it give you delta P
    #     q_co2 = abs(energy_co2(p_co2_in, p_co2_out, t_co2_in, t_co2_out, constants.m_co2_per_tube))
    #     t_air_out = q_co2 / (constants.m_air_segment * constants.cp_air) + t_air_in
    #     return t_air_out
