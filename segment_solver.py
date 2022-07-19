from constants import constants
import numpy as np
from utils import drop_pressure, energy_co2, lmtd
from design import compute_ohtc
from collections import defaultdict

### TODO stop passing everything in as arguments
class Simulator:
    def __init__(self, tube, verbose: int = 0):
        self.tube = tube
        self.properties = defaultdict(list)
        self.converged = False
        self.n_segments = 100
        self._verbose = verbose

    def _solve_segment(
        self,
        p_co2_in: float,
        t_co2_in: float,
        t_air_in: float,
        upstream: bool = False,
        max_temp: float = constants.t_co2_inlet,
        min_temp: float = constants.t_co2_outlet,
    ):
        """
        Solve temperature of output CO2 in a given segment. 
        Can go both upstream or downstream of sCO2 flow.

        Args:
            p_co2_in (float): 
                pressure of the initial CO2.
            t_co2_in (float): 
                temperature of the initial CO2.
            t_air_in (float): 
                temperature of the initial air. This is always in the direction of air flow.
            upstream (bool, optional): 
                Whether the solver runs upstream (opposite) of CO2 flow. Defaults to False.
            max_temp (float, optional): 
                Maximum allowed temperature. When going upstream this should be CO2 inlet temperature,
                otherwise, it should be current input temperature. Defaults to constants.t_co2.intlet.
            min_temp (float, optional): 
                Minimum allowed temperature. CO2 outlet tempreature when going downstream, 
                therwise, it should be the current input tempreature. Defaults to constants.t_co2.outlet.
        """
        # initial guess for t_c02_out
        if upstream:
            # going upstream co2 out gets hotter
            lower_bound = t_co2_in
            upper_bound = max_temp
        else:
            # going downstream co2 out gets colder
            lower_bound = min_temp
            upper_bound = t_co2_in
        # binary search for t_c02_out
        t_co2_out = self.binary_search(
            p_co2_in=p_co2_in,
            t_air_in=t_air_in,
            t_co2_in=t_co2_in,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            upstream=upstream,
        )
        # solve for t_air_out (AGAIN) # TODO remove this double calculations
        p_co2_out = drop_pressure(p_co2_in)  # TODO make it give you delta P
        q_co2 = abs(
            energy_co2(
                p_co2_in, p_co2_out, t_co2_in, t_co2_out, constants.m_co2_per_tube
            )
        )
        t_air_out = q_co2 / (constants.m_air_segment * constants.cp_air) + t_air_in
        return t_air_out, t_co2_out, p_co2_out

    def _solve_tube(
        self,
        p_co2_init: float,
        t_co2_init: float,
        t_air_init: float = None,
        upstream: bool = False,
    ):
        """
        Solve the temperature of an entire tube.

        Args:
            p_co2_init (float): 
                Pressure of the initial CO2.
            t_co2_init (float): 
                Temperature of the initial CO2.
            t_air_init (float): 
                Temperature of the initial air. This is always in the direction of air flow. If not given, 
                the solver uses the last tube air temperatures.
            upstream (bool, optional): 
                True if solver is going upstream of CO2 flow, False if going downstream. Defaults to False.

        Returns:
            _type_: _description_
        """
        tube_t_co2 = []
        tube_t_air = []
        tube_p_co2 = []

        t_co2_out = t_co2_init
        p_co2_out = p_co2_init

        if t_air_init is None:
            t_air_in_list = self.properties["t_air"][-1]
        else:
            t_air_in_list = [t_air_init] * self.n_segments
        if not upstream:
            t_air_in_list.reverse()

        for segment_num in range(self.n_segments):
            t_air_in = t_air_in_list[segment_num]

            t_air_out, t_co2_out, p_co2_out = self._solve_segment(
                p_co2_in=p_co2_out,
                t_co2_in=t_co2_out,
                t_air_in=t_air_in,
                upstream=upstream,
            )
            if self._verbose > 1:
                print(
                    f"Finished segment {segment_num+1}",
                    f"t_co2_out: {t_co2_out}, t_air_out: {t_air_out}",
                )
            if self.converged:
                break

            tube_t_co2.append(t_co2_out)
            tube_t_air.append(t_air_out)
            tube_p_co2.append(p_co2_out)

        self.properties["t_co2"].append(tube_t_co2)
        self.properties["t_air"].append(tube_t_air)
        self.properties["p_co2"].append(tube_p_co2)
        return tube_t_air, tube_t_co2, tube_p_co2

    def binary_search(
        self,
        p_co2_in: float,
        t_air_in: float,
        t_co2_in: float,
        lower_bound: float,
        upper_bound: float,
        max_depth: int = 200,
        upstream: bool = False,
    ) -> float:
        """
        Binary search for the temperature of the final CO2. Depends on direction of the flow.

        Args:
            t_air_in (float): 
                Temperature of the initial air. This is always in the direction of air flow.
            t_co2_in (float): 
                Temperature of the initial co2. Depends on the direction of the flow.
            lower_bound (float): 
                Left bound of the search.
            upper_bound (float): 
                Right bound of the search.
            max_depth (int, optional): 
                Maximum number of iterations. Defaults to 20.
            upstream (bool, optional): True 
                If solver is going upstream of CO2 flow, False if going downstream. Defaults to False.

        Returns:
            float: temperature of the final CO2.
        """
        midpoint = (lower_bound + upper_bound) / 2
        if max_depth == 0:
            return midpoint
        if upstream:
            out = self.energy_balance(
                p_co2_in=p_co2_in,
                t_co2_in=t_co2_in,
                t_co2_out=midpoint,
                t_air_in=t_air_in,
            )
        else:
            out = self.energy_balance(
                p_co2_in=p_co2_in,
                t_co2_in=midpoint,
                t_co2_out=t_co2_in,
                t_air_in=t_air_in,
            )

        if out == 1:
            return midpoint  # TODO also return the temperature of the Air.
        elif out == 2:  # too little heat from co2
            if upstream:
                # increase co2 out temp
                return self.binary_search(
                    p_co2_in=p_co2_in,
                    t_air_in=t_air_in,
                    t_co2_in=t_co2_in,
                    lower_bound=midpoint,
                    upper_bound=upper_bound,
                    max_depth=max_depth - 1,
                    upstream=upstream,
                )
            else:
                # decrease co2 out temp
                return self.binary_search(
                    p_co2_in=p_co2_in,
                    t_air_in=t_air_in,
                    t_co2_in=t_co2_in,
                    lower_bound=lower_bound,
                    upper_bound=midpoint,
                    max_depth=max_depth - 1,
                    upstream=upstream,
                )
        elif out == 3:  # too much heat from co2
            if not upstream:
                # increase co2 out temp
                return self.binary_search(
                    p_co2_in=p_co2_in,
                    t_air_in=t_air_in,
                    t_co2_in=t_co2_in,
                    lower_bound=midpoint,
                    upper_bound=upper_bound,
                    max_depth=max_depth - 1,
                    upstream=upstream,
                )
            else:
                # decrease co2 out temp
                return self.binary_search(
                    p_co2_in=p_co2_in,
                    t_air_in=t_air_in,
                    t_co2_in=t_co2_in,
                    lower_bound=lower_bound,
                    upper_bound=midpoint,
                    max_depth=max_depth - 1,
                    upstream=upstream,
                )

    def energy_balance(
        self, p_co2_in: float, t_co2_in: float, t_co2_out: float, t_air_in: float
    ) -> int:
        """
        Check if the energy balance is satisfied. Compute the heat transfer coefficient from the geometry of the tube.
        Compare against the heat lost/gained by the CO2. This is used in the binary search.

        Args:
            p_co2_in (float): 
                Pressure of the initial co2 (depends on the direction of the flow).
            t_co2_in (float): 
                Temperature of the initial co2 (depends on the direction of the flow).
            t_co2_out (float): 
                Temperature of the final co2 (depends on the direction of the flow). This is what you guess.
            t_air_in (float): 
                Temperature of the initial air. This is always in the direction of air flow.

        Returns:
            int: 1 If the energy balance is satisfied, 2 if co2 emitted too little heat, 3 if co2 emitted too much heat.
        """
        p_co2_out = drop_pressure(p_in=p_co2_in)  # TODO make it give you delta P
        q_co2 = abs(
            energy_co2(
                p_co2_in, p_co2_out, t_co2_in, t_co2_out, constants.m_co2_per_tube
            )
        )
        t_air_out = q_co2 / (constants.m_air_segment * constants.cp_air) + t_air_in

        if t_air_out > min(t_co2_in, t_co2_out):
            # air out can't be hotter than sCO2 used to heat it
            # air is too hot -> decrease out temp
            # -> set current middle to be the new max
            # if self._verbose > 2:
            #     print("t_co2_out:", t_co2_out, "t_air_out:", t_air_out),
            right = t_co2_out
            return 3

        delta_t_m = lmtd(
            t_air_in=t_air_in,
            t_air_out=t_air_out,
            t_co2_in=t_co2_in,
            t_co2_out=t_co2_out,
        )
        # compute the ohtc
        t_air = (
            t_air_in + t_air_out
        ) / 2  # TODO: check if this is correct (LMTD instead?)
        p_co2 = (p_co2_in + p_co2_out) / 2
        t_co2 = (t_co2_in + t_co2_out) / 2
        p_air = constants.p_air_in  # TODO: compute air pressure drop
        ohtc = compute_ohtc(
            t_air=t_air,
            t_s=t_co2,
            p_air=p_air,
            p_s=p_co2,
            m_air=constants.m_air_segment,
            m_s=constants.m_co2_per_tube,
            tube=self.tube,
        )
        q_htc = ohtc * (delta_t_m)
        if np.isclose(q_htc, q_co2, rtol=constants.tolerance):
            return 1
        else:
            if q_htc > q_co2:
                return 2
            else:
                return 3

    def run(self):
        p_co2_init = constants.p_co2_outlet
        t_co2_init = constants.t_co2_outlet
        t_air_init = constants.t_air_inlet
        if self._verbose > 0:
            print("Initial conditions:")
            print("t_co2_in:", t_co2_init, "t_air_in:", t_air_init)
        # Start the solver from the CO2 outlet and air inlet and go upstream of the CO2 flow.
        tube_t_air, tube_t_co2, tube_p_co2 = self._solve_tube(
            p_co2_init=p_co2_init,
            t_co2_init=t_co2_init,
            t_air_init=t_air_init,
            upstream=True,
        )

        for _ in range(0):
            _ = self._solve_tube(
                p_co2_init=tube_p_co2[-1],
                t_co2_init=tube_t_co2[-1],
                t_air_init=None,
                upstream=False,
            )


        return tube_t_air, tube_t_co2, tube_p_co2


if __name__ == "__main__":
    from design import Tube
    from time import time

    tube = Tube()
    t0 = time()
    simulator = Simulator(tube, verbose=0)
    simulator.run()
    t1 = time()
    print("Time:", t1 - t0)

