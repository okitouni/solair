from typing import Dict, Iterable, List, Tuple, Union
from constants import constants
import numpy as np
from utils import drop_pressure, energy_co2, lmtd
from design import compute_ohtc, Tube
from collections import defaultdict

class Simulator:
    def __init__(
        self,
        tube: Tube,
        verbose: int = 0,
        max_iterations: int = 100,
        n_rows: int = constants.n_rows,
        fast: bool = False,
    ):
        """
        Simulator for the full cooler (all Sub-Heat Exchangers SHX). Initialized with a tube object which contains
        the properties of an individual tube in the heat exchanger.
        To use the simulator:
        simulator = Simulator(tube)
        simulator.run()
        results = simulator.results

        Args:
            tube (Tube): The tube object to use in the simulator. See design.py for details. 
            verbose (int, optional): 
                Level of verbosity of the simulator. Options are 0 for quiet, 1 SHX information, 2 tube level information,
                and 3 for full debug at every segement.
                Defaults to 0.
            max_iterations (int, optional): Maximum number of iterations for the binary search. Defaults to 100.
            n_rows (int, optional): Number of rows of tubes per heat exchanger. Defaults to constants.n_rows.
            fast (bool, optional): Whether to use fast scipy interpolation or CoolProp. Defaults to False.
        """    
        self.tube = tube
        self.results = defaultdict(list)
        self.converged = False
        self.n_segments = 100
        self._verbose = verbose
        self.max_iterations = max_iterations
        self.n_rows = n_rows
        self.fast = fast

    def _solve_segment(
        self,
        p_co2_in: float,
        t_co2_in: float,
        t_air_in: float,
        upstream: bool = False,
        max_temp: float = constants.t_co2_inlet + 100,
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
            max_depth=self.max_iterations,
        )
        # solve for t_air_out (AGAIN) # TODO remove this double calculations
        delta_p = drop_pressure(
            p_co2_in, t_co2_in, constants.m_co2_segment, self.tube
        )  # TODO calculation uses incorrect z
        delta_p = delta_p if upstream else -delta_p
        p_co2_out = p_co2_in + delta_p

        q_co2 = abs(
            energy_co2(
                p_co2_in,
                p_co2_out,
                t_co2_in,
                t_co2_out,
                constants.m_co2_segment,
                fast=self.fast,
            )
        )
        t_air_out = q_co2 / (constants.m_air_segment * constants.cp_air) + t_air_in
        return t_air_out, t_co2_out, p_co2_out

    def _solve_tube(
        self,
        p_co2_init: float,
        t_co2_init: float,
        t_air_init: Union[float, Iterable[float]] = None,
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
            t_air_in_list = self.results["t_air"][-1]
        elif isinstance(t_air_init, Iterable):
            t_air_in_list = t_air_init
            if len(t_air_in_list) != self.n_segments:
                raise ValueError(
                    "Length of t_air_init must be equal to n_segments. Got {} instead.".format(
                        len(t_air_in_list)
                    )
                )
        else:
            t_air_in_list = [t_air_init] * self.n_segments

        if upstream:
            # going upstream so reverse the list which is going downstream
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
        if upstream:
            # arrays are stored in the same direction as the flow
            tube_t_co2.reverse()
            tube_t_air.reverse()
            tube_p_co2.reverse()
        return tube_t_air, tube_t_co2, tube_p_co2

    def binary_search(
        self,
        p_co2_in: float,
        t_air_in: float,
        t_co2_in: float,
        lower_bound: float,
        upper_bound: float,
        max_depth: int = None,
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
                Maximum number of iterations. Defaults to self.max_iterations.
            upstream (bool, optional): True 
                If solver is going upstream of CO2 flow, False if going downstream. Defaults to False.

        Returns:
            float: temperature of the final CO2.
        """
        delta_p = drop_pressure(
            p_co2_in, t_co2_in, constants.m_co2_segment, self.tube
        )  # TODO calculation uses incorrect z
        delta_p = delta_p if upstream else -delta_p
        p_co2_out = p_co2_in + delta_p

        midpoint = (lower_bound + upper_bound) / 2
        if max_depth == 0:
            return midpoint
        # if np.isclose(lower_bound, upper_bound, rtol=1e-6):
        #     self.converged = True
        #     return midpoint

        if upstream:
            out = self.energy_balance(
                p_co2_in=p_co2_in,
                p_co2_out=p_co2_out,
                t_co2_in=t_co2_in,
                t_co2_out=midpoint,
                t_air_in=t_air_in,
            )
        else:
            out = self.energy_balance(
                p_co2_in=p_co2_in,
                p_co2_out=p_co2_out,
                t_co2_in=midpoint,
                t_co2_out=t_co2_in,
                t_air_in=t_air_in,
            )
        # THIS IS A HACK
        if isinstance(out, Iterable):
            (q_htc, q_co2) = out
        else:
            # too much heat from co2
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

        if self._verbose > 2:
            print(
                f"depth: {max_depth}, q_co2: {q_co2:.2f} q_htc: {q_htc:.2f}, t_co2_out: {midpoint:.2f}, p_co2_in: {p_co2_in:.2e}, p_co2_out: {p_co2_out:.2e}"
            )

        if np.isclose(q_htc, q_co2, rtol=constants.tolerance):
            return midpoint
        else:
            if q_htc > q_co2:
                # too little heat from co2
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
            else:
                # too much heat from co2
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
        self,
        p_co2_in: float,
        p_co2_out: float,
        t_co2_in: float,
        t_co2_out: float,
        t_air_in: float,
    ) -> int:
        """
        Check if the energy balance is satisfied. Compute the heat transfer coefficient from the geometry of the tube.
        Compare against the heat lost/gained by the CO2. This is used in the binary search.

        Args:
            p_co2_in (float): 
                Pressure of the initial co2 (depends on the direction of the flow).
            p_co2_out (float): 
                Pressure of the final co2 (depends on the direction of the flow).
            t_co2_in (float): 
                Temperature of the initial co2 (depends on the direction of the flow).
            t_co2_out (float): 
                Temperature of the final co2 (depends on the direction of the flow). This is what you guess.
            t_air_in (float): 
                Temperature of the initial air. This is always in the direction of air flow.

        Returns:
            int: 1 If the energy balance is satisfied, 2 if co2 emitted too little heat, 3 if co2 emitted too much heat.
        """

        q_co2 = abs(
            energy_co2(
                p_co2_in,
                p_co2_out,
                t_co2_in,
                t_co2_out,
                constants.m_co2_segment,
                fast=self.fast,
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
            m_s=constants.m_co2_segment,
            tube=self.tube,
        )
        q_htc = ohtc * (delta_t_m)
        return (q_htc, q_co2)

    def _start_shx(self, n_rows: int = None) -> Tuple[List]:
        """
        Simulate the first Sub-Heat Exchanger (SHX) in the system 
        (this will be the closes SHX to the CO2 outlet).

        Args:
            n_rows (int, optional): Number of rows of tubes in the SHX. Defaults to constants.n_rows.

        Returns:
            Tuple[List]: Tuple containing the temperature of air and 
            the temperature and pressure of the CO2 in the last tube. 
        """
        n_rows = n_rows or self.n_rows
        for i in range(n_rows):
            if i == 0:
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
            else:
                tube_t_air, tube_t_co2, tube_p_co2 = self._solve_tube(
                    p_co2_init=tube_p_co2[0],
                    t_co2_init=tube_t_co2[0],
                    t_air_init=None,
                    upstream=False,
                )
            # At each row, save the temperature and pressure of the air and CO2 in the tube.
            self.results["t_co2"].append(tube_t_co2)
            self.results["t_air"].append(tube_t_air)
            self.results["p_co2"].append(tube_p_co2)

        return self.results

    def _intermediate_shx(
        self, n_rows: int = None, guess_t_co2: float = None
    ) -> Dict:
        """
        Solve an intermediate Sub-Heat Exchanger (SHX). 
        The difference between this and the start SHX is that the initial conditions of the SHX 
        (outlet temperature and pressure) have to be guessed and corrected for in an iterative way.

        Args:
            n_rows (int, optional): 
                Number of rows of tubes in the SHX. Defaults to constants.n_rows.
            guess_t_co2 (float, optional): 
                Guess for the outlet temperature of the last row before the next SHX. 
                Defaults to the temperature at the previous row - 1.

        Returns:
            Dict: Dictionary with the CO2 and air properties at the SHX.
        """
        n_rows = n_rows or self.n_rows
        results = defaultdict(list)
        for i in range(n_rows):
            if i == 0:
                # guess the same temperature as the input (flow-wise) to previous shx
                guess_t_co2 = guess_t_co2 or self.results["t_co2"][-1][0] - 1
                guess_p_co2 = self.results["p_co2"][-1][0]
                t_air_init = self.results["t_air"][-1]
                # Begin solving from the outlet of current shx and go upstream of the CO2 flow.
                # if self._verbose > 0:
                #     print("Solving next SHX. Guess for initial conditions:")
                #     print("t_co2_in:", guess_t_co2, "p_co2_in:", guess_p_co2)
                # Start the solver from the CO2 outlet and air inlet and go upstream of the CO2 flow.
                tube_t_air, tube_t_co2, tube_p_co2 = self._solve_tube(
                    p_co2_init=guess_p_co2,
                    t_co2_init=guess_t_co2,
                    t_air_init=t_air_init,
                    upstream=True,
                )
            else:
                tube_t_air, tube_t_co2, tube_p_co2 = self._solve_tube(
                    p_co2_init=tube_p_co2[0],
                    t_co2_init=tube_t_co2[0],
                    t_air_init=None,
                    upstream=False,
                )
            # At each row, save the temperature and pressure of the air and CO2 in the tube.
            results["t_co2"].append(tube_t_co2)
            results["t_air"].append(tube_t_air)
            results["p_co2"].append(tube_p_co2)

        return results

    def run(self):
        """ Run the entire simulation.
        """
        results = self._start_shx()
        if self._verbose > 0:
            print(f"Solving SHX 1.")
            print(f"Input air temp: {constants.t_air_inlet:.2f}")
            print(f"Target CO2 temp: {constants.t_co2_outlet:.2f}")
            shx_t_co2_outlet = [
                    tube_temperatures[-1] for tube_temperatures in results["t_co2"]
                ]
            print(
                    "CO2 outlet temp max/mean/min: %.2f/ %.2f/ %.2f"
                    % (
                        np.max(shx_t_co2_outlet),
                        np.mean(shx_t_co2_outlet),
                        np.min(shx_t_co2_outlet),
                    )
                )

        for i in range(1, 4):
            if i == 1:
                results = self._start_shx()
                average_air_t_in = constants.t_air_inlet
                required_shx_t_co2_outlet = constants.t_co2_outlet
            else:
                converged = False
                left = required_shx_t_co2_outlet - 2
                right = required_shx_t_co2_outlet
                midpoint = required_shx_t_co2_outlet
                while not converged:
                    results = self._intermediate_shx(guess_t_co2=midpoint)
                    shx_t_co2_outlet = [
                        tube_temperatures[-1] for tube_temperatures in results["t_co2"]
                    ]
                    mean_shx_t_co2_outlet = np.mean(shx_t_co2_outlet)
                    converged = np.isclose(
                        required_shx_t_co2_outlet,
                        mean_shx_t_co2_outlet,
                        rtol=constants.tolerance,
                    )
                    if required_shx_t_co2_outlet > mean_shx_t_co2_outlet:
                        left = midpoint
                    else:
                        right = midpoint
                    midpoint = (left + right) / 2
                    if self._verbose > 1:
                        print("SHX t_co2_outlet:", shx_t_co2_outlet)
                        print("Required SHX t_co2_outlet:", required_shx_t_co2_outlet)
                        print("Mean SHX t_co2_outlet:", mean_shx_t_co2_outlet)
                        print(
                            "Percent Error:",
                            (required_shx_t_co2_outlet - mean_shx_t_co2_outlet)
                            / (required_shx_t_co2_outlet + mean_shx_t_co2_outlet)
                            * 200,
                        )
                    
                average_air_t_in = np.mean(results["t_air"][-1])
                required_shx_t_co2_outlet = results["t_co2"][-1][0]

            if self._verbose > 0:
                print(f"SHX {i} converged.")
                print(f"Average input air temp: {average_air_t_in}")
                print(f"Target outlet CO2 temp: {required_shx_t_co2_outlet:.2f}")
                print(
                    "CO2 outlet temp max/mean/min: %.2f/ %.2f/ %.2f"
                    % (
                        np.max(shx_t_co2_outlet),
                        np.mean(shx_t_co2_outlet),
                        np.min(shx_t_co2_outlet),
                    )
                )
            for key, value in results.items():
                for property_list in value:
                    self.results[key].append(property_list)
        if self._verbose > 0:
            print("Simulation complete.")
            print("CO2 inlet tempreature:", self.results["t_co2"][-1][0])
            print("Average air outlet temperature:", np.mean(self.results["t_air"][-1]))


if __name__ == "__main__":
    from time import time

    tube = Tube()
    t0 = time()
    simulator = Simulator(tube, verbose=1, max_iterations=100)
    simulator.run()
    t1 = time()
    print("Time:", t1 - t0)

