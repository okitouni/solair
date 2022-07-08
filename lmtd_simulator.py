# Great job guys this looks good :)
import numpy as np
from typing import Iterable, Tuple
from CoolProp.CoolProp import PropsSI

# Designing for 100 MW


def get_enthalpy(
    p: Iterable[float], t: Iterable[float], fluid: str = "CO2"
) -> Iterable[float]:
    """Compute the enthalpy of a fluid at a given pressure and temperature.

    Args:
        p (Iterable[float]): Iterable containing the pressures of the fluid in Pa.
        t (Iterable[float]): Iterable containing the temperatures of the fluid in Kelvin.
        fluid (str, optional): Fluid to use. Defaults to "CO2".

    Returns:
        Iterable[float]: Iterable containing the enthalpy of the fluid at the given temperature and pressures in J/kg.
    """
    enthalpies = PropsSI("H", "P", p, "T", t, fluid)
    return enthalpies


def compute_air_flow_rate():
    pass


class constants:
    tolerance = 0.01  # tolerance in heat error to use in temperature search
    t_co2_inlet = 89 + 273.15
    t_co2_outlet = 35 + 273.15
    t_air_inlet = 25 + 273.15
    p_co2_inlet = 8e6
    p_co2_outlet = 7.48e6

    # calculate air mass flow rate
    # WRITE m_co2 calculation
    m_co2 = 950  # [kg/s] for 100MW
    cp_air = 1.005  # [kJ/kg-K]

    # calculate total Q exchanged
    # NEED ENTHALPY OF CO2
    h_co2_inlet = get_enthalpy(p_co2_inlet, t_co2_inlet)
    h_co2_outlet = get_enthalpy(p_co2_outlet, t_co2_outlet)
    total_heat_transferred = abs(m_co2 * (h_co2_outlet - h_co2_inlet))

    # find air mass flow rate
    t_air_outlet = 80 + 273.15  # guess outlet air T, cannot be higher than T_CO2_in
    delta_t_air = t_air_outlet - t_air_inlet
    m_air = total_heat_transferred / (
        cp_air * delta_t_air
    )  # Assuming something between 6.75 and 30 kg/s which correlates to about 2-5 m/s air velocity

    OHTC = 700  # [W/m2-K]

    n_tubes_in_row = 30
    n_rows = 4
    n_tubes_tot = n_tubes_in_row * n_rows


class Tube:
    def __init__(self):
        self.d_ot = 0.025  # [m]          outside tube diameter
        self.L_t = 9  # [m]  ##TOCHECK     length of tube
        self.t_f = 5e-4  #  [m]           fin thickness
        self.p_f = 2.8e-3  # [m]        fin pitch (distance between fins)
        self.n_f = np.floor(
            self.L_f / (self.t_f + self.p_f)
        )  #         number of fins (per tube at the moment)
        self.d_f = 57.15e-3  #  [m]     fin outside diameter
        self.d_r = 0.028  # [m] fin inside diameter

    def calculate_cross_section_air(self):
        cross_section_tube_s = np.pi * self.d_ot * (
            self.L_t - self.t_f * self.n_f
        ) + np.pi(
            (self.d_f ** 2 - self.d_ot ** 2) / 2
        )  # [m^2]
        return cross_section_tube_s

    def calculate_surface_area_air(self):
        surface_area_tube = np.pi * self.d_ot * (
            self.L_t - self.t_f * self.n_f
        ) + np.pi * ((d_f ** 2 - d_ot ** 2) / 2)
        return surface_area_tube


def calculate_re_co2(t, p, m):
    rho = PropsSI("D", "T", t, "P", p, "CO2")  # density
    mu = PropsSI("V", "T", t, "P", p, "CO2")  #
    di = 0.02  # internal diameter of single tube
    r = di / 2
    u = m / (rho * (np.pi * r ** 2))
    L = di
    Re = (rho * u * L) / mu
    return Re


def calculate_pr_co2(t, p):
    Pr = PropsSI("PRANDTL", "T", t, "P", p, "CO2")
    return Pr


def calculate_pr_air(t, p):
    Pr = PropsSI("PRANDTL", "T", t, "P", p, "AIR")
    return Pr


# WORKING ON THIS
def calculate_re_air(t, p, m, cross_section_air, tube: Tube):
    rho = PropsSI("D", "T", t, "P", p, "AIR")  # density
    mu = PropsSI("V", "T", t, "P", p, "AIR")  #
    wetter_perimeter = (
        2 * tube.L_t
        + 2 * (tube.t_f - tube.d_ot)
        + (tube.d_f - tube.d_ot) * 2 * tube.n_f
    )
    d_f = (
        4 * cross_section_air / wetter_perimeter
    )  #  hydraulic diameter defined in page 582 pdf Fundamentals of...
    a_ff = cross_section_air
    re_a = m * d_f / (a_ff * mu)
    return re_a


def calculate_htc_a(t, p, m, tube: Tube):
    """Compute OHTC
    Args:
        pr_a (float): Prandtl number.
        re_a (float): Reynolds number.
        tube (class element): class describing the geometry of the tube
    Returns:
        OHTC: overall heat transfer coefficient for air .
    """
    cross_section_air = tube.calculate_cross_section_air()
    pr_a = calculate_pr_air(t, p)
    re_a = calculate_re_air(t, p, m, cross_section_air, tube)
    k_a = 0.025  # W/(mK) thermal conductivity of air #TODO this varies over temperature and pressure
    htc_a = k_a / tube.d_ot(
        0.134
        * pr_a ** (1 / 3)
        * re_a ** 0.681
        * (2 * (tube.p_f - tube.t_f) ** 0.2 / (tube.d_f - tube.d_r))
        * ((tube.p_f - tube.t_f) / t_f) ** 0.1134
    )
    return htc_a


def calculate_htc_s(t, p, m):
    re_s = calculate_re_co2(p, t, m)
    pr_s = calculate_pr_co2(p, t)
    rho_pc = None  # change
    rho_s = None  # change
    t_pc = (
        -122.6 + 6.124 * p - 0.1657 * p ** 2 + 0.01773 * p ** 2.5 - 0.0005608 * p ** 3
    )
    if t > t_pc:
        a = 0.14
        b = 0.69
        c = 0.66
        n = 0
    else:
        a = 0.013
        b = 1.0
        c = -0.05
        n = 1.6
    htc_s = a * (re_s ** b) * (pr_s ** c) * (rho_pc / rho_s) ** n
    return htc_s


def compute_ohtc(t, p, re_s, pr_s, tube):
    # TODO do you need to multiply by n tubes in the z direction?
    surface_area_co2 = tube.calculate_surface_area_co2()
    surface_area_air = tube.calculate_surface_area_air()  # area of the tube
    htc_a = calculate_htc_a(tube)
    htc_s = calculate_htc_s(t, p)
    r_1 = 1 / (surface_area_air * htc_a)  # resistance of cold
    r_2 = 0  # negliged resistance of wall
    r_3 = 1 / (surface_area_co2 * htc_s)  # resistance of hot
    ohtc = 1 / (r_1 + r_2 + r_3)
    pass


def lmtd(t_air_in: float, t_air_out: float, t_co2_in: float, t_co2_out: float) -> float:
    """ " Compute the LMTD for a given set of temperatures.

    Args:
        th (Iterable[float]): Iterable containing the initial and final temperatures of the hot fluid.
        tc (Iterable[float]): Iterable containing the initial and final temperatures of the cold fluid.

    Returns:
        float: The log mean temperature difference.
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


# def lmtd2(t_air_in, t_air_out, t_co2_in, t_co2_out):
#     delta_t_in = t_co2_in - t_air_in
#     delta_t_out = t_co2_out - t_air_out
#     return (delta_t_in - delta_t_out) / np.log(delta_t_in / delta_t_out)


def energy_transferred(htc: float, th: Iterable[float], tc: Iterable[float]) -> float:
    """Compute the energy transferred between the hot and cold fluids.

    Args:
        htc (float): Heat transfer coefficient for the current iteration.
        th (Iterable[float]): Iterable containing the initial and final temperatures of the hot fluid.
        tc (Iterable[float]): Iterable containing the initial and final temperatures of the cold fluid.

    Returns:
        float: The energy transferred between the hot and cold fluids.
    """
    return htc * lmtd(th, tc)


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


def drop_pressure(p_in: float) -> float:
    """Compute the pressure drop of the sCO2 across an element and returns the outlet pressure.

    Args:
        p_in (_type_): Initial pressure of the sCO2.

    Returns:
        float: The output pressure of the sCO2.
    """
    # TODO figure out constants
    rho = 1
    f = 0.112  # from Eq. 16
    u = 1
    L = 1
    d = 1
    delta_p = (rho * f * u ** 2 * L) / (2 * d)
    p_out = p_in - delta_p
    return p_in  # TODO figure out how to compute pressure drop


def temperature_search(
    p_co2_in: float, t_co2_in: float, t_air_in: float, max_depth: int = 20, debug: bool = False
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
        m_air_segment = constants.m_air / n_segments
        # compute air out temp from heat transfer into air
        t_air_out = q_co2 / (m_air_segment * constants.cp_air) + t_air_in

        if t_air_out > t_co2_in:
            # air out can't be hotter than sCO2 used to heat it
            # air is too hot -> decrease out temp
            # -> set current middle to be the new max
            right = t_co2_out
            return binary_search_temps(left, right, max_depth - 1)

        delta_t_m = lmtd(t_air_in, t_air_out, t_co2_in, t_co2_out)
        q_htc = constants.OHTC * (delta_t_m)  # TODO OHTC is not constant
        if debug:
            print(
            f"t_co2_out: {t_co2_out}, t_air_out: {t_air_out}, q_co2: {q_co2}, q_htc: {q_htc}"
        )
        if abs((q_htc - abs(q_co2)) / q_co2) > constants.tolerance:
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
    t_co2_out, t_air_out = binary_search_temps(
        t_co2_in, constants.t_co2_inlet, max_depth
    )
    return t_co2_out, t_air_out


# print(
#     temperature_search(
#         p_co2_in=constants.p_co2_outlet, # start from the end of the tube
#         t_co2_in=constants.t_co2_inlet,
#         t_air_in=constants.t_air_inlet,
#     )
# )
#################


#### iteration over i (x-direction)

n_segments = 100

p_co2_in_new = constants.p_co2_outlet
t_co2_in_new = constants.t_co2_outlet
t_air_in_new = constants.t_air_inlet
print("Initial conditions:")
print("Mass flow rate of CO2:", constants.m_co2)
print("Mass flow rate of air:", constants.m_air)
print("n_segments:", n_segments)
print("p_co2_in:", p_co2_in_new, "t_co2_in:", t_co2_in_new, "t_air_in:", t_air_in_new)
for i in range(n_segments):
    t_co2_out, t_air_out = temperature_search(
        p_co2_in=p_co2_in_new, t_co2_in=t_co2_in_new, t_air_in=t_air_in_new, max_depth=50
    )

    p_co2_in_new = p_co2_in_new  # no pressure drop
    t_air_in_new = t_air_in_new  # same temperature of air at inlet (change later for other direction)
    t_co2_in_new = t_co2_out
    print("Finished segment:", i)
    print(f"t_co2_out: {t_co2_out}, t_air_out: {t_air_out}")
