from ast import Raise
from typing import Iterable
from CoolProp.CoolProp import PropsSI
import numpy as np
from design import Tube, calculate_re_co2
from constants import get_enthalpy
import pickle
import warnings


def lmtd(t_air_in: float, t_air_out: float, t_co2_in: float, t_co2_out: float) -> float:
    """
    Compute the LMTD in a segment.

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
    """
    Compute the energy transferred out of the sCO2.

    Args:
        p_in (float): 
            Initial pressure of the sCO2.
        p_out (float): 
            Final pressure of the sCO2.
        t_in (float): 
            Initial temperature of the sCO2.
        t_out (float): 
            Final temperature of the sCO2.

    Returns:
        float: The energy transferred out of the sCO2.
    """
    h_in = get_enthalpy(p_in, t_in)
    h_out = get_enthalpy(p_out, t_out)
    return m_co2 * abs(h_out - h_in)


# def energy_air(mdot: float, t0: float, t1: float) -> float:
#     """Compute the energy transferred into the air.

#     Args:
#         mdot (float): Mass flow rate of the air.
#         t0 (float): Temperature of the inlet air.
#         t1 (float): Temperature of the outlet air.

#     Returns:
#         float: The energy transferred into the air.
#     """
#     return constants.cp_air * mdot * (t1 - t0)

get_enthalpy_pickle = pickle.load(open("sco2_enthalpies.pkl", "rb"))


def get_enthalpy(p: float, t: float, fluid: str = "CO2", fast=True) -> float:
    """
    Compute the enthalpy of a fluid at a given pressure and temperature. Using scipy interpolation from a pickle file.

    Args:
        p (float): 
            Pressure of the fluid in Pa.
        t (float): 
            Temperature of the fluid in Kelvin.
        fluid (str, optional): 
            Fluid to use. Defaults to "CO2". Not actually used here but kept for consistency.
        fast (bool, optional): 
            Use the scipy interpolation version of the function. Defaults to True. 
            If False, the PropsSI function is used.

    Returns:
        float: The enthalpy of the fluid at the given temperature and pressure in J/kg.
    """
    if fast:
        if fluid.lower() != "co2":
            Raise(ValueError("Only CO2 is supported for fast enthalpy. Got {}".format(fluid)))
        if t > 345.98 or t < 306.15:
            warnings.warn(
                f"Temperature {t} outside of range of enthalpy table. Acceptable range is 306.15 to 345.98 K."
            )
        if p > 8e6 or p < 7.48e6:
            warnings.warn(
                f"Pressure {p} outside of range of enthalpy table. Acceptable range is 7.48e6 to 8e6 Pa."
            )
        h = get_enthalpy_pickle(p, t)
        if len(h) == 1:
            return h[0]
    else:
        h = PropsSI("H", "P", p, "T", t, fluid)
    return h


def drop_pressure(p_in: float) -> float:
    """
    Compute the pressure drop of the sCO2 across an element and returns the outlet pressure.

    Args:
        p_in (_type_): Initial pressure of the sCO2.

    Returns:
        float: The output pressure of the sCO2.
    """
    # if True:
    # return p_in
    # TODO add pressure constants
    rho = PropsSI("D", "T", t, "P", p_in, "CO2")
    Re_co2 = calculate_re_co2(t, p_in, m)
    z = 1  # TODO check reference (can't find it)
    f = 8 * (
        (8 / Re_co2) ** (12)
        + (
            2.457 * np.log(1 / ((7 / Re_co2) ** 0.9) + 0.27 * z) ** 16
            + (37530 / Re_co2) ** 16
        )
        ** (-3 / 2)
    ) ** (1 / 12)
    di = tube.d_i  # internal diameter of single tube
    r = di / 2
    u = m / (rho * (np.pi * r ** 2))
    L = tube.L_t
    d = tube.d_i
    delta_p = (rho * f * u ** 2 * L) / (2 * d)
    p_out = p_in - delta_p
    return p_out
