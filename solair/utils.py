from .design import Tube, calculate_re_co2
from .constants import get_enthalpy
from CoolProp.CoolProp import PropsSI
import numpy as np


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
            "have: t_co2_out = {}, t_co2_in = {}, t_air_out = {}, t_air_in = {}\n"
            "This error is usually thrown when the model converges and t_co2_out == t_co2_in.".format(
                t_co2_out, t_co2_in, t_air_out, t_air_in
            )
        )
    return delta / np.log((t_co2_in - t_air_out) / (t_co2_out - t_air_in))


def energy_co2(
    p_in: float, p_out: float, t_in: float, t_out: float, m_co2: float, fast=True
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
        fast (bool, optional):
            Use the fast enthalpy function. Defaults to True.

    Returns:
        float: The energy transferred out of the sCO2.
    """
    h_in = get_enthalpy(p_in, t_in, fast=fast)
    h_out = get_enthalpy(p_out, t_out, fast=fast)
    q_co2 = m_co2 * abs(h_out - h_in)
    return q_co2


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


def drop_pressure(p_in: float, t: float, m: float, tube: Tube) -> float:
    """Compute the pressure drop of the sCO2 across an element and returns the outlet pressure.

    Args:
        p_in (_type_): Initial pressure of the sCO2.

    Returns:
        float: The output pressure of the sCO2.
    """
    # if True:
    # return p_in
    # _delta_pressure = 1 - np.exp(
    # np.log(1 - 0.375) / 400
    # )   # drop pressure by 37.5% across the whole length 100 segments * 4 SHX.

    # return 0.00001 * p_in  # _delta_pressure * p_in
    # Things are still not working below here
    # TODO explain constants and update docstring
    epsilon = 0.045
    rho = PropsSI("D", "T", t, "P", p_in, "CO2")
    Re_co2 = calculate_re_co2(t, p_in, m, tube)
    z = epsilon / tube.tube_in_diameter
    f = 8 * (
        (8 / Re_co2) ** (12)
        + (
            2.457 * np.log(1 / ((7 / Re_co2) ** 0.9) + 0.27 * z) ** 16
            + (37530 / Re_co2) ** 16
        )
        ** (-3 / 2)
    ) ** (1 / 12)
    di = tube.tube_in_diameter  # internal diameter of single tube
    r = di / 2
    u = m / (rho * (np.pi * r ** 2))
    L = tube.segment_length # TODO check that this is correct
    delta_p = (rho * f * u ** 2 * L) / (2 * di)
    return delta_p
