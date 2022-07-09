from typing import Iterable
from CoolProp.CoolProp import PropsSI

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

def drop_pressure(p_in: float) -> float:
    """Compute the pressure drop of the sCO2 across an element and returns the outlet pressure.

    Args:
        p_in (_type_): Initial pressure of the sCO2.

    Returns:
        float: The output pressure of the sCO2.
    """
    if True:
        return p_in
    # TODO add pressure constants
    rho = 1
    f = 0.112  # from Eq. 16
    u = 1
    L = 1
    d = 1
    delta_p = (rho * f * u ** 2 * L) / (2 * d)
    p_out = p_in - delta_p
    return p_out  # TODO how to compute pressure drop