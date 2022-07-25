from dataclasses import dataclass
import pickle
import warnings
from CoolProp.CoolProp import PropsSI
from pathlib import Path

def get_m_air(
    p_co2_inlet: float,
    p_co2_outlet: float,
    t_co2_inlet: float,
    t_co2_outlet: float,
    t_air_inlet: float,
    t_air_outlet: float,
    m_co2: float,
    cp_air: float,
) -> float:
    """Compute the mass flow rate of air for the entire cooler from the properties of the thermodynamic cycle.

    Args:
        p_co2_inlet (float): Inlet pressure of the sCO2.
        p_co2_outlet (float): Outlet pressure of the sCO2.
        t_co2_inlet (float): Inlet temperature of the sCO2.
        t_co2_outlet (float): Outlet temperature of the sCO2.
        t_air_inlet (float): Inlet temperature of the air.
        t_air_outlet (float): Outlet temperature of the air.
        m_co2 (float): Mass flow rate of the sCO2.
        cp_air (float): Specific heat capacity of the air.

    Returns:
        float: The mass flow rate of air for the entire cooler.
    """
    # calculate total Q exchanged
    _h_co2_inlet = get_enthalpy(p_co2_inlet, t_co2_inlet)
    _h_co2_outlet = get_enthalpy(p_co2_outlet, t_co2_outlet)
    _total_heat_transferred = abs(m_co2 * (_h_co2_outlet - _h_co2_inlet))

    _delta_t_air = t_air_outlet - t_air_inlet
    # Assuming something between 6.75 and 30 kg/s which correlates to about 2-5 m/s air velocity
    # compute air out temp from heat transfer into air
    m_air = _total_heat_transferred / (cp_air * _delta_t_air)
    return m_air


def get_m_air_segment(m_air: float, n_segments: int, n_tubes_in_row: int) -> float:
    """Compute the mass flow rate of air for a segment.

    Args:
        m_air (float): Mass flow rate of air for the entire cooler.
        n_segments (int): Number of segments along the tube.
        n_tubes_in_row (int): Number of tubes in a row.

    Returns:
        float: The mass flow rate of air for a segment.
    """
    # TODO is m_air for segment corrent?
    m_air_segment = (
        m_air / n_segments / n_tubes_in_row
    )  # divided by number of tubes in a row as we are considering a single tube for all calculations
    return m_air_segment


# path = Path(__file__).parent / "sco2_enthalpies.pkl"
# get_enthalpy_pickle = pickle.load(open(path, "rb"))


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
        raise NotImplementedError("Fast enthalpy interpolation not implemented yet.")
    if fluid.lower() != "co2" and fast:
        warnings.warn(
            f"Provided fluid: {fluid}. Only CO2 is supported for fast enthalpy calculation. Using PropsSI. Set fast to False to silence this warning."
        )
    if fast:
        if t > 345.98 + 50 or t < 306.15:  # TODO make this a constant
            warnings.warn(
                f"Temperature outside of range of enthalpy table. Acceptable range is 306.15 to 345.98 K.",
                stacklevel=2
            )
        if p > 8e6 or p < 7.48e6:
            warnings.warn(
                f"Pressure outside of range of enthalpy table. Acceptable range is 7.48e6 to 8e6 Pa.",
                stacklevel=2
            )
        # h = get_enthalpy_pickle(p, t)
        if len(h) == 1:
            return h[0]
    else:
        h = PropsSI("H", "P", p, "T", t, fluid)
    return h


@dataclass
class constants:
    # solver tolerance
    tolerance: float = 0.01  # tolerance in heat error to use in temperature search
    # design constants
    n_segments: int = 30  # done
    n_tubes_in_row: int = 47.5  # done
    n_rows: int = 4  # done
    n_tubes_tot: int = n_tubes_in_row * n_rows  # done 190

    # thermodynamic constants
    t_co2_inlet: float = 71 + 273.15  # done
    t_co2_outlet: float = 40.3 + 273.15  # done
    p_co2_inlet: float = 7.5e6  # done
    p_co2_outlet: float = 7.4999e6  # done
    t_air_inlet: float = 20 + 273.15  # done
    t_air_outlet: float = 43.4 + 273.15  # done
    p_air_in: float = 99695  # done

    # find air mass flow rate
    cp_air: float = 1005.0  # [J/kg-K]
    m_co2: float = 406.6  # done
    m_air: float = get_m_air(  # must equal 1107
        p_co2_inlet,
        p_co2_outlet,
        t_co2_inlet,
        t_co2_outlet,
        t_air_inlet,
        t_air_outlet,
        m_co2,
        cp_air,
    )
    m_air_segment: float = get_m_air_segment(m_air, n_segments, n_tubes_in_row)
    m_co2_segment: float = m_co2 / (n_tubes_in_row * n_rows)
