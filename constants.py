from utils import get_enthalpy
from dataclasses import dataclass


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


@dataclass
class constants:
    # solver tolerance
    tolerance: float = 0.01  # tolerance in heat error to use in temperature search
    # design constants
    n_segments: int = 100
    n_tubes_in_row: int = 30
    n_rows: int = 4
    n_tubes_tot: int = n_tubes_in_row * n_rows

    # thermodynamic constants
    t_co2_inlet: float = 72.83 + 273.15
    t_co2_outlet: float = 33 + 273.15
    p_co2_inlet: float = 8e6
    p_co2_outlet: float = 7.48e6
    t_air_inlet: float = 25 + 273.15
    t_air_outlet: float = 60 + 273.15  # guess outlet air T, cannot be higher than T_CO2_in
    p_air_in: float = 101000

    # find air mass flow rate
    cp_air: float = 1005.0  # [J/kg-K]
    m_co2: float = 58.34625  # 443.177  # [kg/s] for 100MW
    m_air: float = get_m_air(
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
    m_co2_segment: float = m_co2 /(n_tubes_in_row*n_rows)
