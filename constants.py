from utils import get_enthalpy

class constants:
    # solver tolerance
    tolerance = 0.01  # tolerance in heat error to use in temperature search
    # design constants
    n_segments = 100
    n_tubes_in_row = 30
    n_rows = 4
    n_tubes_tot = n_tubes_in_row * n_rows
    
    OHTC = 1000  # [W/m2-K]

    # thermodynamic constants
    t_co2_inlet = 72.83 + 273.15
    t_co2_outlet = 33 + 273.15
    p_co2_inlet = 8e6
    p_co2_outlet = 7.48e6
    t_air_inlet = 25 + 273.15
    t_air_outlet = 80 + 273.15  # guess outlet air T, cannot be higher than T_CO2_in

    # calculate air mass flow rate
    # WRITE m_co2 calculation
    m_co2 = 443.177  # [kg/s] for 100MW
    cp_air = 1.005  # [kJ/kg-K]

    # calculate total Q exchanged
    # NEED ENTHALPY OF CO2
    _h_co2_inlet = get_enthalpy(p_co2_inlet, t_co2_inlet)
    _h_co2_outlet = get_enthalpy(p_co2_outlet, t_co2_outlet)
    _total_heat_transferred = abs(m_co2 * (_h_co2_outlet - _h_co2_inlet))

    # find air mass flow rate
    _delta_t_air = t_air_outlet - t_air_inlet
    m_air = _total_heat_transferred / (
        cp_air * _delta_t_air
    )  # Assuming something between 6.75 and 30 kg/s which correlates to about 2-5 m/s air velocity

