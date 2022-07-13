from utils import get_enthalpy
from dataclasses import dataclass

@dataclass
class constants:
    # solver tolerance
    tolerance = 0.01  # tolerance in heat error to use in temperature search
    # design constants
    n_segments = 100
    n_tubes_in_row = 30
    n_rows = 4
    n_tubes_tot = n_tubes_in_row * n_rows

    # OHTC = 100  # [W/m2-K]

    # thermodynamic constants
    t_co2_inlet = 72.83 + 273.15
    t_co2_outlet = 33 + 273.15
    p_co2_inlet = 8e6
    p_co2_outlet = 7.48e6
    t_air_inlet = 25 + 273.15
    t_air_outlet = 60 + 273.15  # guess outlet air T, cannot be higher than T_CO2_in
    p_air_in = 101000



    # find air mass flow rate
    cp_air = 1005.0  # [J/kg-K]
    m_co2 = 443.177  # [kg/s] for 100MW

    @property
    def m_air(self):
        # calculate total Q exchanged
        # NEED ENTHALPY OF CO2
        _h_co2_inlet = get_enthalpy(self.p_co2_inlet, self.t_co2_inlet)
        _h_co2_outlet = get_enthalpy(self.p_co2_outlet, self.t_co2_outlet)
        _total_heat_transferred = abs(self.m_co2 * (_h_co2_outlet - _h_co2_inlet))

        _delta_t_air = self.t_air_outlet - self.t_air_inlet
        m_air = _total_heat_transferred / (
            self.cp_air * _delta_t_air
        )  # Assuming something between 6.75 and 30 kg/s which correlates to about 2-5 m/s air velocity
        # compute air out temp from heat transfer into air
        self._m_air = m_air
        return self.m_air

    def get_m_air_segment(self, m_air):
        # TODO is m_air for segment corrent?
        m_air_segment = (
        constants.m_air / constants.n_segments / constants.n_tubes_in_row
        )  # divided by number of tubes in a row as we are considering a single tube for all calculations   
        return m_air_segment  

if __name__ == "__main__":
    print(constants.m_air)