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


def get_m_air_segment(
    m_air: float, n_segments: int, n_tubes_in_cross_section: float
) -> float:
    """Compute the mass flow rate of air for a segment.

    Args:
        m_air (float): Mass flow rate of air for the entire cooler.
        n_segments (int): Number of segments along the tube.
        n_tubes_in_cross_section (float): Number of tubes facing the air inlet (surface perpendicular to Air inlet).

    Returns:
        float: The mass flow rate of air for a segment.
    """
    m_air_segment = (
        m_air / n_segments / n_tubes_in_cross_section
    )  # divided by number of tubes in a row as we are considering a single tube for all calculations
    return m_air_segment


# path = Path(__file__).parent / "sco2_enthalpies.pkl"
# get_enthalpy_pickle = pickle.load(open(path, "rb"))


def get_enthalpy(p: float, t: float, fluid: str = "CO2") -> float:
    """
    Compute the enthalpy of a fluid at a given pressure and temperature. Using scipy interpolation from a pickle file.

    Args:
        p (float): 
            Pressure of the fluid in Pa.
        t (float): 
            Temperature of the fluid in Kelvin.
        fluid (str, optional): 
            Fluid to use. Defaults to "CO2". Not actually used here but kept for consistency.

    Returns:
        float: The enthalpy of the fluid at the given temperature and pressure in J/kg.
    """
    h = PropsSI("H", "P", p, "T", t, fluid)
    return h


@dataclass
class constants:
    # TODO fix inconsistent dataclass and init
    def __init__(self, t_air_inlet=20, t_air_out=20):
        # solver tolerance
        self.tolerance: float = 0.0001  # tolerance in heat error to use in temperature search
        # design constants
        self.n_segments: int = 30  # done
        self.n_tubes_in_row: float = 47.5  # done
        self.n_rows: int = 4  # done
        self.n_bundles: int = 49  # done
        self.n_tubes_in_cross_section: float = self.n_tubes_in_row * self.n_bundles  # TODO this is different for different designs
        self.n_tubes_tot: float = self.n_tubes_in_row * self.n_rows * self.n_bundles  # done
        self.n_sub_heat_exchangers: int = 3 # TODO not used in the simulator currently


        # thermodynamic constants
        self.t_co2_inlet: float = 71 + 273.15  # done
        self.t_co2_outlet: float = 40.3 + 273.15  # done
        self.p_co2_inlet: float = 7.5e6  # done
        self.p_co2_outlet: float = 7.4999e6  # done
        self.t_air_inlet: float = t_air_inlet + 273.15  # done
        self.t_air_outlet: float = t_air_out + self.t_air_inlet #43.4 + 273.15  # done
        self.p_air_in: float = 99695  # done

        # find air mass flow rate
        self.cp_air: float = 1005.0  # [J/kg-K]
        self.m_co2: float = 406.6  # done
        self.m_air: float = get_m_air(  # must equal 1107
            self.p_co2_inlet,
            self.p_co2_outlet,
            self.t_co2_inlet,
            self.t_co2_outlet,
            self.t_air_inlet,
            self.t_air_outlet,
            self.m_co2,
            self.cp_air,
        )
        self.m_air_segment: float = get_m_air_segment(
            self.m_air, self.n_segments, n_tubes_in_cross_section=self.n_tubes_in_cross_section
        )
        self.m_co2_segment: float = self.m_co2 / (self.n_tubes_tot)

        # parameters for cost calculation
        self.cost_alu: float = 4.2   # [$/kg] material cost aluminium
        self.cost_steel: float = 0.8 # [$/kg] material cost steel (ASTM A214)
        self.rho_alu: float = 2750   # [kg/m^3] density aluminium
        self.rho_steel: float = 7950 # [kg/m^3] density steel (ASTM A214)

        self.fixed_cost_tube: float = 2 # [$/m] fixed cost of the tube per unit length (Kroger vol 2)
        self.weighting_factor: float = 2 # [$/m] weighting factor (kroger vol 2)

        self.f_header = 0.8     # weighting factor for heat exchanger header cost
        self.f_labor = 0.7      # weighting factor for heat exchanger labor cost
        self.f_HX = 1.2         # ?
        self.lifetime_years: float = 25
        self.LCOE_fanpower_cents: float = 0.05 
