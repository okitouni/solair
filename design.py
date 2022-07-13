import numpy as np
from typing import Iterable
from CoolProp.CoolProp import PropsSI
from constants import constants

class Tube:
    def __init__(self):
        self.d_ot = 0.025  # [m]        outside tube diameter
        self.L_t = 9 / constants.n_segments # [m]  length of tube / number of segments = length of segment
        self.t_f = 5e-4    # [m]        fin thickness
        self.p_f = 2.8e-3  # [m]        fin pitch (distance between fins)
        self.tp = 0.058    # [m]       transversal tube pitch
        self.n_f = np.floor(
            self.L_t / (self.t_f + self.p_f)
        )                  #         number of fins (per segment)
        self.d_f = 57.15e-3  #  [m]     fin outside diameter
        self.d_r = 0.028 # [m] fin inside diameter
        self.d_i = 20e-3 # [m] tube inside diameter
        #self.h_a = 
        #self.k_f = 
    def calculate_cross_section_air(self):
        cross_section_tube_s = (self.tp-self.d_ot)*self.L_t-(self.d_f-self.d_ot)*self.t_f*self.n_f
        return cross_section_tube_s

    def calculate_surface_area_air(self):   # surface area of entire tube
        A_1f = np.pi * ((self.d_f ** 2 - self.d_ot ** 2) / 2) # surface area of one fin
        A_f = self.n_f * A_1f   # surface area of all fins along one tube
        A_t = np.pi * self.d_ot * (
            self.L_t - self.t_f * self.n_f) # outside tube area between fins for entire tube
        return [A_t,A_f]

    def calculate_surface_area_co2(self):
        surface_area_co2 = 2*np.pi*self.d_i*self.L_t 
        return surface_area_co2


def calculate_pr_co2(t, p):
    Pr = PropsSI("PRANDTL", "T", t, "P", p, "CO2")
    return Pr


def calculate_pr_air(t, p):
    Pr = PropsSI("PRANDTL", "T", t, "P", p, "AIR")
    return Pr


def calculate_re_air(t, p, m, cross_section_air, tube: Tube):
    rho = PropsSI("D", "T", t, "P", p, "AIR")  # density
    mu = PropsSI("V", "T", t, "P", p, "AIR")  # viscosity
    wetted_perimeter = (
        2 * tube.L_t
        + 2 * (tube.t_f - tube.d_ot)
        + (tube.d_f - tube.d_ot) * 2 * tube.n_f
    )
    d_f = (
        4 * cross_section_air / wetted_perimeter
    )  #  hydraulic diameter defined in page 582 pdf Fundamentals of...
    a_ff = cross_section_air
    re_a = m * d_f / (a_ff * mu)
    return re_a


def calculate_re_co2(t, p, m):
    rho = PropsSI("D", "T", t, "P", p, "CO2")  # density
    mu = PropsSI("V", "T", t, "P", p, "CO2")  #
    di = 0.02  # internal diameter of single tube
    r = di / 2
    u = m / (rho * (np.pi * r ** 2))
    L = di
    Re = (rho * u * L) / mu     # from heat transfer book page 487:
    return Re



def get_enthalpy(p: Iterable[float], t: Iterable[float], fluid: str = "CO2") -> Iterable[float]:
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



def get_thermal_conductivity(p: Iterable[float], t: Iterable[float], fluid: str = "CO2") -> Iterable[float]:
    """Compute the thermal conductivity of a fluid at a given pressure and temperature.

    Args:
        p (Iterable[float]): Iterable containing the pressures of the fluid in Pa.
        t (Iterable[float]): Iterable containing the temperatures of the fluid in Kelvin.
        fluid (str, optional): Fluid to use. Defaults to "CO2".

    Returns:
        Iterable[float]: Iterable containing the thermal conductivity of the fluid at the given temperature and pressures in W/m/K.
    """
    thermal_conductivity = PropsSI("conductivity", "P", p, "T", t, fluid)
    return thermal_conductivity



def calculate_htc_a(t, p, m, tube: Tube):
    """Compute HTC for air for given temperature, pressure, mass flow rate and tube+fin design
    Args:
        pr_a (float): Prandtl number.
        re_a (float): Reynolds number.
        tube (class element): class describing the geometry of the tube
    Returns:
        HTC_air: heat transfer coefficient for air .
    """
    cross_section_air = tube.calculate_cross_section_air()
    pr_a = calculate_pr_air(t, p)
    re_a = calculate_re_air(t, p, m, cross_section_air, tube)
    k_a = 0.025  # W/(mK) thermal conductivity of air #TODO this varies over temperature and pressure
    htc_a = k_a/ tube.d_ot*(
        0.134
        * pr_a ** (1 / 3)
        * re_a ** 0.681
        * (2 * (tube.p_f - tube.t_f) / (tube.d_f - tube.d_r))**0.2
        * ((tube.p_f - tube.t_f) / tube.t_f) ** 0.1134
    )
    return htc_a


def calculate_htc_s(t, p, m):
    re_s = calculate_re_co2(t, p, m)
    pr_s = calculate_pr_co2(t, p)
    rho_pc = 800  # change / look up
    rho_s =  PropsSI("D", "T", t, "P", p, "CO2")  # [kg/m^3]
    
    t_pc = ( 273.15 +
        -122.6 + 6.124 * (p * 1e-5) - 0.1657 * (p * 1e-5) ** 2 + 0.01773 * (p * 1e-5) ** 2.5 - 0.0005608 * (p * 1e-5) ** 3 
    )
    # pressure in eq for t_pc has to be in bar
    # t_pc is calculated in Celsius, added +273.15 term to get Kelvin

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
    Nu_s = a * (re_s ** b) * (pr_s ** c) * (rho_pc / rho_s) ** n       # Nusselt number for sCO2 side
    k_s = get_thermal_conductivity(p,t,"CO2")                          # thermal conductivity of sCO2 at givem T, p
    htc_s = Nu_s * k_s / tube.d_i                                      # heat transfer coefficient
    
    return htc_s


### T and P are different for air and CO2, OHTC can not be calculated with one pair of values
### 
def compute_ohtc(t_air, t_s, p_air, p_s, m_air, m_s, tube):
    """ Computes overall heat transfer coefficient of the heat exchanger for one tube. 
    This is done by calculating HTC_air and HTC_sco2 and taken the inverse of their sum. 
    Thermal resistance of the tube wall r2 is neglected.

    Args:
        t_air (_type_): _description_
        t_s (_type_): _description_
        p_air (_type_): _description_
        p_s (_type_): _description_
        m_air (_type_): _description_
        m_s (_type_): _description_
        tube (_type_): _description_

    Returns:
        _type_: _description_
    """    
    surface_area_co2 = tube.calculate_surface_area_co2()
    [A_t,A_f] = tube.calculate_surface_area_air()   # areas of the segment
    htc_a = calculate_htc_a(t_air, p_air, m_air, tube)
    htc_s = calculate_htc_s(t_s, p_s, m_s)
    
    #calculation of the fin efficiency
    k_f = 236 # [W*m*K from reference 27] thermal conductivity fin TOCHECK
    m_efficiency = np.sqrt(2*htc_a/(k_f*tube.t_f))
    O_efficiency = (tube.d_f/tube.d_ot -1)*(1+0.35*np.log(tube.d_f/tube.d_ot))
    efficiency_fin = np.tanh(m_efficiency*O_efficiency*(tube.d_ot/2))/(m_efficiency*O_efficiency*(tube.d_ot/2))
    surface_area_air = A_t + efficiency_fin*A_f

    r_1 = 1 / (surface_area_air * htc_a)  # resistance of cold
    r_2 = 0  # negliged resistance of wall
    r_3 = 1 / (surface_area_co2 * htc_s)  # resistance of hot
    ohtc = 1 / (r_1 + r_2 + r_3)
    return ohtc


if __name__  == "__main__":
    
    
    
    # TESTING:    

    tube = Tube()
    m_s =  58.34625 / 120 #443 / (120*100)
    m_a = 6.7 /constants.n_segments

    htc_s = calculate_htc_s(320, 8e6, m_s)
    print('htc_s', htc_s)

    htc_a = calculate_htc_a(300, 1e5, m_a, tube)
    print('htc_a',htc_a)

    ohtc = compute_ohtc(300, 320, 1e5, 8e6, m_a, m_s , tube)
    print('ohtc', ohtc)

    cross_section_air = tube.calculate_cross_section_air()
    re_air = calculate_re_air(300, 1e5, m_a, cross_section_air, tube )
    print('re_air',re_air)

    re_CO2 = calculate_re_co2(320, 8e6, m_s)
    print('re_CO2',re_CO2)

### htc of air does not change between tube and segment if mass flow rate is adjusted
### ohtc now scales proportionally between n segments = 1 (ohtc = 335) and n_segments = 100 (ohtc=3.38)