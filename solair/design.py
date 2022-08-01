from .constants import constants, get_m_air
from CoolProp.CoolProp import PropsSI
import numpy as np
from typing import Iterable


class Tube:
    """Class describing the parameters of tube using Eshan2018 reference.
    
         Attributes:
         ----------
            tube_out_diameter (float): Tube outside diameter [m]
            tube_in_diameter (float): Tube inside diameter [m]
            tube_length (float): Length of tube [m]
            tube_transverse_pitch (float): Transversal tube pitch (distance between the center of two tubes in the vertical direction) [m]
            fin_out_diameter (float): Fin outside diameter [m]
            fin_in_diameter (float): Fin inside diameter [m]  [0.01113,0.04089]
            fin_pitch (float): Fin pitch (distance between fins) [m] [0.0013,0.00406]
            fin_thickness (float): Fin root thickness [m]  [3.3e-4,2.02e-3]
            
            Addional information on bounds:
            - Bound on half fin diameter (d_f-d_r)/2 [m] [0.00142,0.01657]
            - Bound on Reynolds number for air [1000,18000]
        """
    def __init__(
        self,
        tube_out_diameter: float = 25e-3,
        tube_in_diameter: float = 20e-3,
        tube_segment_length: float = 0.2,
        n_segments: int = constants.n_segments,
        tube_transverse_pitch=58e-3,
        fin_out_diameter: float = 57e-3,
        fin_in_diameter: float = 28e-3,
        fin_pitch: float = 2.8e-3,
        fin_thickness: float = 7.5e-4,
    ):
        
        self.tube_out_diameter = tube_out_diameter 
        self.tube_in_diameter = tube_in_diameter 
        self.length = tube_segment_length * n_segments
        self.n_segments: int = n_segments
        self.segment_length = tube_segment_length # length of tube / number of segments = length of segment
        self.tube_transverse_pitch = tube_transverse_pitch   
        self.fin_pitch = fin_pitch  
        self.fin_thickness = fin_thickness   

        self.n_f = np.floor(
            self.segment_length / (self.fin_thickness + self.fin_pitch)
        )                  #         number of fins (per segment)
        self.fin_out_diameter = fin_out_diameter  
        self.fin_in_diameter = fin_in_diameter 
 
    def calculate_cross_section_air(self):
        """
        Calculates cross section of air as in equation 5 of ref2 (TODO check this)
        
        Returns
        -------
        cross section of the tube segment
        """
        # TODO replace tube_out_diameter with fin_in_diameter ?
        # cross_section_tube_s = (self.tube_transverse_pitch-self.tube_out_diameter)*self.segment_length-(self.fin_out_diameter-self.tube_out_diameter)*self.fin_thickness*self.n_f
        # Now this calculates the area between tubes where air can flow
        cross_section_tube_s = (self.tube_transverse_pitch-self.fin_in_diameter)*self.segment_length-(self.fin_out_diameter-self.fin_in_diameter)*self.fin_thickness*self.n_f * 2
        assert cross_section_tube_s > 0, "Cross section of tube segment is negative"
        return cross_section_tube_s

    def calculate_surface_area_air(self):
        """
        Calculates surface area of air as in equation 6 of ref2
        
        Returns
        -------
        A_t : outside tube area between fins for entire tube
        A_f : surface area of all fins along one tube
        """# surface area of entire tube
        A_1f = np.pi * ((self.fin_out_diameter ** 2 - self.tube_out_diameter ** 2) / 2) # surface area of one fin
        A_f = self.n_f * A_1f   
        A_t = np.pi * self.tube_out_diameter * (
            self.segment_length - self.fin_thickness * self.n_f) 
        return [A_t,A_f]

    def calculate_surface_area_co2(self):
        surface_area_co2 = np.pi*self.tube_in_diameter*self.segment_length          
        return surface_area_co2


def calculate_pr_co2(t, p):
    Pr = PropsSI("PRANDTL", "T", t, "P", p, "CO2")
    return Pr


def calculate_pr_air(t, p):
    Pr = PropsSI("PRANDTL", "T", t, "P", p, "AIR")
    return Pr


def calculate_re_air(t, p, m, cross_section_air, tube: Tube):
    mu = PropsSI("V", "T", t, "P", p, "AIR")  # viscosity
    # Perimeter around the area the air can flow through
    wetted_perimeter = (
        2 * tube.segment_length
        + 2 * (tube.tube_transverse_pitch - tube.fin_in_diameter)
        + (tube.fin_out_diameter - tube.fin_in_diameter) * 2 * 2 * tube.n_f # two sides per fin and two for fin on either side.
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
    di = 20e-3  # internal diameter of single tube
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
    htc_a = k_a/ tube.tube_out_diameter*(
        0.134
        * pr_a ** (1 / 3)
        * re_a ** 0.681
        * (2 * (tube.fin_pitch - tube.fin_thickness) / (tube.fin_out_diameter - tube.fin_in_diameter))**0.2
        * ((tube.fin_pitch - tube.fin_thickness) / tube.fin_thickness) ** 0.1134
    )
    return htc_a


def calculate_htc_s(t, p, m,tube):
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
def compute_ohtc(t_air, t_s, p_air, p_s, m_air, m_s, tube: Tube):
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
    htc_s = calculate_htc_s(t_s, p_s, m_s,tube)
    
    #calculation of the fin efficiency
    k_f = 204 #[W*m*K]       #done (adjust to Ehsan)  thermal conductivity fin 
    m_efficiency = np.sqrt(2*htc_a/(k_f*tube.fin_thickness))
    O_efficiency = (tube.fin_out_diameter/tube.tube_out_diameter -1)*(1+0.35*np.log(tube.fin_out_diameter/tube.tube_out_diameter))
    efficiency_fin = np.tanh(m_efficiency*O_efficiency*(tube.tube_out_diameter/2))/(m_efficiency*O_efficiency*(tube.tube_out_diameter/2))
    surface_area_air = A_t + efficiency_fin*A_f

    r_1 = 1 / (surface_area_air * htc_a)  # resistance of cold
    r_2 = 0  # negliged resistance of wall
    r_3 = 1 / (surface_area_co2 * htc_s)  # resistance of hot
    ohtc = 1 / (r_1 + r_2 + r_3)
    if np.isnan(ohtc):
        raise ValueError("ohtc is nan. Check design of the tube.")
    return ohtc


if __name__  == "__main__":
    
    
    
    # TESTING:    

    tube = Tube()
    m_s =  406.6 / 190 / 49 # #443 / (120*100)
    #m_a = 1107 /constants.n_segments
    m_a = get_m_air(
        constants.p_co2_inlet,
        constants.p_co2_outlet,
        constants.t_co2_inlet,
        constants.t_co2_outlet,
        constants.t_air_inlet,
        constants.t_air_outlet,
        constants.m_co2,
        constants.cp_air
    )
    print('m_air = ', round(m_a,2), 'ref: 1107')

    m_a = m_a /47.5 /49  # for htc_air we are looking at air flow across one tube

    htc_s = calculate_htc_s(constants.t_co2_inlet, constants.p_co2_inlet, m_s, tube)
    print('htc_s', round(htc_s,2), 'ref: ca 750 (fig11)')

    htc_a = calculate_htc_a(constants.t_air_inlet, constants.p_air_in, m_a, tube)
    print('htc_a',round(htc_a,2), 'ref: 31.3 avg in first row')

    ohtc = compute_ohtc(300, 320, 1e5, 8e6, m_a, m_s , tube)
    print('ohtc', ohtc)

    cross_section_air = tube.calculate_cross_section_air()
    re_air = calculate_re_air(300, 1e5, m_a, cross_section_air, tube )
    print('re_air',re_air)

    re_CO2 = calculate_re_co2(320, 8e6, m_s)
    print('re_CO2',re_CO2)

    area_air = tube.calculate_surface_area_air()
    print('area_air', (area_air[0]+area_air[1])*190*49*30)

    area_CO2 = tube.calculate_surface_area_co2()
    print('area_co2', area_CO2*190*49*30)



    '''
    enthalpy = get_enthalpy(836, 320)
    print('enthalpy', enthalpy)

    print('m_air',constants.m_air)


    #c_p_air = 

    C_min = 1005 * 6.7 #10.21   # 10291
    C_max = 525e3 * 48 #3.7    # 1.944e3
    print ('C_max', C_max)

    UA =  4.4 * 100 * 360              #  = OHTC = OHTC_segmant * n_segments * n_tubes

    C = C_min / C_max # =  0.005293
    N = UA / C_min #  0.0004276 for one segment

    e = (1 - np.exp(-C*(1-np.exp(-N))))/C   # cross flow with Cmax fluid mixed

    e2 = 1 - np.exp( N**(0.22) * (np.exp(-C*N**(0.78))-1)/C)
    
    print('N', N)
    print('effectivity of our heat exchanger:',e2)
'''
### htc of air does not change between tube and segment if mass flow rate is adjusted
### ohtc now scales proportionally between n segments = 1 (ohtc = 335) and n_segments = 100 (ohtc=3.38)