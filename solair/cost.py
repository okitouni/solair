from .design import Tube
from .constants import constants
import numpy as np


def calculate_tube_cost_per_length(rho_tube, c_tm, tube: Tube):
    """Compute cost of tube core per unit length [$/m]
    Args:
        rho_tube (float):   density of tube material [kg/m^3]
        c_tm (float):       cost of tube material per unit mass [$/kg]
        tube (class element): class describing the geometry of the tube
    Returns:
        cost_tube: cost of tube per unit length
    """

    cost_tube = np.pi * rho_tube * (tube.d_ot**2 - tube.d_i**2) * c_tm / 4
    return cost_tube

def calculate_fin_cost_per_length(rho_fin, c_fm, tube: Tube):
    """Compute cost of fins per unit tube length [$/m]
    Args:
        rho_fin (float):   density of fin material [kg/m^3]
        c_fm (float):       cost of fin material per unit mass [$/kg]
        tube (class element): class describing the geometry of the tube
    Returns:
        cost_fin: cost of fins per unit tube length
    """

    cost_fin = np.pi * rho_fin / (4 * tube.p_f )  * ((tube.d_f**2 - tube.d_r**2) * tube.t_f + (tube.d_r**2 - tube.d_ot**2)*(tube.p_f - tube.t_f) ) * c_fm
    return cost_fin


def calculate_total_cost_finned_tubes(rho_tube, c_tm, rho_fin, c_fm, tube: Tube):
    """Compute total cost of the finned tubes [$/m]
    Args:
        rho_tube (float):   density of tube material [kg/m^3]
        c_tm (float):       cost of tube material per unit mass [$/kg]
        rho_fin (float):   density of fin material [kg/m^3]
        c_fm (float):       cost of fin material per unit mass [$/kg]
        tube (class element): class describing the geometry of the tube
    Returns:
        total_cost_finned_tube: total cost of finned tube
    """
    
    cost_tube = calculate_tube_cost_per_length(rho_tube, c_tm, tube)
    cost_fin = calculate_fin_cost_per_length(rho_fin, c_fm, tube)

    total_cost_finned_tubes = ( constants.weighting_factor * (cost_tube + cost_fin) + constants.fixed_cost_tube ) * ( tube.segment_length * tube.n_segments ) * constants.n_tubes_in_row * constants.n_rows * 3 # TODO: implement number heat exchangers
    return total_cost_finned_tubes



def calculate_total_cost_air_cooler(rho_tube, c_tm, rho_fin, c_fm, tube: Tube):
    """Compute total cost of the finned tubes [$/m]
    Args:
        rho_tube (float):   density of tube material [kg/m^3]
        c_tm (float):       cost of tube material per unit mass [$/kg]
        rho_fin (float):   density of fin material [kg/m^3]
        c_fm (float):       cost of fin material per unit mass [$/kg]
        tube (class element): class describing the geometry of the tube
    Returns:
        total_cost_finned_tube: total cost of finned tube
    """
    
    # cost_tube = calculate_tube_cost_per_length(rho_tube, c_tm, tube)
    # cost_fin = calculate_fin_cost_per_length(rho_fin, c_fm, tube)
    total_cost_finned_tubes = calculate_total_cost_finned_tubes(rho_tube, c_tm, rho_fin, c_fm, tube)

    total_cost_air_cooler =  total_cost_finned_tubes * ( 1 + constants.f_header ) * (1 + constants.f_labor ) * constants.f_HX 

    ### add cost of fans!

    
    return total_cost_air_cooler



### Testing: references are for Khatoon et al paper with 10MW plant
'''
tube = Tube()
cost_t = calculate_tube_cost_per_length(constants.rho_steel, constants.cost_steel, tube)
print(cost_t, 'ref: 1.12 $/m')

cost_f = calculate_fin_cost_per_length(constants.rho_alu, constants.cost_alu, tube)
print(cost_f, 'ref: 5.21 $/m')

cost_total_finned_tubes = calculate_total_cost_finned_tubes(constants.rho_steel, constants.cost_steel,constants.rho_alu, constants.cost_alu, tube)
print(cost_total_finned_tubes, 'ref: 47494.72$')

total_cost_air_cooler = calculate_total_cost_air_cooler(constants.rho_steel, constants.cost_steel,constants.rho_alu, constants.cost_alu, tube)
print(total_cost_air_cooler, 'ref: 174400.61$')
'''