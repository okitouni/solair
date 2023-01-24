from .design import Tube
#from .constants import constants
import numpy as np
from .design import calculate_pressure_drop_air
from .utils import drop_pressure


def calculate_tube_cost_per_length(rho_tube, c_tm, tube: Tube):
    """Compute cost of tube core per unit length [$/m]
    Args:
        rho_tube (float):   density of tube material [kg/m^3]
        c_tm (float):       cost of tube material per unit mass [$/kg]
        tube (class element): class describing the geometry of the tube
    Returns:
        cost_tube: cost of tube per unit length
    """

    cost_tube = (
        np.pi
        * rho_tube
        * (tube.tube_out_diameter ** 2 - tube.tube_in_diameter ** 2)
        * c_tm
        / 4
    )
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

    cost_fin = (
        np.pi
        * rho_fin
        / (4 * tube.fin_pitch)
        * (
            (tube.fin_out_diameter ** 2 - tube.fin_in_diameter ** 2)
            * tube.fin_thickness
            + (tube.fin_in_diameter ** 2 - tube.tube_out_diameter ** 2)
            * (tube.fin_pitch - tube.fin_thickness)
        )
        * c_fm
    )
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

    total_cost_finned_tubes = (
        (
            tube.constants_t.weighting_factor * (cost_tube + cost_fin)
            + tube.constants_t.fixed_cost_tube
        )
        * (tube.segment_length * tube.n_segments)
        * tube.constants_t.n_tubes_in_row
        * tube.constants_t.n_rows
        * tube.constants_t.n_sub_heat_exchangers
    )  # TODO: implement number heat exchangers
    return total_cost_finned_tubes


def calculate_cost_air_cooler_no_fans(rho_tube, c_tm, rho_fin, c_fm, tube: Tube):
    """Compute total cost of the finned tubes [$/m]
    Args:
        rho_tube (float):   density of tube material [kg/m^3]
        c_tm (float):       cost of tube material per unit mass [$/kg]
        rho_fin (float):   density of fin material [kg/m^3]
        c_fm (float):       cost of fin material per unit mass [$/kg]
        pressure_drop (float):      air pressure drop across heat exchanger (Pa)
        lifetime_years (float):     assumed lifetime for cost calculation (years), standard: 25
        LCOE_fanpower_cents (float):    assumed cost of kWh elecricity consumed by fans (cents), standard: 5
        tube (class element): class describing the geometry of the tube
    Returns:
        total_cost_finned_tube: total cost of finned tube
    """
    # cost_tube = calculate_tube_cost_per_length(rho_tube, c_tm, tube)
    # cost_fin = calculate_fin_cost_per_length(rho_fin, c_fm, tube)
    total_cost_finned_tubes = calculate_total_cost_finned_tubes(
        rho_tube, c_tm, rho_fin, c_fm, tube
    )

    cost_air_cooler_no_fans = (
        total_cost_finned_tubes
        * (1 + constants.f_header)
        * (1 + constants.f_labor)
        * constants.f_HX
    )
    # this is assuming fan cost is not increasing header, labor and HX costs
    return cost_air_cooler_no_fans


# fan: HC-100 6/H: Best average between cost and power
class Fan_HC_100_6H:
    def __init__(self,):
        self.airflow = 37000  # [m^3/h]
        self.cost = 1255.55  # [Euro]
        self.power = 1500  # [Watt]
        self.fan_speed = 955  # [rpm]

    def get_airflow_for_pressure(self, pressure_drop: float):
        # pressure drop in Pascal
        if pressure_drop < 320:
            airflow_for_pressure = 36e3 - 26e3 * (
                pressure_drop / 248.84
            )  # conversion Pa to inch H20
        else:
            airflow_for_pressure = 1
        return airflow_for_pressure
        # validation DONE


# fan: HFW-90 4/T-7.5: Optimal for airflow/cost
class Fan_HFW_90_4_T_75:
    def __init__(self,):
        self.airflow = 46150  # [m^3/h]
        self.cost = 1259.05  # [Euro]
        self.power = 5500  # [Watt]
        self.fan_speed = 1440  # [rpm]

    def get_airflow_for_pressure(self, pressure_drop: float):
        # pressure drop in Pascal
        if pressure_drop < 300:
            airflow_for_pressure = 46e3 - 8.5e3 * (
                pressure_drop / 248.84
            )  # conversion Pa to inch H20
        elif pressure_drop < 700:
            airflow_for_pressure = 37.5e3 - 75 * (
                pressure_drop - 250
            )  # conversion Pa to inch H2
        else:
            airflow_for_pressure = 1
        return airflow_for_pressure
        # validation DONE


# fan: CJHCH-80-8T-0.5 : Optimal for airflow/power
class Fan_CJHCH_80_8T_05:
    def __init__(self,):
        self.airflow = 16600  # [m^3/h]
        self.cost = 1379.25  # [Euro]
        self.power = 370  # [Watt]
        self.fan_speed = 700  # [rpm]

    def get_airflow_for_pressure(self, pressure_drop: float):
        # pressure drop in Pascal
        # Need to find fan curves for this pressure
        if pressure_drop < 100:
            airflow_for_pressure = 17e3 - 23.6e3 * (
                pressure_drop / 248.84
            )  # conversion Pa to inch H20
        else:
            airflow_for_pressure = 1

        return airflow_for_pressure
        # validation DONE


# fan: CJHCH-100-8T-1.5 :
class Fan_CJHCH_100_8T_15:
    def __init__(self,):
        self.airflow = 32500  # [m^3/h]
        self.cost = 1850  # [Euro]
        self.power = 1100  # [Watt]
        self.fan_speed = 850  # [rpm]

    def get_airflow_for_pressure(self, pressure_drop: float):
        # pressure drop in Pascal
        # Need to find fan curves for this pressure
        if pressure_drop < 175:
            airflow_for_pressure = 4e3 + 2116 * np.sqrt(
                175 - pressure_drop
            )  # conversion Pa to inch H20
        else:
            airflow_for_pressure = 1

        return airflow_for_pressure
        # validation DONE


# fan: HC-63-6T/H : Optimal for airflow//cost/power
class Fan_HC_63_6T_H:
    """Class describing the parameters of fan HC-100 6/H
    
         Attributes:
         ----------

    """

    def __init__(self,):
        self.airflow = 12350  # [m^3/h]
        self.cost = 535.85  # [Euro]
        self.power = 370  # [Watt]
        self.fan_speed = 900  # [rpm]

    def get_airflow_for_pressure(self, pressure_drop: float):
        # pressure drop in Pascal
        if pressure_drop < 220:
            airflow_for_pressure = 12e3 - 13e3 * (
                pressure_drop / 248.84
            )  # conversion Pa to inch H20
        else:
            airflow_for_pressure = 1

        return airflow_for_pressure
        # validation DONE


# def fan_lifetime_cost(pressure_drop: float, lifetime_years: float, LCOE_fanpower_cents: float, fan1: Fan_CJHCH_100_8T_15, fan2: Fan_CJHCH_80_8T_05, fan3: Fan_HC_100_6H, fan4: Fan_HFW_90_4_T_75, fan5: Fan_HC_63_6T_H):


def calculate_fan_lifetime_cost(
    pressure_drop: float, lifetime_years: float, LCOE_fanpower_cents: float, tube: Tube
):
    # 1

    fan_CJHCH_100_8T_15 = Fan_CJHCH_100_8T_15()
    airflow_Fan_CJHCH_100_8T_15 = fan_CJHCH_100_8T_15.get_airflow_for_pressure(
        pressure_drop
    )
    n_fans_required_CJHCH_100_8T_15 = (
        tube.constants_t.m_air * 3000 / airflow_Fan_CJHCH_100_8T_15
    )  # *3000 for conversion kg/s to m^3/h
    fan_lifetime_cost_CJHCH_100_8T_15 = n_fans_required_CJHCH_100_8T_15 * (
        fan_CJHCH_100_8T_15.cost
        + fan_CJHCH_100_8T_15.power
        / 1000  # /1000 for conversion power in W to kW
        * LCOE_fanpower_cents
        / 100
        * 24
        * 365
        * lifetime_years
    )  # /100 for conversion cost in cents to $
    # 2
    fan_HC_100_6H = Fan_HC_100_6H()
    airflow_Fan_HC_100_6H = fan_HC_100_6H.get_airflow_for_pressure(pressure_drop)
    n_fans_required_HC_100_6H = (
        tube.constants_t.m_air * 3000 / airflow_Fan_HC_100_6H
    )  # *3000 for conversion kg/s to m^3/h
    fan_lifetime_cost_HC_100_6H = n_fans_required_HC_100_6H * (
        fan_HC_100_6H.cost
        + fan_HC_100_6H.power
        / 1000  # /1000 for conversion power in W to kW
        * LCOE_fanpower_cents
        / 100
        * 24
        * 365
        * lifetime_years
    )  # /100 for conversion cost in cents to $
    # 3
    fan_HFW_90_4_T_75 = Fan_HFW_90_4_T_75()
    airflow_Fan_HFW_90_4_T_75 = fan_HFW_90_4_T_75.get_airflow_for_pressure(
        pressure_drop
    )
    n_fans_required_HFW_90_4_T_75 = (
        tube.constants_t.m_air * 3000 / airflow_Fan_HFW_90_4_T_75
    )  # *3000 for conversion kg/s to m^3/h
    fan_lifetime_cost_HFW_90_4_T_75 = n_fans_required_HFW_90_4_T_75 * (
        fan_HFW_90_4_T_75.cost
        + fan_HFW_90_4_T_75.power
        / 1000  # /1000 for conversion power in W to kW
        * LCOE_fanpower_cents
        / 100
        * 24
        * 365
        * lifetime_years
    )  # /100 for conversion cost in cents to $
    # 4
    fan_CJHCH_80_8T_05 = Fan_CJHCH_80_8T_05()
    airflow_Fan_CJHCH_80_8T_05 = fan_CJHCH_80_8T_05.get_airflow_for_pressure(
        pressure_drop
    )
    n_fans_required_CJHCH_80_8T_05 = (
        tube.constants_t.m_air * 3000 / airflow_Fan_CJHCH_80_8T_05
    )  # *3000 for conversion kg/s to m^3/h
    fan_lifetime_cost_CJHCH_80_8T_05 = n_fans_required_CJHCH_80_8T_05 * (
        fan_CJHCH_80_8T_05.cost
        + fan_CJHCH_80_8T_05.power
        / 1000  # /1000 for conversion power in W to kW
        * LCOE_fanpower_cents
        / 100
        * 24
        * 365
        * lifetime_years
    )  # /100 for conversion cost in cents to $
    # 5
    fan_HC_63_6T_H = Fan_HC_63_6T_H()
    airflow_Fan_HC_63_6T_H = fan_HC_63_6T_H.get_airflow_for_pressure(pressure_drop)
    n_fans_required_HC_63_6T_H = (
        tube.constants_t.m_air * 3000 / airflow_Fan_HC_63_6T_H
    )  # *3000 for conversion kg/s to m^3/h
    fan_lifetime_cost_HC_63_6T_H = n_fans_required_HC_63_6T_H * (
        fan_HC_63_6T_H.cost
        + fan_HC_63_6T_H.power
        / 1000  # /1000 for conversion power in W to kW
        * LCOE_fanpower_cents
        / 100
        * 24
        * 365
        * lifetime_years
    )  # /100 for conversion cost in cents to $

    cost_list = [
        fan_lifetime_cost_CJHCH_100_8T_15,
        fan_lifetime_cost_HC_100_6H,
        fan_lifetime_cost_HFW_90_4_T_75,
        fan_lifetime_cost_CJHCH_80_8T_05,
        fan_lifetime_cost_HC_63_6T_H,
    ]
    cost_list_round = [round(i, 1) for i in cost_list]
    minimum_cost = min(cost_list)
    # print('cost list for 5 fan types: ' , cost_list_round)
    return minimum_cost
    # return cost_list


def calculate_total_cost_air_cooler(tube_cost, fan_cost, tube):
    """Compute total cost of the finned tubes [$/m]
    Args:
        tube_cost (float): cost of the finned tube
        fan_cost (float): cost of the fan 
    Returns:
        total_cost_finned_tube: total cost of finned tube
    """

    total_cost_air_cooler = (
        tube_cost
        * (1 + tube.constants_t.f_header)
        * (1 + tube.constants_t.f_labor)
        * tube.constants_t.f_HX
        + fan_cost
    )
    # this is assuming fan cost is not increasing header, labor and HX costs
    return total_cost_air_cooler

def calculate_sub_cost_air_cooler(
    rho_tube, c_tm, rho_fin, c_fm, lifetime_years, LCOE_fanpower_cents, tube: Tube
):
    """Compute sub costs of the finned tubes [$/m]
    Args:
        rho_tube (float):   density of tube material [kg/m^3]
        c_tm (float):       cost of tube material per unit mass [$/kg]
        rho_fin (float):   density of fin material [kg/m^3]
        c_fm (float):       cost of fin material per unit mass [$/kg]
        pressure_drop (float):      air pressure drop across heat exchanger (Pa)
        lifetime_years (float):     assumed lifetime for cost calculation (years), standard: 25
        LCOE_fanpower_cents (float):    assumed cost of kWh elecricity consumed by fans (cents), standard: 5
        tube (class element): class describing the geometry of the tube
    Returns:
        total_cost_finned_tube: total cost of finned tube
    """

    # cost_tube = calculate_tube_cost_per_length(rho_tube, c_tm, tube)
    # cost_fin = calculate_fin_cost_per_length(rho_fin, c_fm, tube)
    total_cost_finned_tubes = calculate_total_cost_finned_tubes(
        rho_tube, c_tm, rho_fin, c_fm, tube
    )

    pressure_drop = calculate_pressure_drop_air(tube)
    fan_lifetime_cost = calculate_fan_lifetime_cost(
        pressure_drop, lifetime_years, LCOE_fanpower_cents, tube,
    )

    # this is assuming fan cost is not increasing header, labor and HX costs
    return [total_cost_finned_tubes, fan_lifetime_cost]

### PLOTTING:       need to replace line 'return minimum_cost' with return cost_list
"""
fan1 = []
fan2 = []
fan3 = []
fan4 = []
fan5 = []

for i in range(0, 170, 1):
    run = calculate_fan_lifetime_cost(i, 25, 5)
    fan1.append(run[0])
    fan2.append(run[1])
    fan3.append(run[2])
    fan4.append(run[3])
    fan5.append(run[4])

x = np.arange(0, 170, 1)
plt.plot(x, fan1,  label='fan1')
plt.plot(x, fan2,  label='fan2')
plt.plot(x, fan3, label='fan3')
plt.plot(x, fan4, label='fan4')
plt.plot(x, fan5, label='fan5')
plt.xlabel('Pressure drop across heat exchanger (Pa)')
plt.ylabel('Lifetime cost of fans required ($)')
plt.legend()
plt.show()
"""


"""   
    1. air flow             is proportional to      fan speed 
    2. static pressure      is proportional to     (fan speed)^2
    3. power consumption    is proportional to     (fan speed)^3

procedure:  a) look at required Static pressure (= air pressure drop)
            b) find required fan speed 
            c) use fan performance graphs to find airflow and power consumption
            d) divide m_air_total / m_air_fan to find number of fans and total power consumed by fans


        v)  write # occurences of different temperatures across year into a file and run fan cost calc for varying temperatures
            (get m_air for different ambient temp, switch off some fans --> adapt running cost of fans)

        vi) how do we consider 40+ C temperatures and not being able to cool CO2 enough...loosing efficiency of power cycle, but how consider it?
"""

### Testing: references are for Khatoon et al paper with 10MW plant
"""
tube = Tube()

cost_t = calculate_tube_cost_per_length(constants.rho_steel, constants.cost_steel, tube)
print(cost_t, 'ref: 1.12 $/m')

cost_f = calculate_fin_cost_per_length(constants.rho_alu, constants.cost_alu, tube)
print(cost_f, 'ref: 5.21 $/m')

cost_total_finned_tubes = calculate_total_cost_finned_tubes(constants.rho_steel, constants.cost_steel,constants.rho_alu, constants.cost_alu, tube)
print(cost_total_finned_tubes, 'ref: 47494.72$')

cost_air_cooler_no_fans = calculate_cost_air_cooler_no_fans(constants.rho_steel, constants.cost_steel,constants.rho_alu, constants.cost_alu,  tube)
print('cost air cooler no fans', cost_air_cooler_no_fans)
total_cost_air_cooler = calculate_total_cost_air_cooler(constants.rho_steel, constants.cost_steel,constants.rho_alu, constants.cost_alu, 25, 5,  tube)
print('total cost air cooler', total_cost_air_cooler) #, 'ref: 174400.61$')
"""

