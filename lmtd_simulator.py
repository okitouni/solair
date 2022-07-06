import numpy as np
from typing import Iterable

# Designing for 100 MW


def get_enthalpy(p, t):
    # this is a placeholder until we parametrize REFPROP
    return 800 * np.cos(2 * np.pi * p / t)


class constants:
    T_CO2_in = 89
    T_CO2_out = 35
    delta_t_CO2 = T_CO2_in - T_CO2_out  # 89 is inlet CO2, 35 is outlet CO2

    T_air_in = 25

    # calculate air mass flow rate
    # WRITE m_CO2 calculation
    m_CO2 = 950  # [kg/s] for 100MW
    Q_CO2_tot = m_CO2 * delta_t_CO2
    cp_air = 1

    enthalpy_table = {
        (p, t): get_enthalpy(p, t)
        for p in np.arange(1, 100, 10)
        for t in np.arange(33, 89)
    }

    OHTC = 1


def compute_htc(*params) -> float:
    pass


def lmtd(th: Iterable[float], tc: Iterable[float]) -> float:
    return ((th[0] - tc[1]) - (th[1] - tc[0])) / np.log(
        (th[0] - tc[1]) / (th[1] - tc[0])
    )


def energy_transferred(htc: float, th: Iterable[float], tc: Iterable[float]):
    return htc * lmtd(th, tc)


def energy_co2(p1, p0, t1, t0):
    h0 = get_enthalpy(p0, t0)
    h1 = get_enthalpy(p1, t1)
    return constants.m_CO2 * (h1 - h0)


def energy_air(mdot, t1, t0):
    return constants.cp_air * mdot * (t1 - t0)

def p_out(p_in):
    rho = 1
    f = 1 # from Eq. 16
    u = 1
    L = 1
    d = 1
    delta_p = (rho*f*u**2*L)/(2*d)
    p_out = p_in - delta_p
    return p_out

# calculate total Q exchanged
# NEED ENTHALPY OF CO2
h_CO2_in = constants.enthalpy_table[(1, 33)]
h_CO2_out = constants.enthalpy_table[(1, 34)]
Q = constants.m_CO2 * (h_CO2_out - h_CO2_in)

# find air mass flow rate
T_air_out = (constants.T_CO2_in + constants.T_CO2_out) / 2  # guess outlet air T
delta_t_air = T_air_out - constants.T_air_in
m_air = Q / (constants.cp_air * delta_t_air)


##################################

P_11 = 8 # from previous step
P_21 = p_out(P_11) 
Th_11 = 89 # known from cycle

# guess Th_i_j+1
Th_21 = 89 # make better guess

################ LOOP

# calculate Q1

h_11 = constants.enthalpy_table[(P_11, Th_11)]
h_21 = constants.enthalpy_table[(P_21, Th_21)]
Q1_11 = constants.m_CO2*(h_11 - h_21) # find h from tables

Tc_11 = constants.T_air_in # change if inside loop

Tc_12 = Q1_11 /(m_air*constants.cp_air) + Tc_11

# calculate delta_T_m
th = [Th_11,Th_21]
tc = [Tc_11,Tc_12]
delta_T_m = lmtd(th,tc)

Q2_11 = constants.OHTC*(delta_T_m)

# compare Q1 Q2

# update Th_1_2

################# 


# express Q1 -Q2 as fct of guessed temperature th_i+1_j

# Q1-Q2 = (OHTC * ( ) ) / () 