import numpy as np
from CoolProp.CoolProp import PropsSI

class Tube:
    def __init__(self):
        self.d_ot = 0.025  # [m]          outside tube diameter
        self.L_t = 9  # [m]  ##TOCHECK     length of tube
        self.t_f = 5e-4  #  [m]           fin thickness
        self.p_f = 2.8e-3  # [m]        fin pitch (distance between fins)
        self.n_f = np.floor(
            self.L_f / (self.t_f + self.p_f)
        )  #         number of fins (per tube at the moment)
        self.d_f = 57.15e-3  #  [m]     fin outside diameter
        self.d_r = 0.028  # [m] fin inside diameter

    def calculate_cross_section_air(self):
        cross_section_tube_s = np.pi * self.d_ot * (
            self.L_t - self.t_f * self.n_f
        ) + np.pi(
            (self.d_f ** 2 - self.d_ot ** 2) / 2
        )  # [m^2]
        return cross_section_tube_s

    def calculate_surface_area_air(self):
        surface_area_tube = np.pi * self.d_ot * (
            self.L_t - self.t_f * self.n_f
        ) + np.pi * ((d_f ** 2 - d_ot ** 2) / 2)
        return surface_area_tube


def calculate_re_co2(t, p, m):
    rho = PropsSI("D", "T", t, "P", p, "CO2")  # density
    mu = PropsSI("V", "T", t, "P", p, "CO2")  #
    di = 0.02  # internal diameter of single tube
    r = di / 2
    u = m / (rho * (np.pi * r ** 2))
    L = di
    Re = (rho * u * L) / mu
    return Re


def calculate_pr_co2(t, p):
    Pr = PropsSI("PRANDTL", "T", t, "P", p, "CO2")
    return Pr


def calculate_pr_air(t, p):
    Pr = PropsSI("PRANDTL", "T", t, "P", p, "AIR")
    return Pr


# WORKING ON THIS
def calculate_re_air(t, p, m, cross_section_air, tube: Tube):
    rho = PropsSI("D", "T", t, "P", p, "AIR")  # density
    mu = PropsSI("V", "T", t, "P", p, "AIR")  #
    wetter_perimeter = (
        2 * tube.L_t
        + 2 * (tube.t_f - tube.d_ot)
        + (tube.d_f - tube.d_ot) * 2 * tube.n_f
    )
    d_f = (
        4 * cross_section_air / wetter_perimeter
    )  #  hydraulic diameter defined in page 582 pdf Fundamentals of...
    a_ff = cross_section_air
    re_a = m * d_f / (a_ff * mu)
    return re_a


def calculate_htc_a(t, p, m, tube: Tube):
    """Compute OHTC
    Args:
        pr_a (float): Prandtl number.
        re_a (float): Reynolds number.
        tube (class element): class describing the geometry of the tube
    Returns:
        OHTC: overall heat transfer coefficient for air .
    """
    cross_section_air = tube.calculate_cross_section_air()
    pr_a = calculate_pr_air(t, p)
    re_a = calculate_re_air(t, p, m, cross_section_air, tube)
    k_a = 0.025  # W/(mK) thermal conductivity of air #TODO this varies over temperature and pressure
    htc_a = k_a / tube.d_ot(
        0.134
        * pr_a ** (1 / 3)
        * re_a ** 0.681
        * (2 * (tube.p_f - tube.t_f) ** 0.2 / (tube.d_f - tube.d_r))
        * ((tube.p_f - tube.t_f) / t_f) ** 0.1134
    )
    return htc_a


def calculate_htc_s(t, p, m):
    re_s = calculate_re_co2(p, t, m)
    pr_s = calculate_pr_co2(p, t)
    rho_pc = None  # change
    rho_s = None  # change
    t_pc = (
        -122.6 + 6.124 * p - 0.1657 * p ** 2 + 0.01773 * p ** 2.5 - 0.0005608 * p ** 3
    )
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
    htc_s = a * (re_s ** b) * (pr_s ** c) * (rho_pc / rho_s) ** n
    return htc_s


def compute_ohtc(t, p, re_s, pr_s, tube):
    # TODO do you need to multiply by n tubes in the z direction?
    surface_area_co2 = tube.calculate_surface_area_co2()
    surface_area_air = tube.calculate_surface_area_air()  # area of the tube
    htc_a = calculate_htc_a(tube)
    htc_s = calculate_htc_s(t, p)
    r_1 = 1 / (surface_area_air * htc_a)  # resistance of cold
    r_2 = 0  # negliged resistance of wall
    r_3 = 1 / (surface_area_co2 * htc_s)  # resistance of hot
    ohtc = 1 / (r_1 + r_2 + r_3)
    pass