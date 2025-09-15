import math
import numpy as np


def saturation_vapor_pressure(T: float) -> float:
    e0 = 611.2
    T0 = 273.15
    Lv = 2.501e6
    Rv = 461.5
    es = e0 * np.exp((Lv / Rv) * (1 / T0 - 1 / T))
    return float(es)


def vapor_concentration_RH(T: float, RH: float = 100.0) -> float:
    M_w = 18.015
    R = 8.314
    e_s = 0.6108 * math.exp((17.27 * T) / (T + 237.3)) * 1000
    e = RH / 100 * e_s
    C = (e * M_w) / (R * (T + 273.15))
    return float(C)


def vapor_concentration(es: float, T: float) -> float:
    M = 18.015
    R = 8.314
    rho = (es * M) / (R * T)
    return float(rho)


def e_sat(T: float) -> float:
    e_s = 0.6108 * math.exp((17.27 * T) / (T + 237.3)) * 1000
    return float(e_s)
