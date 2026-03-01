import math
import numpy as np


def saturation_vapor_pressure(T):
    """
    Saturation vapor pressure over liquid water (Pa) via Clausius–Clapeyron
    with constant latent heat.

    Parameters
    ----------
    T : float or array-like
        Temperature in Kelvin.

    Returns
    -------
    float or np.ndarray
        Saturation vapor pressure in Pa (same shape as T).
    """
    T = np.asarray(T, dtype=float)
    e0 = 611.2     # Pa at T0
    t0 = 273.15    # K
    rv = 461.5     # J kg^-1 K^-1
    lv = 2.501e6   # J kg^-1  (latent heat of vaporization, constant)
    es = e0 * np.exp((lv / rv) * (1.0 / t0 - 1.0 / T))
    return es if es.shape else float(es)



def vapor_concentration_RH_Tetens(T_c, RH=100.0):
    T_c = np.asarray(T_c, dtype=float)
    RH  = np.clip(np.asarray(RH, dtype=float), 0.0, 100.0)

    # Tetens (T in °C): e_s in kPa → ×1000 to Pa
    e_s = 0.6108 * np.exp((17.27 * T_c) / (T_c + 237.3)) * 1000.0  # Pa
    e   = (RH / 100.0) * e_s

    M_w = 18.015  # g mol^-1
    R   = 8.314   # J mol^-1 K^-1
    C = (e * M_w) / (R * (T_c + 273.15))  # g m^-3

    return C if C.shape else float(C)



def vapor_concentration_RH(T_c, RH=100.0):
    """
    Absolute humidity (g m^-3) at temperature T_c (°C) and relative humidity RH (%).
    Uses saturation_vapor_pressure(T_K) and ideal gas law.
    """
    T_c = np.asarray(T_c, dtype=float)

    # saturation vapor pressure (Pa) using CC function (expects Kelvin)
    e_s = saturation_vapor_pressure(T_c + 273.15)  # Pa

    # partial pressure of water vapor (Pa)
    e = (RH / 100.0) * e_s

    # absolute humidity C = rho_v (g m^-3) using ideal gas:  rho = (e * M_w) / (R * T)
    M_w = 18.015    # g mol^-1
    R   = 8.314     # J mol^-1 K^-1
    C = (e * M_w) / (R * (T_c + 273.15))   # g m^-3

    return C if C.shape else float(C)



def vapor_concentration(es: float, T: float) -> float:
    M = 18.015
    R = 8.314
    rho = (es * M) / (R * T)
    return float(rho)


def e_sat(T: float) -> float:
    e_s = 0.6108 * math.exp((17.27 * T) / (T + 237.3)) * 1000
    return float(e_s)
