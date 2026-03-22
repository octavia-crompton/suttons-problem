"""
Physical constants used throughout the sutton package.

Import from here rather than defining ad-hoc module-level globals.
"""

# Stefan–Boltzmann constant (W m^-2 K^-4)
SIGMA_SB: float = 5.670374419e-8

# Dry-air properties (standard conditions, ~20 °C, 1 atm)
CP_AIR: float = 1005.0    # J kg^-1 K^-1  (specific heat at constant pressure)
RHO_AIR: float = 1.20     # kg m^-3

# Latent heat of vaporisation (~30 °C)
LV: float = 2.43e6        # J kg^-1
LV_G: float = 2430.0      # J g^-1   (= LV / 1000)

# Water-vapour gas constant
RV: float = 461.5          # J kg^-1 K^-1
