import math
import numpy as np

def thomas(aa, bb, cc, dd):
    """
    """
    n = len(bb)
    bet = np.zeros(n)
    gam = np.zeros(n)
    q = np.zeros(n)

    # Forward sweep
    bet[0] = bb[0]
    gam[0] = dd[0] / bb[0]

    for i in range(1, n):
        bet[i] = bb[i] - (aa[i] * cc[i-1] / bet[i-1])
        gam[i] = (dd[i] - aa[i] * gam[i-1]) / bet[i]

    # Back substitution
    q[n-1] = gam[n-1]

    for i in range(n-2, -1, -1):
        q[i] = gam[i] - (cc[i] * q[i+1] / bet[i])

    return q

def thomas_alt(aa, bb, cc, dd):
    """
    Solves for a vector q of length n 
    """
    n = len(bb)
    gam = np.zeros(n)
    q = np.zeros(n)

    # Forward sweep
    bet = bb[0]
    q[0] = dd[0]/bet


    for i in range(1, n):
        gam[i] = cc[i-1]/bet
        bet = bb[i]-aa[i]*gam[i]
        if bet ==0:
            break
        q[i] = (dd[i]-aa[i]*q[i-1])/bet

    for i in range(n-2, -1, -1):
        q[i] = q[i] - gam[i+1]*q[i+1]

    return q

def our_central_difference(s, dz):
    """
    Computes gradients of a data vector using central differencing. At edges, 
    forward and  backward differences are used
    """
    m = len(s)
    dsdz = np.zeros(m)
    
    # Central difference for interior points
    dsdz[1:m-1] = (s[2:m] - s[0:m-2]) / (2 * dz)
    
    # Forward difference at the first node
    dsdz[0] = (s[1] - s[0]) / dz
    
    # Backward difference at the last node
    dsdz[m-1] = (s[m-1] - s[m-2]) / dz
    
    return dsdz


#### Vapor pressure functions

def saturation_vapor_pressure(T):
    # Constants
    e0 = 611.2  # Reference vapor pressure in Pa at T0
    T0 = 273.15  # Reference temperature in K (0°C)
    Lv =  2.501e6  # Latent heat of vaporization for water in J/kg
    Rv = 461.5  # Specific gas constant for water vapor in J/(kg·K)
    
    # Saturation vapor pressure calculation
    es = e0 * np.exp((Lv / Rv) * (1/T0 - 1/T))
    return es



def vapor_concentration_RH( T, RH = 100):
    """
    Calculate the vapor concentration (in g/m³) for a given relative humidity
    (RH) and temperature (T in °C).
    
    Parameters:
    - RH (float): Relative humidity as a fraction (e.g., 0.5 for 50%).
    - T (float): Temperature in degrees Celsius.
    
    Returns:
    - float: Vapor concentration in g/m³.
    """
    # Constants
    M_w = 18.015  # molar mass of water vapor in g/mol
    R = 8.314  # universal gas constant in J/(mol*K)
    
    # Calculate the saturation vapor pressure (Pa)
    e_s = 0.6108 * math.exp((17.27 * T) / (T + 237.3))*1000
    
    # Actual vapor pressure (Pa)
    e = RH/100 * e_s
    
    # Vapor concentration (g/m³)
    C = (e * M_w) / (R * (T + 273.15))
    
    return C

def vapor_concentration(es, T):
    # Constants
    M = 18.015  # Molar mass of water in g/mol
    R = 8.314  # Universal gas constant in J/(mol·K)
    
    # Water vapor concentration calculation
    rho = (es * M) / (R * T)
    return rho


def e_sat(T):
    # Calculate the saturation vapor pressure (Pa)
    e_s = 0.6108 * math.exp((17.27 * T) / (T + 237.3))*1000    
    return e_s


def integrate_H2O_implicit(n, m, dx, dz, A, B, C, Qup, Qs, Qa):
    """
    This function computes the coefficients
      of the second-order system of variable s(t+dt)
      formed by the (known) profile values at s(t).
      
      The system is AA1 d2s/dz2 + AA2 ds/dz + AA3 s=AA4
      The coefficients in AA4 are used from s(t)
      The coefficients in AA1 are mainly the diffusion
      The coefficients in AA2 include advection
      The coefficients in AA3 include advection + s(t+dt)  
    """
    # Setup the tridiagonal solver for mean H2O concentration
    AA1 = - A * B
    AA2 = - C * B
    AA3 = 1 / dx
    AA4 = Qup / dx
    
    upd = (AA1 / (dz ** 2) + AA2 / (2 * dz))
    dia = (-2 * AA1 / (dz ** 2) + AA3)
    lod = (AA1 / (dz ** 2) - AA2 / (2 * dz))
    
    co = np.zeros(m)
    co[:] = AA4

    # Ensure the boundary conditions are state, not flux
    lod[0] = 0
    lod[m-1] = 0
    dia[0] = 1
    dia[m-1] = 1
    upd[0] = 0
    upd[m-1] = 0

    # Enforce surface and upper H2O concentration
    co[0] = Qs
    co[m-1] = Qa

    # Call the tridiagonal solver
    Q1 = thomas(lod, dia, upd, co)
    dQdz = our_central_difference(Q1, dz)

    # Fq = - A * dQdz
    Fq = np.zeros_like(Q1)    
    Fq[1:] = - 0.5 * (A[1:] + A[:-1]) * (0*dQdz[:-1] + dQdz[1:]) * 0.5

    Fq[1] = Fq[2]
    Fq[0] = Fq[1]    
    
    return Q1, Fq


def integrate_T_implicit(n, m, dx, dz, A, B, C, Tup, Ts, Ta):
    """
    This function computes the coefficients
      of the second-order system of variable s(t+dt)
      formed by the (known) profile values at s(t).
      
      The system is AA1 d2s/dz2 + AA2 ds/dz + AA3 s=AA4
      The coefficients in AA4 are used from s(t)
      The coefficients in AA1 are mainly the diffusion
      The coefficients in AA2 include advection
      The coefficients in AA3 include advection + s(t+dt)  
    """
    # Setup the tridiagonal solver for mean H2O concentration
    AA1 = - A * B
    AA2 = - C * B
    AA3 = 1 / dx
    AA4 = Tup / dx
    
    upd = (AA1 / (dz ** 2) + AA2 / (2 * dz))
    dia = (-2 * AA1 / (dz ** 2) + AA3)
    lod = (AA1 / (dz ** 2) - AA2 / (2 * dz))
    
    co = np.zeros(m)
    co[:] = AA4

    # Ensure the boundary conditions are state, not flux
    lod[0] = 0
    lod[m-1] = 0
    dia[0] = 1
    dia[m-1] = 1
    upd[0] = 0
    upd[m-1] = 0

    # Enforce surface and upper H2O concentration
    co[0] = Ts
    co[-1] = Ta

    # Call the tridiagonal solver
    T1 = thomas(lod, dia, upd, co)
    dTdz = our_central_difference(T1, dz)
    
    FT = np.zeros_like(T1)
    FT[1:] = - 0.5 * (A[1:] + A[:-1]) * (0*dTdz[:-1] + dTdz[1:]) * 0.5
    
    FT[1] = FT[2]
    FT[0] = FT[1]
    

    # FT = - A * dTdz
    return T1, FT


##### Stability functions

def stability(xi):
    """
    Stability functions from Brutsaert for stable and unstable conditions
    """
    xi = np.asarray(xi)  
    phi = np.zeros_like(xi)

    inds = np.where(xi <= 1/16)
    phi[inds]  =  (1 - 16 * xi[inds]) ** -0.5
    inds = np.where(xi > 1/16)
    phi[inds]  = 1 + 5.2 * xi[inds]

    return phi


def padit(x,  nz, upwind, array):
    """
    """
    xx = np.concatenate([np.arange(-40, 0, 2), x])

    up = upwind.reshape(1, nz)
    up = np.concatenate([up, up, up, up, up,up, up, up, up, up])

    array = np.concatenate([up, up, array])

    return xx, array


#### Data functions

from dataclasses import dataclass, asdict
import numpy as np

from dataclasses import dataclass, asdict
import numpy as np

@dataclass
class Params:
    """
    Container for simulation parameters with derived fields.

    Notes:
    - `fallow_size` and `field_size` are derived (in grid cells) from
      `fallow_length`, `dx`, and `fallow_fraction` 
    - `Q_c`, `Q_f`, and `Q_a` (water vapor concentrations) are added in `to_dict()`,
      using external helper functions:
          saturation_vapor_pressure(T_K)  # expects Kelvin
          vapor_concentration_RH(T_C, RH_percent)
    """

    # --- Physical constants / surface roughness ---
    k: float = 0.4                 # von Karman constant
    zom_f: float = 0.005           # Momentum roughness length (m), fallow
    zom_c: float = 0.05            # Momentum roughness length (m), cultivated

    # --- Canopy/obstacle height etc. ---
    h: float = 2.3                 # Canopy height or characteristic obstacle height (m)

    # --- Domain and grid ---
    Lx: float = 2000.0             # Domain length in x-direction (m)
    Hmax: float = 100.0            # Maximum height (m)
    dz: float = 0.5                # Vertical grid spacing (m)
    dx: float = 5.0                # Horizontal grid spacing (m)
    xmin: float = 0.0              # Domain start in x (m)

    # --- Patch geometry (inputs that define derived sizes) ---
    fallow_fraction: float = 0.5   # Fraction of fallow area relative to cultivated
    fallow_length: float = 1000.0  # Length of a single fallow patch in meters (physical)

    # --- Thermodynamic state (temperatures in °C, RH in %) ---
    T_sc: float = 30.0             # Cultivated surface temperature (°C)
    T_sf: float = 50.0             # Fallow surface temperature (°C)
    T_a: float  = 18.0             # Atmosphere (upwind) air temperature (°C)
    RH_c: float = 43.0             # Cultivated relative humidity (%)
    RH_f: float = 20.0             # Fallow relative humidity (%)
    RH_a: float = 18.0             # Atmospheric/background relative humidity (%)

    # --- Radiation / emissivity / albedo ---
    SW_in: float = 400.0           # Incoming shortwave radiation (W m^-2) [if used elsewhere]
    e_f: float = 0.95              # Surface emissivity (fallow)
    e_c: float = 0.95              # Surface emissivity (cultivated)
    e_a: float = 0.80              # Atmospheric emissivity
    alpha_c: float = 0.22          # Surface albedo (cultivated)
    alpha_f: float = 0.22          # Surface albedo (fallow)

    # --- Optional model switch(es) ---
    lm_option: str | None = 'kz'   # Mixing-length / closure scheme flag (e.g., 'kz')

    # -----------------------------
    # Derived helpers (not stored)
    # -----------------------------
    def _grid(self):
        """
        Build grid arrays and counts from domain extents and spacings.
        """
        zmin = self.zom_f               # lowest reference height ~ roughness length
        zmax = self.Hmax
        xmax = self.Lx

        # Standard grid definitions (match prior notebooks):
        x = np.arange(self.xmin, xmax + self.dx, self.dx)
        z = np.arange(zmin + self.dz, zmax + self.dz, self.dz)

        return {
            "x": x,
            "z": z,
            "nx": len(x),
            "nz": len(z),
            "xmax": xmax,
            "zmin": zmin,
            "zmax": zmax,
        }

    def _ustars(self):
        """
        Compute friction velocities (m s^-1) from a reference mean wind at 4 m
        using the log-law, consistent with older code:
            u*_f = k * Ubar_4 / ln(z_ref / zom_f)
            u*_c = k * Ubar_4 / ln(z_ref / zom_c)
        """
        Ubar_4 = 3.76                 # Reference mean wind speed at 4 m (m s^-1)
        z_ref  = 4.0                  # Reference height (m)
        ustar_f = self.k * Ubar_4 / np.log(z_ref / self.zom_f)   # Friction velocity (m/s), fallow
        ustar_c = self.k * Ubar_4 / np.log(z_ref / self.zom_c)   # Friction velocity (m/s), cultivated
        return {"ustar_f": ustar_f, "ustar_c": ustar_c}

    def _block_sizes(self):
        """
        Compute canonical block sizes (grid cells), preserving old logic:
            fallow_size = int(fallow_length / dx)
            field_size  = int(fallow_size * (1 - fallow_fraction) / fallow_fraction) 
                            if fallow_fraction > 0 else 0
        """
        fallow_size = int(self.fallow_length / self.dx)
        if self.fallow_fraction > 0:
            field_size = int(fallow_size * (1 - self.fallow_fraction) / self.fallow_fraction)
        else:
            field_size = 0
        
        return {"fallow_size": fallow_size, "field_size": field_size}

    def to_dict(self):
        """
        Export a flat dict like your old `get_params`, including derived and Q-fields.

        Adds:
          - 'ustar_f', 'ustar_c' : Friction velocities (m/s)
          - 'x','z','nx','nz','xmax','zmin','zmax' : Grid and extents
          - 'fallow_size','field_size' : Block sizes (grid cells; canonical naming)
          - 'Q_c','Q_f','Q_a' : Water vapor concentrations
                Q_c: Surface water vapor concentration for cultivated areas
                Q_f: Surface water vapor concentration for fallowed areas
                Q_a: Upwind background atmospheric water vapor concentration

        Requires the helper functions to be in scope:
            saturation_vapor_pressure(T_K) and vapor_concentration_RH(T_C, RH_%)
        """
        d = asdict(self)
        d.update(self._grid())
        d.update(self._ustars())
        d.update(self._block_sizes())

        # Moisture-related diagnostics (same as your old code)
        # Convert °C -> K for saturation_vapor_pressure
        Tsc_K = self.T_sc + 273.15
        Tsf_K = self.T_sf + 273.15

        # These functions are expected to exist in the runtime (as in your notebooks)
        es_c = saturation_vapor_pressure(Tsc_K)         # (Pa) Saturation vapor pressure at cultivated surface
        es_f = saturation_vapor_pressure(Tsf_K)         # (Pa) Saturation vapor pressure at fallow surface
        Q_c  = vapor_concentration_RH(self.T_sc, self.RH_c)  # (g/m^3) cultivated
        Q_f  = vapor_concentration_RH(self.T_sf, self.RH_f)  # (g/m^3) fallow
        Q_a  = vapor_concentration_RH(self.T_a,  self.RH_a)  # (g/m^3) atmospheric/upwind

        # Store Q-fields with the same interpretations as before:
        d["Q_c"] = Q_c   # Surface water vapor concentration for cultivated areas
        d["Q_f"] = Q_f   # Surface water vapor concentration for fallowed areas
        d["Q_a"] = Q_a   # Upwind background atmospheric water vapor concentration

        # (Optionally keep these around if you use them downstream)
        d["es_c"] = es_c # Saturation vapor pressure at cultivated surface (Pa)
        d["es_f"] = es_f # Saturation vapor pressure at fallow surface (Pa)

        return d
