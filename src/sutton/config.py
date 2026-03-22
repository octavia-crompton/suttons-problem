"""
Parameter dataclasses for Sutton's advection problem.

``AdvectionParams``
    Extended parameter class with per-patch canopy geometry, energy
    partition, radiation closure, and boundary-condition methods for
    the implicit advection solver.

``Params``
    Legacy parameter class (kept for backward compatibility with older
    notebooks).  New code should use ``AdvectionParams``.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, fields, InitVar
from functools import cached_property
from typing import Optional

import numpy as np

from .constants import SIGMA_SB, CP_AIR, RHO_AIR, LV_G, RV
from .physics import saturation_vapor_pressure, vapor_concentration_RH


# ---------------------------------------------------------------------------
# AdvectionParams  (migrated from *Advection - general.ipynb*, March 2026)
# ---------------------------------------------------------------------------

@dataclass
class AdvectionParams:
    """Extended parameter dataclass for the two-patch advection model.

    Stores all grid, surface, meteorology, and energy-partition parameters.
    Surface moisture is stored canonically as ``Q_c`` / ``Q_f`` (g m^-3);
    ``RH_c`` / ``RH_f`` are read-only ``@property`` accessors derived from
    (T, Q) and are never stored.

    Typical usage
    -------------
    >>> p = AdvectionParams(ustar_f=0.8, ustar_c=0.8, le_factor=1.4,
    ...                     avail_ratio=1.4, Hmax=2, dz=0.005, dx=0.2)
    >>> p.solve_surface_radiation_inplace(fix="alpha_c")
    >>> p.update_surface_BCs_from_reference(update_T=True, update_Q=True)
    >>> p.sync_surface_moisture_inplace("RH_from_Q")
    >>> p.reconcile_heat_with_prescribed_surface_T(mode="update_ustar",
    ...                                            ref="z_ref")
    >>> p.update_surface_BCs_from_reference(update_T=True, update_Q=True)
    """

    # --- core physics / numerics ---
    k: float = 0.4
    zom_f: float = 0.0005            # m  (fallow / upwind momentum roughness)
    zom_c: float = 0.0005            # m  (crop / downwind momentum roughness)

    # canopy & displacement inputs (per patch)
    # (If None -> derive from z0m via alpha_m = 0.1, disp_frac)
    h_f_opt: Optional[float] = None
    h_c_opt: Optional[float] = None
    d_f_opt: Optional[float] = None
    d_c_opt: Optional[float] = None

    alpha_m: float = 0.10            # z0m ~ alpha_m * h
    disp_frac: float = 0.67          # d   ~ disp_frac * h

    G: float = 0.0                   # ground heat flux  (W m^-2)

    # domain / grid
    Lx: float = 100.0                # m
    Hmax: float = 10.0               # m
    dz: float = 0.1                  # m
    dx: float = 1.0                  # m
    xmin: float = 0.0                # m

    # friction velocities (optional; None => compute via resolve_ustars)
    ustar_f: Optional[float] = 0.15
    ustar_c: Optional[float] = 0.15

    # block pattern
    fallow_fraction: float = 0.5
    fallow_length: float = 1000.0    # m

    # temperatures (deg C)
    T_sc: float = 28.5
    T_sf: float = 40.0
    T_a:  float = 30.0

    # ---- surface moisture (choose ONE per surface at init) ----
    # Provide either RH_* (%) or Q_* (g m^-3).  If neither is given,
    # defaults (RH_c=60 %, RH_f=10 %) are used.
    # RH_* are InitVar: accepted at __init__ but NOT stored as fields.
    # After construction, p.RH_c / p.RH_f are computed @property.
    RH_c: InitVar[Optional[float]] = None
    RH_f: InitVar[Optional[float]] = None
    Q_c:  Optional[float] = None     # g m^-3  (canonical stored field)
    Q_f:  Optional[float] = None     # g m^-3  (canonical stored field)

    Q_a: float = 6.0
    tie_Qa_to_Qf: bool = False

    # radiation inputs
    SW_in: float = 400.0
    epsilon_f: float = 0.90
    epsilon_c: float = 0.95
    epsilon_a: float = 0.80
    alpha_c: float = 0.22
    alpha_f: float = 0.22

    lm_option: str | None = "kz"
    lm_zshift: float = 0.0005        # m, additive shift; tuned to match Rao (1974)

    # reference wind for u* inference (if computed from log law)
    U_ref: float = 3.76              # m s^-1 at z_ref_wind
    z_ref_wind: float = 4.0          # m

    # upstream / base flux: Rao base case uses (Rn-G)_up = H_f
    H_f: float = 355.64              # ~ 356 W m^-2
    LE_f: float = 0.0                # upstream LE (W m^-2)

    # downstream available-energy ratio & partition
    avail_ratio: float = 1.0         # (Rn-G)_down / (Rn-G)_up
    le_factor: float = 1.0           # LE_c = le_factor * (Rn-G)_down

    # ------------------------------------------------------------------
    #  helpers & grids
    # ------------------------------------------------------------------

    @staticmethod
    def _K(Tc: float) -> float:
        """deg C -> K."""
        return Tc + 273.15

    # derived canopy height & displacement (per patch) ----------------

    @cached_property
    def h_f(self) -> float:
        if self.h_f_opt is not None:
            return float(self.h_f_opt)
        if self.d_f_opt is not None:
            return float(self.d_f_opt) / self.disp_frac
        return float(self.zom_f) / self.alpha_m

    @cached_property
    def h_c(self) -> float:
        if self.h_c_opt is not None:
            return float(self.h_c_opt)
        if self.d_c_opt is not None:
            return float(self.d_c_opt) / self.disp_frac
        return float(self.zom_c) / self.alpha_m

    @cached_property
    def d_f(self) -> float:
        if self.d_f_opt is not None:
            return float(self.d_f_opt)
        if self.h_f_opt is not None:
            return self.disp_frac * float(self.h_f_opt)
        return self.disp_frac * (float(self.zom_f) / self.alpha_m)

    @cached_property
    def d_c(self) -> float:
        if self.d_c_opt is not None:
            return float(self.d_c_opt)
        if self.h_c_opt is not None:
            return self.disp_frac * float(self.h_c_opt)
        return self.disp_frac * (float(self.zom_c) / self.alpha_m)

    # grid -------------------------------------------------------------

    @cached_property
    def zmin(self) -> float:
        """Conservative global lower bound: z > max(z0m_f, z0m_c)."""
        return max(self.zom_f, self.zom_c)

    @cached_property
    def z(self) -> np.ndarray:
        return np.arange(self.zmin, self.Hmax, self.dz)

    @cached_property
    def x(self) -> np.ndarray:
        return np.arange(self.xmin, self.Lx + self.dx, self.dx)

    @property
    def nz(self) -> int:
        return len(self.z)

    @property
    def nx(self) -> int:
        return len(self.x)

    # friction velocity ------------------------------------------------

    def resolve_ustars(self) -> "AdvectionParams":
        """Compute u* from log law if not already set."""
        if self.ustar_f is None:
            assert self.z_ref_wind > self.zom_f > 0
            self.ustar_f = self.k * self.U_ref / np.log(self.z_ref_wind / self.zom_f)
        if self.ustar_c is None:
            assert self.z_ref_wind > self.zom_c > 0
            self.ustar_c = self.k * self.U_ref / np.log(self.z_ref_wind / self.zom_c)
        return self

    # block sizes ------------------------------------------------------

    @cached_property
    def fallow_size(self) -> int:
        return int(self.fallow_length / self.dx)

    @cached_property
    def field_size(self) -> int:
        return (
            int(self.fallow_size
                * (1 - self.fallow_fraction) / self.fallow_fraction)
            if self.fallow_fraction > 0 else 0
        )

    # mixing length ----------------------------------------------------

    @cached_property
    def lm(self) -> np.ndarray:
        d = 2 / 3 * self.h_c
        z = self.z
        k = self.k
        h = self.h_c

        if self.lm_option == "kz":
            lm = k * (z + self.lm_zshift)
        else:
            a = 1
            lm = k * z.copy()
            lm[z < a * (h - d)] = k * a * h / 3
        return lm

    # saturation vapour pressure (cached) ------------------------------

    @cached_property
    def es_c(self) -> float:
        return saturation_vapor_pressure(self._K(self.T_sc))

    @cached_property
    def es_f(self) -> float:
        return saturation_vapor_pressure(self._K(self.T_sf))

    # RH: always derived from canonical (T, Q) state -------------------

    @property
    def RH_c(self) -> float:
        """Relative humidity at cultivated surface (%), from (T_sc, Q_c)."""
        return self.RH_from_Q(self.T_sc, self.Q_c)

    @property
    def RH_f(self) -> float:
        """Relative humidity at fallow surface (%), from (T_sf, Q_f)."""
        return self.RH_from_Q(self.T_sf, self.Q_f)

    # energy partition -------------------------------------------------

    @property
    def RNmG_up(self) -> float:
        """Upstream available energy (W m^-2)."""
        return self.H_f + self.LE_f

    @property
    def RNmG_down(self) -> float:
        """Downstream available energy = avail_ratio * RNmG_up (W m^-2)."""
        return self.RNmG_up * self.avail_ratio

    @property
    def LE_c(self) -> float:
        """Downstream latent heat flux (W m^-2).

        ``le_factor > 1`` is valid and gives negative ``H_c``
        (evaporative cooling, T_sc < T_a).
        """
        return self.le_factor * self.RNmG_up * self.avail_ratio

    @property
    def H_c(self) -> float:
        """Downstream sensible heat flux (W m^-2)."""
        return (1.0 - self.le_factor) * self.RNmG_up * self.avail_ratio

    # ------------------------------------------------------------------
    #  to_dict (legacy call-sites)
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        d = asdict(self)
        d.update({
            "x": self.x,
            "z": self.z,
            "nx": self.nx,
            "nz": self.nz,
            "zmin": self.zmin,
            "zmax": self.Hmax,
            "xmax": self.Lx,
            "ustar_f": self.ustar_f,
            "ustar_c": self.ustar_c,
            "fallow_size": self.fallow_size,
            "field_size": self.field_size,
            "Q_c": self.Q_c,
            "Q_f": self.Q_f,
            "Q_a": self.Q_a,
            "RH_c": self.RH_c,
            "RH_f": self.RH_f,
            "es_c": self.es_c,
            "es_f": self.es_f,
            "LE_c": self.LE_c,
            "H_c": self.H_c,
            "RNmG_up": self.RNmG_up,
            "RNmG_down": self.RNmG_down,
            "h_f": self.h_f,
            "h_c": self.h_c,
            "d_f": self.d_f,
            "d_c": self.d_c,
        })
        return d

    # ------------------------------------------------------------------
    #  Radiation closure
    # ------------------------------------------------------------------

    def solve_surface_radiation_inplace(self, *, fix: str = "alpha_f") -> dict:
        """Solve the two-surface net-radiation balances, updating in place.

        The system is::

            Rn_f = SW*(1 - alpha_f) + eps_a * sigma * T_a^4 - eps_f * sigma * T_sf^4
            Rn_c = SW*(1 - alpha_c) + eps_a * sigma * T_a^4 - eps_c * sigma * T_sc^4

        with  Rn_f = RNmG_up + G  and  Rn_c = RNmG_down + G.

        One albedo is kept fixed (``fix``); the other and SW are solved.

        Parameters
        ----------
        fix : {"alpha_f", "alpha_c"}
            Which albedo to treat as prescribed.

        Returns
        -------
        dict with SW, alpha_f, alpha_c, Rn_f, Rn_c.
        """
        if fix not in ("alpha_f", "alpha_c"):
            raise ValueError("fix must be 'alpha_f' or 'alpha_c'.")

        Rn_f = float(self.RNmG_up + self.G)
        Rn_c = float(self.RNmG_down + self.G)

        TaK  = self._K(self.T_a)
        TsfK = self._K(self.T_sf)
        TscK = self._K(self.T_sc)

        sigma = SIGMA_SB
        ea, ef, ec = self.epsilon_a, self.epsilon_f, self.epsilon_c

        LW_down = ea * sigma * TaK ** 4
        LWup_f  = ef * sigma * TsfK ** 4
        LWup_c  = ec * sigma * TscK ** 4

        # Direct algebraic solve (2x2 linear system, no sympy needed):
        #   SW * (1 - alpha_fixed) = Rn_fixed - LW_down + LWup_fixed
        #   alpha_free  = 1 - (Rn_free - LW_down + LWup_free) / SW

        if fix == "alpha_f":
            denom = 1.0 - self.alpha_f
            if abs(denom) < 1e-15:
                raise RuntimeError("alpha_f ~ 1; cannot solve for SW.")
            SW_val = (Rn_f - LW_down + LWup_f) / denom
            if abs(SW_val) < 1e-15:
                raise RuntimeError("Implied SW ~ 0; cannot solve for alpha_c.")
            alpha_c_val = 1.0 - (Rn_c - LW_down + LWup_c) / SW_val
            alpha_c_val = min(max(alpha_c_val, 0.0), 1.0)

            self.SW_in   = SW_val
            self.alpha_c = alpha_c_val
        else:
            denom = 1.0 - self.alpha_c
            if abs(denom) < 1e-15:
                raise RuntimeError("alpha_c ~ 1; cannot solve for SW.")
            SW_val = (Rn_c - LW_down + LWup_c) / denom
            if abs(SW_val) < 1e-15:
                raise RuntimeError("Implied SW ~ 0; cannot solve for alpha_f.")
            alpha_f_val = 1.0 - (Rn_f - LW_down + LWup_f) / SW_val
            alpha_f_val = min(max(alpha_f_val, 0.0), 1.0)

            self.SW_in   = SW_val
            self.alpha_f = alpha_f_val

        # Dynamic attributes (not dataclass fields -- for diagnostic access)
        self.Rn_c    = Rn_c
        self.Rn_f    = Rn_f
        self.LW_down = LW_down

        return {
            "SW": SW_val,
            "alpha_f": float(self.alpha_f),
            "alpha_c": float(self.alpha_c),
            "Rn_f": Rn_f,
            "Rn_c": Rn_c,
        }

    # ------------------------------------------------------------------
    #  Boundary-condition updates
    # ------------------------------------------------------------------

    def update_surface_BCs_from_reference(
        self,
        *,
        update_T: bool = True,
        update_Q: bool = True,
    ) -> dict:
        """Update wall BCs from fluxes using a two-point neutral log profile.

        Reference air at ``z_ref = z[0]``, top at ``z_h = Hmax``::

            T_s = T_a + H / (rho cp u* k) * ln(z_h / z_ref)
            Q_s = Q_a + LE / (Lv u* k) * ln(z_h / z_ref)
        """
        z_h      = float(self.Hmax)
        z_ref    = float(self.z[0])
        kappa    = float(self.k)
        log_frac = np.log(z_h / z_ref)

        # cultivated (c)
        T_sc_new = self.T_a + (self.H_c  / (RHO_AIR * CP_AIR * self.ustar_c * kappa)) * log_frac
        Q_c_new  = self.Q_a + (self.LE_c / (LV_G * self.ustar_c * kappa)) * log_frac

        if update_T:
            self.T_sc = float(T_sc_new)
            self.__dict__.pop("es_c", None)
        if update_Q:
            self.Q_c = float(Q_c_new)

        # fallow (f)
        T_sf_new = self.T_a + (self.H_f  / (RHO_AIR * CP_AIR * self.ustar_f * kappa)) * log_frac
        Q_f_new  = self.Q_a + (self.LE_f / (LV_G * self.ustar_f * kappa)) * log_frac

        if update_T:
            self.T_sf = float(T_sf_new)
            self.__dict__.pop("es_f", None)
        if update_Q:
            self.Q_f = float(Q_f_new)
            self._apply_Qa_tie()

        return {
            "z_ref": z_ref,
            "T_sc": float(self.T_sc),
            "T_sf": float(self.T_sf),
            "Q_c":  float(self.Q_c),
            "Q_f":  float(self.Q_f),
        }

    def update_surface_BCs_from_zom(
        self,
        *,
        update_T: bool = True,
        update_Q: bool = True,
    ) -> dict:
        """Like ``update_surface_BCs_from_reference`` but uses z0h = z0m/5."""
        zmax  = float(self.Hmax)
        kappa = float(self.k)

        if self.ustar_c is None or self.ustar_f is None:
            raise ValueError("ustar_f/ustar_c is None; set them or call resolve_ustars().")

        z0h_c = self.zom_c / 5
        z0q_c = self.zom_c / 5
        z0h_f = self.zom_f / 5
        z0q_f = self.zom_f / 5

        # cultivated
        T_sc_new = self.T_a + (self.H_c  / (RHO_AIR * CP_AIR * self.ustar_c * kappa)) * np.log(zmax / z0h_c)
        Q_c_new  = self.Q_a + (self.LE_c / (LV_G * self.ustar_c * kappa)) * np.log(zmax / z0q_c)

        if update_T:
            self.T_sc = float(T_sc_new)
            self.__dict__.pop("es_c", None)
        if update_Q:
            self.Q_c = float(Q_c_new)

        # fallow
        T_sf_new = self.T_a + (self.H_f  / (RHO_AIR * CP_AIR * self.ustar_f * kappa)) * np.log(zmax / z0h_f)
        Q_f_new  = self.Q_a + (self.LE_f / (LV_G * self.ustar_f * kappa)) * np.log(zmax / z0q_f)

        if update_T:
            self.T_sf = float(T_sf_new)
            self.__dict__.pop("es_f", None)
        if update_Q:
            self.Q_f = float(Q_f_new)
            self._apply_Qa_tie()

        return {
            "z_ref": zmax,
            "T_sc": float(self.T_sc),
            "T_sf": float(self.T_sf),
            "Q_c":  float(self.Q_c),
            "Q_f":  float(self.Q_f),
        }

    # ------------------------------------------------------------------
    #  Heat reconciliation
    # ------------------------------------------------------------------

    def reconcile_heat_with_prescribed_surface_T(
        self,
        *,
        mode: str = "update_H",
        ref: str = "z_ref",
        eps_dT: float = 1e-6,
    ) -> dict:
        """Enforce consistency between (T_sf, T_sc), u*, and fluxes.

        Parameters
        ----------
        mode : {"update_H", "update_ustar"}
            ``"update_H"``    -- keep u* fixed, update H_f and le_factor.
            ``"update_ustar"`` -- keep H_f and le_factor fixed, update u*.
        ref : {"z_ref", "z0h"}
            ``"z_ref"`` uses ln(Hmax / z[0]) (matches uniform_T).
            ``"z0h"``   uses ln(Hmax / (z0m/5)) per surface.
        """
        z_h   = float(self.Hmax)
        kappa = float(self.k)

        def log_fac(surface: str) -> float:
            if ref == "z_ref":
                z_ref = float(self.z[0])
                if z_h <= z_ref:
                    raise ValueError(f"Hmax={z_h} must be > z_ref={z_ref}.")
                return float(np.log(z_h / z_ref))
            elif ref == "z0h":
                z0h = (self.zom_f / 5.0) if surface == "f" else (self.zom_c / 5.0)
                if z_h <= z0h:
                    raise ValueError(f"Hmax={z_h} must be > z0h={z0h} for surface={surface}.")
                return float(np.log(z_h / z0h))
            else:
                raise ValueError("ref must be 'z_ref' or 'z0h'.")

        dTf = float(self.T_sf - self.T_a)
        dTc = float(self.T_sc - self.T_a)
        Lf  = log_fac("f")
        Lc  = log_fac("c")

        if mode == "update_H":
            if self.ustar_f is None:
                raise ValueError("ustar_f is None; set it before mode='update_H'.")
            Hf_new = RHO_AIR * CP_AIR * float(self.ustar_f) * kappa * dTf / Lf
            self.H_f = float(Hf_new)

            if self.ustar_c is None:
                raise ValueError("ustar_c is None; set it before mode='update_H'.")
            Hc_target = RHO_AIR * CP_AIR * float(self.ustar_c) * kappa * dTc / Lc
            RNmG_down = float(self.RNmG_up * self.avail_ratio)
            if abs(RNmG_down) < 1e-12:
                raise ValueError("RNmG_down is ~0; cannot solve for le_factor.")
            self.le_factor = float(1.0 - (Hc_target / RNmG_down))

            return {
                "mode": mode, "ref": ref,
                "H_f_new": float(self.H_f),
                "H_c_target": float(Hc_target),
                "le_factor_new": float(self.le_factor),
                "RNmG_up": float(self.RNmG_up),
                "RNmG_down": float(RNmG_down),
            }

        elif mode == "update_ustar":
            # ustar_f
            if abs(dTf) < eps_dT:
                if abs(self.H_f) > 1e-6:
                    raise ValueError("T_sf ~ T_a but H_f is nonzero.")
                ustar_f_new = 0.0
            else:
                ustar_f_new = float(self.H_f) * Lf / (RHO_AIR * CP_AIR * kappa * dTf)
                if ustar_f_new < 0:
                    raise ValueError("Implied ustar_f < 0 (sign mismatch).")
            self.ustar_f = float(ustar_f_new)

            # ustar_c
            Hc_now = float(self.H_c)
            if abs(dTc) < eps_dT:
                if abs(Hc_now) > 1e-6:
                    raise ValueError("T_sc ~ T_a but H_c is nonzero.")
                ustar_c_new = 0.0
            else:
                ustar_c_new = Hc_now * Lc / (RHO_AIR * CP_AIR * kappa * dTc)
                if ustar_c_new < 0:
                    raise ValueError("Implied ustar_c < 0 (sign mismatch).")
            self.ustar_c = float(ustar_c_new)

            return {
                "mode": mode, "ref": ref,
                "ustar_f_new": float(self.ustar_f),
                "ustar_c_new": float(self.ustar_c),
                "H_f": float(self.H_f),
                "H_c": float(Hc_now),
            }
        else:
            raise ValueError("mode must be 'update_H' or 'update_ustar'.")

    # ------------------------------------------------------------------
    #  Moisture initialisation & helpers
    # ------------------------------------------------------------------

    def __post_init__(self, RH_c: Optional[float], RH_f: Optional[float]):
        # The @property descriptors for RH_c / RH_f shadow the InitVar
        # defaults.  When the caller omits RH_c / RH_f, the auto-generated
        # __init__ passes the property object instead of None.
        if isinstance(RH_c, property):
            RH_c = None
        if isinstance(RH_f, property):
            RH_f = None
        self._init_surface_moisture(RH_c=RH_c, RH_f=RH_f)
        self._apply_Qa_tie()

    def _apply_Qa_tie(self):
        """If requested, enforce Q_a = Q_f."""
        if self.tie_Qa_to_Qf:
            self.Q_a = float(self.Q_f)

    def _init_surface_moisture(
        self,
        *,
        RH_c: Optional[float],
        RH_f: Optional[float],
    ):
        self._coerce_surface_pair(label="c", Tc=self.T_sc, RH=RH_c,
                                  Q_name="Q_c", RH_default=60.0)
        self._coerce_surface_pair(label="f", Tc=self.T_sf, RH=RH_f,
                                  Q_name="Q_f", RH_default=10.0)

    def _coerce_surface_pair(
        self,
        *,
        label: str,
        Tc: float,
        RH: Optional[float],
        Q_name: str,
        RH_default: float,
        enforce_RH_bounds: bool = True,
    ):
        Q = getattr(self, Q_name)
        if RH is not None and Q is not None:
            raise ValueError(
                f"Specify only ONE of RH_{label} or {Q_name} (not both)."
            )
        if RH is None and Q is None:
            RH = RH_default

        if RH is not None:
            RH = float(RH)
            if enforce_RH_bounds and not (0.0 <= RH <= 100.0):
                raise ValueError(f"RH_{label} must be in [0,100], got {RH}.")
            Q = float(self.Q_from_RH(Tc, RH))
        else:
            Q = float(Q)
            if Q < 0.0:
                raise ValueError(f"{Q_name} must be >= 0, got {Q}.")
            RH_implied = float(self.RH_from_Q(Tc, Q))
            if enforce_RH_bounds and not (0.0 <= RH_implied <= 100.0):
                raise ValueError(
                    f"{Q_name}={Q} at T={Tc} deg C implies "
                    f"RH_{label}={RH_implied:.1f}% (outside [0,100])."
                )
        setattr(self, Q_name, Q)

    @staticmethod
    def Q_from_RH(Tc: float, RH_pct: float) -> float:
        """Q (g m^-3) from T (deg C) and RH (%)."""
        return float(vapor_concentration_RH(Tc, RH_pct))

    @staticmethod
    def RH_from_Q(Tc: float, Q_gm3: float) -> float:
        """RH (%) from T (deg C) and Q (g m^-3)."""
        Tk    = Tc + 273.15
        rho_v = Q_gm3 / 1000.0          # kg m^-3
        e     = rho_v * RV * Tk          # Pa
        es    = saturation_vapor_pressure(Tk)
        return 100.0 * (e / es)

    def sync_surface_moisture_inplace(
        self,
        direction: str = "RH_from_Q",
        clip_RH: bool = True,
        clip_Q: bool = True,
    ) -> dict:
        """Return current (T, Q)-consistent RH values for inspection.

        ``direction="RH_from_Q"`` is the only supported mode (RH is a
        read-only property; it is always consistent with Q).
        """
        direction = direction.strip()
        if direction == "RH_from_Q":
            RHc = float(np.clip(self.RH_c, 0.0, 100.0)) if clip_RH else self.RH_c
            RHf = float(np.clip(self.RH_f, 0.0, 100.0)) if clip_RH else self.RH_f
            return {"RH_c": RHc, "RH_f": RHf}
        elif direction == "Q_from_RH":
            raise TypeError(
                "sync_surface_moisture_inplace('Q_from_RH') is no longer supported.\n"
                "Use:  p.Q_c = p.Q_from_RH(p.T_sc, new_RH_c)"
            )
        else:
            raise ValueError("direction must be 'RH_from_Q' or 'Q_from_RH'.")

    # ------------------------------------------------------------------
    #  repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        keys = sorted(f.name for f in fields(self))
        body = ",\n  ".join(f"{k}={getattr(self, k)!r}" for k in keys)
        return f"{self.__class__.__name__}(\n  {body}\n)"

    def __str__(self) -> str:
        return self.__repr__()


# ---------------------------------------------------------------------------
# Legacy Params  (backward-compatible with older notebooks)
# ---------------------------------------------------------------------------

@dataclass
class Params:
    """Original (pre-AdvectionParams) parameter class.

    Kept for backward compatibility with ``Baldocchi - 1995``,
    ``Advection - grapex``, and other older notebooks.
    New code should use :class:`AdvectionParams`.
    """

    k: float = 0.4
    zom_f: float = 0.005
    zom_c: float = 0.05
    h: float = 2.3
    Lx: float = 2000.0
    Hmax: float = 100.0
    dz: float = 0.5
    dx: float = 5.0
    xmin: float = 0.0

    fallow_fraction: float = 0.5
    fallow_length: float = 1000.0

    T_sc: float = 30.0
    T_sf: float = 50.0
    T_a: float = 18.0
    RH_c: float = 43.0
    RH_f: float = 20.0
    RH_a: float = 18.0

    SW_in: float = 400.0
    e_f: float = 0.95
    e_c: float = 0.95
    e_a: float = 0.80
    alpha_c: float = 0.22
    alpha_f: float = 0.22

    lm_option: str | None = "kz"

    def _grid(self):
        zmin = self.zom_f
        zmax = self.Hmax
        xmax = self.Lx
        x = np.arange(self.xmin, xmax + self.dx, self.dx)
        z = np.arange(zmin + self.dz, zmax + self.dz, self.dz)
        return {
            "x": x, "z": z, "nx": len(x), "nz": len(z),
            "xmax": xmax, "zmin": zmin, "zmax": zmax,
        }

    def _ustars_from_4m(self):
        Ubar_4 = 3.76
        z_ref = 4.0
        ustar_f = self.k * Ubar_4 / np.log(z_ref / self.zom_f)
        ustar_c = self.k * Ubar_4 / np.log(z_ref / self.zom_c)
        return {"ustar_f": ustar_f, "ustar_c": ustar_c}

    def _block_sizes(self):
        fallow_size = int(self.fallow_length / self.dx)
        field_size = (
            int(fallow_size * (1 - self.fallow_fraction) / self.fallow_fraction)
            if self.fallow_fraction > 0 else 0
        )
        return {"fallow_size": fallow_size, "field_size": field_size}

    def to_dict(self):
        d = asdict(self)
        d.update(self._grid())
        d.update(self._ustars_from_4m())
        d.update(self._block_sizes())
        Tsc_K = self.T_sc + 273.15
        Tsf_K = self.T_sf + 273.15
        es_c = saturation_vapor_pressure(Tsc_K)
        es_f = saturation_vapor_pressure(Tsf_K)
        Q_c = vapor_concentration_RH(self.T_sc, self.RH_c)
        Q_f = vapor_concentration_RH(self.T_sf, self.RH_f)
        Q_a = vapor_concentration_RH(self.T_a, self.RH_a)
        d["Q_c"] = Q_c
        d["Q_f"] = Q_f
        d["Q_a"] = Q_a
        d["es_c"] = es_c
        d["es_f"] = es_f
        return d
