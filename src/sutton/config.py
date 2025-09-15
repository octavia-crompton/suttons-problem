from dataclasses import dataclass, asdict
import numpy as np
from .physics import saturation_vapor_pressure, vapor_concentration_RH


@dataclass
class Params:
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
        return {"x": x, "z": z, "nx": len(x), "nz": len(z), "xmax": xmax, "zmin": zmin, "zmax": zmax}

    def _ustars(self):
        Ubar_4 = 3.76
        z_ref = 4.0
        ustar_f = self.k * Ubar_4 / np.log(z_ref / self.zom_f)
        ustar_c = self.k * Ubar_4 / np.log(z_ref / self.zom_c)
        return {"ustar_f": ustar_f, "ustar_c": ustar_c}

    def _block_sizes(self):
        fallow_size = int(self.fallow_length / self.dx)
        field_size = int(fallow_size * (1 - self.fallow_fraction) / self.fallow_fraction) if self.fallow_fraction > 0 else 0
        return {"fallow_size": fallow_size, "field_size": field_size}

    def to_dict(self):
        d = asdict(self)
        d.update(self._grid())
        d.update(self._ustars())
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
