"""Integration test for AdvectionParams (migrated to src/sutton/config.py)."""
import numpy as np
import pytest

from sutton import AdvectionParams, SIGMA_SB, CP_AIR, RHO_AIR, LV, LV_G, RV


def test_construction_defaults():
    p = AdvectionParams()
    assert p.T_sc == 28.5   # default crop surface T
    assert p.T_sf == 40.0   # default fallow surface T
    assert p.T_a == 30.0   # default air temperature
    assert p.z.shape[0] > 0
    assert p.x.shape[0] > 0


def test_construction_custom():
    p = AdvectionParams(ustar_f=0.8, ustar_c=0.8, le_factor=1.4,
                        avail_ratio=1.4, Hmax=2, dz=0.005, dx=0.2)
    assert p.ustar_f == 0.8
    assert p.le_factor == 1.4
    assert p.Hmax == 2


def test_grid_shapes():
    p = AdvectionParams(Lx=100, dz=0.5, dx=1.0, Hmax=10)
    # z grid starts at lm_zshift, not 0, so length may differ by 1
    assert len(p.z) >= int(10 / 0.5) - 1
    assert len(p.x) >= int(100 / 1.0)


def test_radiation_closure():
    p = AdvectionParams(ustar_f=0.8, ustar_c=0.8, le_factor=1.4,
                        avail_ratio=1.4, Hmax=2, dz=0.005, dx=0.2)
    p.solve_surface_radiation_inplace(fix="alpha_c")
    assert p.SW_in > 0
    assert 0 < p.alpha_f < 1
    assert 0 < p.alpha_c < 1
    assert p.Rn_c is not None
    assert p.Rn_f is not None
    assert p.LW_down is not None


def test_update_surface_BCs():
    p = AdvectionParams(ustar_f=0.8, ustar_c=0.8)
    p.solve_surface_radiation_inplace(fix="alpha_c")
    T_sc_before = p.T_sc
    p.update_surface_BCs_from_reference(update_T=True, update_Q=True)
    # Surface BCs should be updated (may or may not change)
    assert isinstance(p.T_sc, float)
    assert isinstance(p.T_sf, float)


def test_sync_moisture():
    p = AdvectionParams(ustar_f=0.8, ustar_c=0.8)
    p.solve_surface_radiation_inplace(fix="alpha_c")
    p.update_surface_BCs_from_reference(update_T=True, update_Q=True)
    p.sync_surface_moisture_inplace("RH_from_Q")
    assert p.RH_c > 0
    assert p.RH_f > 0


def test_reconcile_heat():
    p = AdvectionParams(ustar_f=0.8, ustar_c=0.8)
    p.solve_surface_radiation_inplace(fix="alpha_c")
    p.update_surface_BCs_from_reference(update_T=True, update_Q=True)
    p.reconcile_heat_with_prescribed_surface_T(mode="update_ustar", ref="z_ref")
    assert p.ustar_f > 0
    # ustar_c can be 0 when H_c ≈ 0 (all energy goes to LE)
    assert p.ustar_c >= 0


def test_constants():
    assert SIGMA_SB == 5.670374419e-8
    assert CP_AIR == 1005.0
    assert RHO_AIR == 1.20
    assert LV == 2.43e6
    assert LV_G == 2430.0
    assert RV == 461.5


def test_lm_cached_property():
    p = AdvectionParams(Hmax=10, dz=0.5)
    lm = p.lm
    assert lm.shape == p.z.shape
    assert np.all(lm > 0)


def test_energy_partition():
    """Check that H + LE = RNmG for both patches."""
    p = AdvectionParams(ustar_f=0.8, ustar_c=0.8, le_factor=1.4,
                        avail_ratio=1.4)
    # Fallow (upwind)
    np.testing.assert_allclose(p.H_f + p.LE_f, p.RNmG_up, rtol=1e-10)
    # Crop (downwind)
    np.testing.assert_allclose(p.H_c + p.LE_c, p.RNmG_down, rtol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
