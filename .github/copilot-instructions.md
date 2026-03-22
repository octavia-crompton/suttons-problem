# Sutton's Problem — AI Agent Instructions

> **Last updated:** 2026-03-01
> Read this file at the start of every new chat session.
> Update this file when new conventions, preferences, or project structure changes are established in chat.

---

## 1. Project Overview

This project numerically solves Sutton's advection problem — a 2-D boundary-layer transport model describing how heat and water vapour are redistributed downwind of an abrupt surface discontinuity (e.g. irrigated field next to fallow land). The Python codebase re-implements and extends the original FORTRAN/MATLAB solutions (see `CLOSET.FOR`, `matlab_versions/`), with the goal of reproducing published results (Rao 1974, Baldocchi 1995) and exploring advective enhancement of evapotranspiration over fragmented agricultural landscapes.

### Workspace layout

```
suttons_problem/
├── .github/copilot-instructions.md   ← THIS FILE
├── src/
│   └── sutton/                        ← core Python package
│       ├── __init__.py                ← public API exports
│       ├── config.py                  ← Params dataclass
│       ├── integrators.py             ← integrate_T_implicit, integrate_H2O_implicit
│       ├── numerics.py                ← thomas, thomas_alt, our_central_difference
│       ├── physics.py                 ← saturation_vapor_pressure, e_sat, vapor_concentration*
│       ├── constants.py               ← SIGMA_SB, CP_AIR, RHO_AIR, LV, LV_G, RV
│       ├── stability.py               ← stability()
│       └── utils.py                   ← padit()
│   └── sutton_functions.py            ← DEPRECATED shim; do not use in new code
├── notebooks/                         ← Jupyter notebooks for analysis & figures
│   └── archive/                       ← old / dated notebooks
├── data/                              ← input data files
├── matlab_versions/                   ← reference MATLAB implementations
├── tests/
│   ├── test_advection_params.py
│   ├── test_numerics.py
│   └── test_params.py
├── CLOSET.FOR                         ← original FORTRAN reference
└── pyproject.toml
```

### Key symbols

| Python name | Meaning |
|---|---|
| `AdvectionParams` | Extended dataclass: per-patch canopy, energy partition, radiation closure, BCs |
| `Params` | Legacy dataclass (grid, surface, meteorology); prefer `AdvectionParams` for new code |
| `SIGMA_SB`, `CP_AIR`, `RHO_AIR`, `LV`, `LV_G`, `RV` | Physical constants (in `constants.py`) |
| `integrate_T_implicit` | Implicit time-stepping for temperature advection |
| `integrate_H2O_implicit` | Implicit time-stepping for water vapour advection |
| `thomas` / `thomas_alt` | Thomas algorithm (tridiagonal solver) |
| `our_central_difference` | Central-difference advection scheme |
| `stability` | Stability correction functions |
| `e_sat` / `saturation_vapor_pressure` | Saturation vapour pressure functions |
| `padit` | Utility to pad/extend arrays |

---

## 2. Package Conventions

### Imports in notebooks

Always import from the `sutton` package, not from `sutton_functions.py` (deprecated shim).

```python
import sys as _sys
_sys.path.insert(0, "/Users/octaviacrompton/Projects/sutton_advection/src")
from sutton import (
    AdvectionParams,
    Params,
    thomas,
    our_central_difference,
    no_central_difference,
    saturation_vapor_pressure,
    vapor_concentration_RH,
    SIGMA_SB, CP_AIR, RHO_AIR, LV, LV_G, RV,
)
```

Or, since the package is installed in editable mode (`pip install -e .`), you can just:

```python
from sutton import AdvectionParams, integrate_T_implicit, integrate_H2O_implicit
```

### Python environment

- Python: ≥ 3.10 (uses `str | None` union syntax)
- Key packages: `numpy`, `matplotlib`, `scipy`, `xarray`, `pandas`, `jupyter`
- Package is installable via `pip install -e .` from the repo root

### `Params` dataclass

`Params` centralises every model parameter. Important fields:

| Field | Default | Meaning |
|---|---|---|
| `Lx` | 2000.0 m | Domain length in x |
| `Hmax` | 100.0 m | Domain height |
| `dz` | 0.5 m | Vertical grid spacing |
| `dx` | 5.0 m | Horizontal grid spacing |
| `fallow_fraction` | 0.5 | Fraction of domain that is fallow |
| `fallow_length` | 1000.0 m | Length of fallow patch |
| `T_sc`, `T_sf`, `T_a` | 30, 50, 18 °C | Crop/fallow/air temperatures |
| `RH_c`, `RH_f`, `RH_a` | 43, 20, 18 % | Crop/fallow/air relative humidities |
| `lm_option` | `"kz"` | Mixing-length option |

Use `Params()` to get all defaults; override specific fields as keyword arguments.

---

## 3. Cell Structure Preferences

- One logical task per cell.
- Markdown cells before major sections.
- Use `df.copy()` before any in-place mutations.
- Use `_private` leading-underscore names for cell-local temporaries.
- Keep parameter overrides close to the cell that uses them; do not scatter `Params` modifications across many cells.

---

## 4. Figures

There is no formal figure registry in this project yet. Current practice:

- Exploratory/scratch figures: save to `notebooks/` or a local path; no strict naming convention required.
- When publication figures are created, they should be named descriptively: `fig<N>_<description>.png`.
- Use `matplotlib` with explicit `fig, ax = plt.subplots(...)` rather than the implicit state machine.
- Label axes with physical units (e.g. `"x (m)"`, `"z (m)"`, `"T (°C)"`).
- For 2-D field plots, use `ax.contourf` or `ax.pcolormesh`; add a colorbar with a labelled unit.

---

## 5. Reference Implementations

- **FORTRAN:** `CLOSET.FOR` — original reference; parameters sourced from Rao et al.
- **MATLAB:** `matlab_versions/Matlab/` — working multi-variable version; `matlab_versions/sutton_only_wv/` — water-vapour-only version.
- When debugging numerical results, cross-check against MATLAB outputs or the Rao (1974) figures in `data/Rao_1974_Figure*.csv`.

---

## 6. Tests

Run with:

```bash
pytest
```

Tests live in `tests/`. Add tests for any new numerical scheme or physics function in `src/sutton/`.

---

## 7. Notebook File Management

- **Active notebooks** live in `notebooks/`.
- **Archived notebooks** go to `notebooks/archive/` with a date suffix if not already dated.
- Do not create new notebooks without being asked — modify existing ones.

---

## 8. Editing Preferences

- When editing `.ipynb` files, use `replace_string_in_file` with exact text matching (include 3–5 lines of context), or `edit_notebook_file` for cell-level edits.
- For complex multi-line notebook cell edits that fail with text matching, fall back to a Python script using `json.load` / `json.dump`.
- Always use `df.copy()` before any operation that modifies a subset in-place.
- Prefer `_private` names (leading underscore) for cell-local temporaries.
- Do not use `sutton_functions.py` in new code — import directly from `sutton`.

---

## 9. Updating This File

**This file should be updated whenever:**
- New modules or public functions are added to `src/sutton/`
- New `Params` fields are added
- A figure or data convention is established
- A new dataset is added to `data/`
- Test or CI conventions change
