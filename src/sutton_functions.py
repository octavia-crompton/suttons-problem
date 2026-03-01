"""
Deprecated shim. Use the package under src/sutton instead, for example:
    from sutton import Params, integrate_T_implicit

This module re-exports the public API to avoid breaking existing notebooks.
"""

from pathlib import Path
import sys
import warnings

# Ensure local 'src' is on sys.path for direct notebook use
_SRC = Path(__file__).parent / "src"
if _SRC.exists():
    sys.path.insert(0, str(_SRC))

warnings.warn(
    "sutton_functions.py is deprecated; import from the 'sutton' package (src/sutton) instead.",
    DeprecationWarning,
    stacklevel=2,
)

# Re-export from the package
from sutton import (
    Params,
    thomas,
    thomas_alt,
    our_central_difference,
    saturation_vapor_pressure,
    vapor_concentration_RH,
    vapor_concentration,
    e_sat,
    integrate_H2O_implicit,
    integrate_T_implicit,
    stability,
    padit,
)

__all__ = [
    "Params",
    "thomas",
    "thomas_alt",
    "our_central_difference",
    "saturation_vapor_pressure",
    "vapor_concentration_RH",
    "vapor_concentration",
    "e_sat",
    "integrate_H2O_implicit",
    "integrate_T_implicit",
    "stability",
    "padit",
]
