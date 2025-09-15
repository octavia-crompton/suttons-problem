# Sutton's Problem

Refactor and modeling workspace for advection/stability experiments.

## Package layout
- `src/sutton/`: modular code
	- `numerics.py`, `integrators.py`, `physics.py`, `stability.py`, `config.py`, `utils.py`
	- `__init__.py` exposes a simple API: `from sutton import Params, integrate_T_implicit, ...`
- `sutton_functions.py`: temporary shim re-exporting from `sutton` (deprecated)
- `tests/`: minimal tests
## Install (editable)
```bash
pip install -e .[dev]
```

## Use
```python
from sutton import Params, integrate_T_implicit, integrate_H2O_implicit, stability
```

## Notebook migration
- Replace old local function copies and `get_params` with `Params().to_dict()`
- Replace `import sutton_functions as sf` with `from sutton import ...`

## Notes
- Large data and media are ignored via `.gitignore`.
