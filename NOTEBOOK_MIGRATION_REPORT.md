# Notebook Migration Report

This report lists notebooks updated to import from the `sutton` package and flags local function definitions that overlap with the package functions.

## Updated notebooks (added import cell)
- notebooks/advection - grapex.ipynb
- notebooks/advection - repeat patches.ipynb

## Overlapping local functions detected

- notebooks/advection - grapex.ipynb
  - def integrate_H2O_implicit(n, m, dx, dz, A, B, C, Qup, Q_s, Q_a)
    - Diff vs src/sutton/integrators.py:
      - Parameter names differ (Q_s/Q_a vs Qs/Qa) but signature is compatible.
      - Body appears equivalent (Thomas solver + central difference + flux).

- notebooks/advection - repeat patches.ipynb
  - def integrate_H2O_implicit(n, m, dx, dz, A, B, C, Qup, Q_s, Q_a)
    - Same deviation as above.

- notebooks/advection - rao 1974.ipynb
  - def our_central_difference(s, dz)
    - Body matches src/sutton/numerics.py central difference.

- notebooks/advection - baldocchi.ipynb
  - commented versions of integrate_H2O_implicit/integrate_T_implicit present; recommend removing or replacing calls with imports.

- notebooks/advection problem - temperature.ipynb
  - def integrate_T_implicit(...)
  - def integrate_H2O_implicit(...)
    - Recommend replacing with imports from sutton.integrators and aligning argument names (Tup/Ts/Ta, Qup/Qs/Qa).

- notebooks/advection problem - stability.ipynb
  - def integrate_H2O_implicit(params, A, B, C, Q_up, Qs, Qa)
  - def integrate_T_implicit(params, A, B, C, T_up, Ts, Ta)
    - Signatures take a params object and compute grid sizes internally; our package functions accept (n, m, dx, dz, ...). If needed, we can provide thin wrappers that accept Params.

## Recommended actions
- Replace local function definitions with imported versions:
  - from sutton import integrate_H2O_implicit, integrate_T_implicit, our_central_difference
- If notebook uses Param-based signatures, add small adapter cells or we can add wrappers in sutton:
  - e.g., def integrate_T_with_params(p, A,B,C,Tup,Ts,Ta): return integrate_T_implicit(p['nx'], p['nz'], p['dx'], p['dz'], A,B,C,Tup,Ts,Ta)
- Remove commented duplicate code blocks after verifying results.

## Next steps
- I can continue batch-inserting import cells for remaining notebooks and optionally remove local duplicates, or add wrappers in `sutton` for Param-based calls.
