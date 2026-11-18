# 2i-OTCA-C v1.0 -- Two-Impulse Orbital Transfer Cost Algorithm -- Circular Case

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Status: implemented, v1.0.** This folder holds the scientific/operational
algorithm itself -- the solver library, its automated tests, a usage
example, an independent-comparison validation suite, and a fast surrogate
module -- as a reusable component, separate from `2i-OTCA-C-v1.0-demo/`
(the browser demonstrator) at the repository root.

## Contents

- `src/otca_c/` -- the Python solver library: `core.py` (cost function and
  the p* root-finder), `hohmann.py` (exact closed form for the Hohmann
  case), `geometry.py` (independent geometric re-derivation, used as one of
  the validation layers), `lambert.py` (an independent universal-variable
  Lambert solver, the second validation layer), `rk4.py` (direct numerical
  propagation, the third validation layer), `derived.py` (periods,
  altitudes, time of flight, propellant mass fraction), `units.py`,
  `errors.py` (`ConvergenceError`, raised rather than ever returning a
  silently wrong or clamped value), `surrogate.py` (a closed-form estimate
  of the transfer cost, valid within 30 degrees of the Hohmann angle, with
  an explicit out-of-domain flag), `budget_inversion.py` (the largest
  angular deviation from Hohmann a given delta-v budget affords, via a
  safeguarded Newton solver in `z := cot^2(phi/2)` with monotone
  bracketing in `phi`), `curvature_regime.py` (classifies a R=2 transfer
  relative to Theorem 4's certified curvature sign-change point; R=2
  only, by construction), and `api.py` (the public entry point tying all
  of the above together).
- `tests/` -- the automated test suite (`pytest`), including
  `test_js_python_parity.py`, which checks that the JavaScript port used by
  `2i-OTCA-C-v1.0-demo/otca.js` produces identical results to this Python
  library across a battery of cases -- the two are kept in sync by this
  test, not by hand.
- `Independent_comparison/` -- a second, independently written validation
  suite (grid comparison, k-convergence, an independent multiprecision
  Lambert solver, a global-minimum check, and the associated plots/results)
  cross-checking the main library's output.
- `visualize_transfer.py` -- a usage example producing
  `transfer_visualization.png`.
- `docs_atbd/` -- the Algorithm Theoretical Basis Document (ATBD) for the
  surrogate module (`ATBD_2i-OTCA-C_surrogate_v1.0.docx`), companion ATBDs
  for the `budget_inversion.py` and `curvature_regime.py` modules
  (`ATBD_2i-OTCA-C_budget_inversion_v1.0.docx`,
  `ATBD_2i-OTCA-C_curvature_regime_v1.0.docx`), the full-suite
  verification summary, and the module/theorem correspondence table.
- `pyproject.toml` -- packaging metadata for the `otca_c` library.

## Scope

### In scope

- Circular-orbit transfers within the mathematical model established in the
  accompanying proofs (two-impulse, coplanar, circular orbits).
- Evaluation of the cost function for any supported input parameters.
- Selection of the provably unique, provably global minimum-cost transfer
  for a chosen transfer angle (the existence/uniqueness result, Theorem 1,
  Theorem 2 -- see `circular/` at the repository root).
- Three independent validation layers on every solved transfer: a geometry
  check, a Lambert cross-check, and RK4 numerical propagation.
- A fast closed-form surrogate for transfers within 30 degrees of the
  Hohmann angle, with an explicit flag whenever a query falls outside that
  measured domain (`surrogate.py`).
- Inversion of a delta-v budget into the largest affordable angular
  deviation from Hohmann (`budget_inversion.py`), for every supported R.
- Classification of a R=2 transfer's position relative to Theorem 4's
  certified curvature sign-change point (`curvature_regime.py`), R=2 only.

### Explicitly out of scope

This component does **not** model or optimize:

- Elliptical-orbit transfers.
- Orbital-plane changes.
- Perturbations such as J2.
- Finite-duration or finite-thrust propulsion.
- Launch-window, phasing, rendezvous, or mission-scheduling constraints.
- Operational mission planning involving constraints not represented in the
  mathematical model.

2i-OTCA-C is a research software component, not a general-purpose mission
planner.

## Mathematical contract

Exact inputs, the angle convention, the exact expression evaluated, the
domain over which the existence/uniqueness/monotonicity theorems guarantee
the result, and the defined behaviour for out-of-domain or invalid input are
all specified in `docs/mathematical_conventions.md` at the repository root.

## Running the tests

```
pip install -e .
pytest tests/
```
