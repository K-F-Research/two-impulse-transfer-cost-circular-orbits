Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

# Verification Summary — otca_c (2i-OTCA-C v1.0, package version 0.1.0)

Generated from a direct execution of the automated test suite. Figures in
this document are taken verbatim from that execution; none are retyped by
hand. Re-generate with:

```
pip install -e .
pytest tests/ -v
```

run from the repository root, with the `2i-OTCA-C-v1.0-algorithm/` and
`2i-OTCA-C-v1.0-demo/` folders both present as siblings (required for the
JavaScript/Python parity tests, Section 2) and with `circular/` present as
a sibling of `2i-OTCA-C-v1.0-algorithm/` (required for the curvature
sign-change pin test, Section 2).

## 1. Scope

This document reports the result of the complete automated test suite for
the `otca_c` package, as of the test run identified in Section 5. It
should be read together with the Algorithm Theoretical Basis Document
(ATBD) for this component, which this summary is attached to as an annex.

Three categories of automated check are distinguished, consistent with
the distinction already established in ATBD Section 9:

- **Unit and reference-value tests**: a specific input/output pair is
  checked against a value stated in the manuscript or computed by an
  independent closed-form derivation.
- **Metamorphic / property tests**: an invariant that must hold across a
  family of inputs by construction of the underlying physics (scaling
  laws, normalization independence) is checked directly, not assumed.
- **Regression tests against an archived snapshot**: the current output
  of `otca_c.api.compute_transfer` is compared against a numerical
  snapshot recorded at a prior date (`Independent_comparison/
  grid_compare_full_results.json`), to detect unintended future changes
  to the implementation. This is distinct from, and does not repeat, the
  independent cross-check documented in
  `Independent_comparison/FINAL_REPORT.md`.

## 2. Result, by test file

| Test file | Test cases collected | Result | Category |
|---|---:|---|---|
| `test_api.py` | 6 | 6 passed | Unit (public API) |
| `test_budget_inversion.py` | 29 | 29 passed | Unit + metamorphic (budget-to-angle inversion; Theorem 1, Theorem 3) |
| `test_core.py` | 12 | 12 passed | Unit (cost function, p* solver) |
| `test_curvature_regime.py` | 6 | 6 passed | Unit (R=2 curvature-regime classifier; Theorem 4), including a pin of the module's certified constants against `circular/theorem_4/03_central_mesh/mesh_certificate_R2.py` |
| `test_derived.py` | 20 | 20 passed | Unit (periods, altitudes, propellant fraction) |
| `test_geometry.py` | 35 | 35 passed | Unit (independent geometry cross-check layer) |
| `test_hohmann_limit.py` | 9 | 9 passed | Unit (phi -> 0 limit) |
| `test_js_python_parity.py` | 33 | 33 passed | Unit (Python/JavaScript parity) |
| `test_lambert.py` | 21 | 21 passed | Unit (independent Lambert cross-check layer) |
| `test_regression_published_grid.py` | 231 | 231 passed | Regression against archived snapshot |
| `test_rk4.py` | 16 | 16 passed | Unit (independent RK4 propagation layer) |
| `test_robustness.py` | 31 | 31 passed | Metamorphic / failure-injection |
| `test_surrogate.py` | 16 | 16 passed | Unit (surrogate reference values) |
| `test_surrogate_adversarial.py` | 20 | 20 passed | Adversarial / numerical stability |
| `test_surrogate_extra.py` | 43 | 43 passed | Unit (surrogate domain boundaries, package exposure) |
| `test_surrogate_properties.py` | 57 | 57 passed | Metamorphic (surrogate scaling laws) |
| **Total** | **585** | **585 passed, 0 skipped, 0 failed** | |

Verbatim console output of this run is archived alongside this document
as `pytest_full_run_output.txt`.

## 3. Regression test against the archived independent-comparison grid

`test_regression_published_grid.py` re-executes `otca_c.api.compute_transfer`
at all 230 (R, angular-deviation) pairs of the archived Run B grid
(`Independent_comparison/grid_compare_full_results.json`: R in
{1.02 ... 20}, deviation from Hohmann in {1 ... 30} deg, r1 = 6678 km,
mu = Earth) and compares the result against the `dv_core` value recorded
in that file. A second test checks the archived file's own internal
consistency (its claimed worst-case relative difference against the
independent method does not exceed the 5.7e-12 figure stated in
`Independent_comparison/FINAL_REPORT.md`).

This is a regression guard against unintended future changes to
`core.py` / `api.py`, distinct in purpose from the independent
cross-check already documented in `Independent_comparison/
FINAL_REPORT.md`: a failure here means the exact solver's own output has
changed since the snapshot was recorded, not that an independent method
disagrees with it.

Result: 231 of 231 cases passed (230 regression comparisons plus one
snapshot-consistency check), 0 failed.

## 4. Coverage of the two newer modules

`budget_inversion.py` (budget-to-angle inversion, relying on Theorem 1 for
existence/uniqueness and Theorem 3 for the safeguarded-Newton-in-`z`
solver's own structure) and `curvature_regime.py` (the R=2-only curvature
sign-change classifier relying on Theorem 4) each have their own dedicated
test file, `test_budget_inversion.py` (29 cases) and
`test_curvature_regime.py` (6 cases). The latter includes a pin test that
reads the certified constants `W_C_CERTIFIED` and `PHI_C_CERTIFIED_DEG`
directly out of `circular/theorem_4/03_central_mesh/mesh_certificate_R2.py`
and asserts they still match the module's own copy, so that a future
re-certification with a different value is caught here rather than
silently diverging. See `module_theorem_correspondence.md` for the full
dependency statement for both modules.

## 5. Run identification

- Suite executed: `pytest tests/ -v`, all 16 files under `tests/`, run
  from the repository root with `2i-OTCA-C-v1.0-demo/` and `circular/`
  present as siblings of `2i-OTCA-C-v1.0-algorithm/`.
- Package under test: `otca-c`, version `0.1.0` (`pyproject.toml`),
  installed in editable mode (`pip install -e .`).
- Result: 585 passed, 0 skipped, 0 failed.
- This document reports a single execution; it is not a statement about
  behaviour across Python versions, operating systems, or dependency
  versions beyond those installed at execution time (see ATBD Section
  9.7.1, Reproducibility, for the pinned/resolved dependency record).
