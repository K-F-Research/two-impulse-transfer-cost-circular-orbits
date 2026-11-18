# Module-to-Theorem Correspondence

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

This document provides the authoritative mapping between the results of
the manuscript, as well as the repository's own supplementary
(non-numbered) results, and the modules of `src/otca_c/` that rely on
them. Its purpose is to maintain this correspondence in a single,
explicit location, rather than relying on module filenames or
information dispersed across individual docstrings.

For the corresponding mapping for the `circular/` proof folders, see
`docs/script_article_correspondence.md` at the repository root. That
document covers the proof folders only and does not include the
`src/otca_c/` software package.

Each module in `src/otca_c/` includes a short "Theorems relied upon"
entry in its docstring, consistent with the dependencies recorded in the
table below. This document is therefore the authoritative reference for
auditing the consistency of these dependencies and should be updated
whenever a module's reliance on the manuscript or supplementary results
changes.

| Module | Theorems/results relied upon | Notes |
|---|---|---|
| `__init__.py` | None directly | Re-exports public names defined elsewhere only. |
| `errors.py` | None directly | Defines exception *types* only; the specific result enforced by a given exception is documented where that exception is raised. |
| `units.py` | None | Physical constants (mu, body radii) and unit conversions; standard astronomical data, not a result proven in this repository. |
| `core.py` | The repository's supplementary existence/uniqueness result (`circular/existence_uniqueness_minimizer/`) | Not a numbered lemma of the manuscript -- see `STATUS.md`. Establishes the sign pattern of `C_p` that makes the root-finder (`find_p_star`) safe. Theorem 1 and Theorem 2 are not needed: this module solves at a single, fixed `phi`. |
| `hohmann.py` | **Theorem 2** of the manuscript | Strict local minimum of the cost at the Hohmann transfer (`phi=0`). The proof folder `circular/hohmann_local_minimum/` files this same result internally as "Theorem 1" -- a documented numbering swap with the article (see `docs/mathematical_conventions.md`); this module's docstring and this table use the article's own number. |
| `surrogate.py` | **Theorem 2** of the manuscript | The closed-form local curvature coefficient `K(r1,r2)` at Hohmann (manuscript Section 4.6-4.8). `_K_normalized` reproduces that closed form term for term. |
| `geometry.py` | None | Independent cross-check via classical two-body relations (conic equation, angular momentum, vis-viva) -- deliberately not a re-derivation of this repository's own results. |
| `lambert.py` | None | Independent classical Lambert-problem solver (Vallado's universal-variable formulation), used only as a cross-check method. |
| `rk4.py` | None | Independent numerical propagation (fixed-step RK4) of the two-body equations of motion, used only as a cross-check method. |
| `derived.py` | None | Classical two-body relations (Kepler's equation, the Tsiolkovsky rocket equation) and altitude bookkeeping; independent of this repository's own proven results. Delegates to `lambert.py` for `lambert_cross_check`, itself theorem-independent. |
| `api.py` | The existence/uniqueness result (via `core.py`) and **Theorem 2** (via `hohmann.py`) | Theorem 1 (global monotonicity in `phi`) is not needed for the correctness of a single `compute_transfer` call -- it describes a separate, global property of how the cost varies across different `phi`, not the correctness of one solve. The evenness of the cost formula in `phi` (used to justify working with `phi_abs` for the cost/geometry solve) is an elementary algebraic property of the `S1`/`S2` primitives themselves (`core.py`'s `_e2`), not a conclusion of Theorem 2 or any other numbered theorem -- see `compute_transfer`'s inline comment. |
| `budget_inversion.py` | **Theorem 1** and **Theorem 3** of the manuscript | Theorem 1 (global monotonicity) establishes existence and uniqueness of the budget-to-angle inversion. Theorem 3 (strict convexity of `V_R` in `z := cot^2(phi/2)`, for every `R>1`) underlies the safeguarded-Newton-in-`z` solver's own structure. |
| `curvature_regime.py` | **Theorem 4** of the manuscript | Certified curvature sign change of `H(p,phi)` at `R=2` only. The module enforces `R=2` and raises `DomainError` otherwise; its extension to general `R>1` is Conjecture 2 (`circular/curvature_sign_change_general_R/`), an OPEN CONJECTURE, never extrapolated to here. |

Any change to a module's actual dependencies must update both that
module's own docstring and this table in the same change.
