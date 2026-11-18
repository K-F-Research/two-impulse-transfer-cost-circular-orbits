# Theorem 4 (R=2) -- Certified Curvature Sign Change

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

## Statement

Let `C(p,phi)` be the two-impulse transfer cost at the radius ratio R=2,
and let the physical stationary branch be the curve `p*(phi)` defined by
`C_p(p,phi) = 0` for `0 < phi < pi`. Define the curvature determinant

    H(p,phi) = C_phiphi(p,phi) * C_pp(p,phi) - C_pphi(p,phi)^2.

Theorem 4 states that `H` changes sign **exactly once** along the
physical branch: there is a unique point `w_c` (`w := sqrt(p)`) such that
`H < 0` for `w < w_c` and `H > 0` for `w > w_c`, with

    w_c = 0.30698081390187630871...,  p_c = 0.09423722010385841476...,
    phi_c = 2.86881604632355948264...  (radians, approx. 164.371 degrees).

These decimal values are given for orientation only; the proof relies on
the certified enclosures stored in the certificate file for each region.

As stated in `docs/mathematical_conventions.md` and in the manuscript
(Section 4.3), this result does not contradict Theorem 3's universal
strict convexity of `V(z)` in the reparametrized variable
`z := cot^2(phi/2)`: the two statements concern the sign of the second
derivative of the same cost function with respect to two different
horizontal variables, and a monotonic but non-affine change of variable
can turn a one-way bend in `z` into a sign change in `phi`.

## Proof architecture

The physical domain `0 < w < w_H = sqrt(4/3)` is covered by six regions,
each with its own certificate:

| Region | Scope | Method |
|---|---|---|
| `01_lower_endpoint/` | `0 < s <= 0.00015` (`s := sqrt(p)`) | exact endpoint identities and a Taylor-Lagrange bound, `H < 0` |
| `02_lower_bridge/` | `0.00015 <= s <= 0.15` | certified adaptive interval mesh, `H < 0` |
| `03_central_mesh/` | `w in [0.15, 1.144701]` | Krawczyk-certified isolation of `w_c`, flanked by a certified interval mesh in both directions |
| `04_hohmann_endpoint/` | `0 < t <= t_archive` (`t := sin(phi/2)^2`, near the Hohmann point `p = 4/3`) | regularized branch map, certified sign of the rescaled quantity `t*H`, `H > 0` |
| `05_upper_bridge/` | `t_archive <= t <= 0.19055...` | certified adaptive interval mesh in `t`, `H > 0` |
| `06_global_audit/` | cross-cutting | sign-convention and transformation checks spanning all five regions above, plus an independent re-certification of one region's archived data |

Each region uses the domain parameter most convenient for its own
construction (`s`, `w`, or `t`). The exact maps between these parameters,
as well as the junction certificates establishing that the five
sign-determined regions form a single continuous branch, are documented
in each region's README and in `THEOREM4_PROOF_MAP.md`.

## Reproducing the computations

Each region's README lists the exact commands and dependencies for its
scripts. The regions depend on one another only through the shared
pickled symbolic objects and junction values documented in their READMEs
(for example, `01_lower_endpoint/build_P.py` must run before
`02_lower_bridge/full_mesh_bridge.py`). There is no single top-level
driver script, since the six regions were certified independently and
then joined at their shared boundaries.
