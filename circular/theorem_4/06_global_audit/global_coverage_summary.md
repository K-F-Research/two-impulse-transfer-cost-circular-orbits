# Theorem 4 (R=2) — Global coverage summary

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

## Purpose

This document assembles the five regional certificates of `theorem_4`
(`01_lower_endpoint` through `05_upper_bridge`) into a single statement
of the domain each one covers, the sign of `H` it certifies there, and
the junction with its neighbor(s). It introduces no new computation: every
claim below is cited from the corresponding region's own certificate.

## Statement

Let `C(p,phi) = sqrt(S1) + sqrt(S2)` be the two-impulse cost function
for `R=2`, let the physical stationary branch be `C_p(p,phi)=0`, and
let

```
H(p,phi) := C_phiphi(p,phi) * C_pp(p,phi) - C_pphi(p,phi)^2
```

be the curvature determinant along that branch. Using `w := sqrt(p)`
(equivalently `s` in some regions' own notation; see each region's
README for the exact correspondence) and `t := sin(phi/2)^2` in the
two regions nearest the Hohmann point, the domain of interest is the
open interval `w in (0, w_H)`, `w_H = 2/sqrt(3)`.

## Region-by-region coverage

| Region | Domain | Sign of H | Method |
|---|---|---|---|
| `01_lower_endpoint` | `0 < w <= 0.00015` | `H < 0` | exact symbolic identities `P(0)=0`, `P'(0)=0`, plus an interval-certified bound `P''<0`, combined via Taylor-Lagrange |
| `02_lower_bridge` | `0.00015 <= w <= 0.15` | `H < 0` | interval mean-value enclosure of `P`, 464 tiles, 463 tile-to-tile stitchings, independently re-evaluated at two separate audit levels |
| `03_central_mesh` | `0.15 <= w <= 1.144701` | `H < 0` for `w < w_c`, `H=0` at `w_c`, `H > 0` for `w > w_c` | interval mean-value enclosure of `H` around a Krawczyk-certified isolated zero `K_1` (`w_c approx 0.306980813901876...`), 5142 tiles total, 5140 tile-to-tile continuity checks |
| `05_upper_bridge` | `1e-4 <= t <= 0.19055605997847746...` (`t := sin(phi/2)^2`) | `H > 0` | adaptive interval mesh in `t`, 98 tiles, 97 internal junctions, 0 failure |
| `04_hohmann_endpoint` | `0 < t <= 1e-4` | `H > 0` | subdivided interval bound on the rescaled quantity `t*H`, `n_a=100` sub-boxes, certified lower bound `~1.061 > 0` |

## Junctions between adjacent regions

| Shared boundary | Regions | Certifying script | Result |
|---|---|---|---|
| `w = 0.00015` | `01_lower_endpoint` / `02_lower_bridge` | documented in `02_lower_bridge`'s certificate (`junction_with_01_lower_endpoint`) | the two domains meet at a single closed point under the same definition of `P` and `H` (shared `P_built.pkl`) |
| `w = 0.15` | `02_lower_bridge` / `03_central_mesh` | documented in `03_central_mesh`'s certificate (`junction_with_02_lower_bridge`) | the two domains meet at a single shared point; `H<0` certified on each side by its own method up to and including that point |
| `w = 1.144701` | `03_central_mesh` / `05_upper_bridge` | `hohmann_central_join.py` (certified in `05_upper_bridge`'s certificate, junction `J_UPPER_BRIDGE_CENTRAL_MESH`) | non-empty overlap of the `p`-boxes and `phi`-boxes from each side's last tile, strict sign of `C_pphi` on the overlap (uniqueness of `phi` at fixed `p`), existence certified on both sides for every `p` in the overlap: same physical `(p,phi)` point on both sides |
| `t = 1e-4` | `05_upper_bridge` / `04_hohmann_endpoint` | `hohmann_raccord_t1e4.py` (certified in both regions' certificates, junction `J_HOHMANN_LOCAL_UPPER_BRIDGE` / `J_HOHMANN_UPPER_BRIDGE`) | exact rational equality of the shared boundary value, strict inclusion of the mesh's first-tile `a`-box inside the local certificate's `a`-box, identity of the frozen symbolic expressions on both sides (shared SHA-256), independent re-verification of existence and uniqueness on both sides |

Each junction above establishes a branch-identity obligation (the two
constructions on either side certify the *same* physical stationary
point, not merely overlapping numerical ranges) in addition to the
sign obligation already certified separately on each region's own
tiles.

## Resulting global picture

Combining the five regions and four junctions above, the sign of `H`
on the open interval `w in (0, w_H)` is:

```
H(w) < 0   for  0 < w < w_c
H(w_c) = 0
H(w) > 0   for  w_c < w < w_H
```

with `w_c approx 0.306980813901876...` the unique, Krawczyk-certified
sign-change point located in `03_central_mesh`. No gap in `w`-coverage
exists between `w=0` and `w=w_H`: every point in the open interval
belongs to exactly one region's certified domain, or to a doubly-
covered junction point where both sides agree.

## What this document does not establish

This summary does not itself recompute, re-derive, or independently
audit any of the cited results; it is an assembly of claims already
certified in each region's own certificate file, which remain the
primary source. In particular:

- No region's certificate has yet received the status `auditor.status:
  "verified"` from a reviewer independent of the certificates'
  authors; every certificate currently records `auditor.status:
  "pending"`.
- This document does not perform the sign-convention or coordinate-
  transformation checks described in `theorem_4/README.md`'s own
  outline of `06_global_audit`'s intended scope; it covers only the
  domain-coverage and junction-consistency assembly.

The result `theorem_4` establishes should accordingly be read as
certified by the computations cited above, pending independent review
of those computations, not yet as independently audited.
