# Theorem 4 (R=2) — Region 02_lower_bridge: sign argument

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

## 1. Scope and notation

This document establishes `H < 0` (strict negativity of the curvature
determinant) on the closed segment `s in [0.00015, 0.15]`, where
`s := sqrt(p)`, continuing the stationary-branch segment covered by
`theorem_4/01_lower_endpoint` (`0 < s <= 0.00015`) up to the point
`s = 0.15`.

`s` is used here exclusively; the certifying scripts use the historical
variable name `w` for the same quantity (`w = s`), consistent with the
project-wide convention `w := sqrt(p)` documented in
`docs/mathematical_conventions.md`. As in `01_lower_endpoint`, this is
a different quantity from the `w := R^(-1/2)` used in
`theorem_4/05_upper_bridge`'s certificate.

The branch parameter `x(s) = -1 + s^2*Y(s)` and the function
`P(s) := s^6*H(s^2, x(s))` are defined exactly as in
`01_lower_endpoint/proof/lower_endpoint_sign_conclusion.md`, Sections 1
and 5, and reuse the same pickled symbolic form of `P`
(`P_built.pkl`), so the two regions certify the same underlying
quantity.

## 2. Method: direct interval bound, away from s=0

Unlike `01_lower_endpoint`, where `s` approaches `0` and the sign of
`H` had to be obtained from the exact identities `P(0)=0`, `P'(0)=0`
together with an interval-certified bound `P''<0` via a
Taylor-Lagrange argument, the segment `[0.00015, 0.15]` is bounded away
from `0`. There, `H<0` is established directly, tile by tile, by a
first-order interval (mean-value) enclosure of `P` itself:

For a tile `[a,b]` with midpoint `w0` and half-width `h = (b-a)/2`:

  1. The stationary branch `Y(w)` (the root in `Y` of the defining
     equation `G(w,Y)=0`, `G := w^2 * C_p` in the regularized
     variables) is located at `w0` by root-finding, then certified to
     lie inside an interval box `Y_B = [Y_lo, Y_hi]` by the same
     sign-bracketing-plus-monotonicity argument as in
     `01_lower_endpoint/scripts/certify_Ppp_v2.py`: `G_Y<0` strictly on
     the box `[a,b] x Y_B`, and `G` changes sign strictly between
     `Y_lo` and `Y_hi` for every `w` in `[a,b]`.
  2. `Q0`, the interval value of `P` at the tile's center `(w0, Y_0)`
     (`Y_0` a certified box around the root at `w0` alone), is
     computed.
  3. `T`, the interval value of the total derivative of `P` along the
     branch (`P_w + P_Y * Y'(w)`, `Y'(w) = -G_w/G_Y`), is computed on
     the full tile box `[a,b] x Y_B`.
  4. The mean-value enclosure `Q_final := Q0 + T * [-h, h]` bounds
     `P` on the whole tile `[a,b]`.
  5. The tile is accepted only if `margin := -sup(Q_final) > MARGIN_TOL`
     (`MARGIN_TOL = 1e-12`), i.e. `P < 0` strictly throughout the tile,
     together with four further safeguards that must hold
     simultaneously: `S1>0` and `S2>0` on the tile (well-definedness of
     the square roots in the construction), the denominator of `G`
     separated from `0` by at least `DENG_SEP_MIN = 0.5`, and
     `width(Y_B) <= YB_WIDTH_MAX = 0.05`.

If any criterion fails, the tile is subdivided (never relaxed) and the
process repeats on the two halves. Since `s = w > 0` throughout this
region (`s >= 0.00015`), `s^6 > 0`, so `P<0` on a tile gives
`H(s^2,x(s)) = P(s)/s^6 < 0` on that same tile directly -- no
Taylor-Lagrange step is needed here, because `s` never approaches `0`.

## 3. Mesh construction and coverage

Starting from `s = 0.00015`, the domain is covered by tiles of initial
relative half-width `0.75%` of the current point, growing geometrically
with `s`; any tile failing the criteria above is bisected (up to a
maximum subdivision depth of 40, never reached). The resulting mesh
covers `[0.00015, 0.15]` exactly, contiguously, with no gap and no
overlap: 464 tiles, 0 unrecoverable failure.

## 4. Tile-to-tile stitching

Beyond certifying each tile in isolation, the construction verifies
that the stationary branch itself is continuous across tile
boundaries: at each shared boundary `w_s`, the hull (union) of the two
adjacent tiles' `Y_B` boxes is re-checked for `G_Y<0` and for the
correct sign change of `G` at `w_s` on that hull. Geometric overlap of
the two boxes alone would not establish this; the hull re-check does.
463 internal boundaries were tested, 0 failure.

## 5. Independent audit (two levels)

Two independently-written programs audit the mesh archive without
sharing code with the generator beyond the common symbolic
construction of `S1`, `S2`, `P`, `G`:

  - `scripts/pilot_mesh_audit.py` ("level 1"): recomputes nothing
    numerically. It re-reads the archived per-tile bounds and the
    archived stitching records, and re-checks every criterion
    (coverage, `G_Y<0`, `S1>0`, `S2>0`, denominator separation, sign
    bracketing at the branch box edges, `Y_B` width, margin) directly
    from those archived bounds.
  - `scripts/pilot_mesh_audit_level2.py` ("level 2"): reconstructs the
    symbolic chain independently and, for every tile, re-evaluates by
    interval arithmetic from the archived `W`, `Y0`, `Y_B` boxes alone
    (not from the archived derived quantities), checking that the
    re-evaluated values satisfy every criterion and match the archived
    values within serialization tolerance (`1e-30`).

Both audits pass on the archive at `archive/lower_bridge_archive.json`:
level 1 and level 2 each report 464/464 tiles and 463/463 stitchings
passing, with 0 violation and 0 archive/re-evaluation discrepancy.

## 6. Junction with 01_lower_endpoint

`01_lower_endpoint`'s certified domain is `0 < s <= 0.00015` (closed at
the upper end). This region's certified domain is `0.00015 <= s <=
0.15` (closed at the lower end), and the first tile's `w_lo` equals
`0.00015` exactly. The two regions therefore share the single point
`s = 0.00015` under an identical definition of `P` and `H` (the same
`P_built.pkl`), with `H<0` established on each side by its own method
(Taylor-Lagrange approaching `0` on one side, direct interval
mean-value bound away from `0` on the other). No separate branch-
identification argument is needed at this junction, since both
constructions certify the same stationary branch of the same function
on a domain that meets at a single shared, doubly-covered point.

## 7. Conclusion

`H(s^2, x(s)) < 0` for every `s in [0.00015, 0.15]`, certified tile by
tile with explicit sign margins, continuous branch stitching across
all 463 internal boundaries, and independently re-verified at two
separate audit levels (Section 5).
