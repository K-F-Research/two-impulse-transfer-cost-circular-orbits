# Theorem 4 (R=2) — Region 03_central_mesh: sign-change argument

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

## 1. Scope and notation

This document establishes, on the stationary branch of the two-impulse
cost function `C(p,phi) = sqrt(S1) + sqrt(S2)` for `R=2`, that the
curvature determinant

```
H(p,phi) := C_phiphi * C_pp - C_pphi^2
```

changes sign **exactly once** over `w in [0.15, 1.144701]`
(`w := sqrt(p)`), with `H < 0` immediately below the sign-change point
and `H > 0` immediately above it. This continues the `H < 0` result
already certified on `[0.00015, 0.15]` (`02_lower_bridge`) and is
itself continued, on the Hohmann side, by `05_upper_bridge` and
`04_hohmann_endpoint`.

This region works directly in the native variables `(p, phi)` of the
manuscript's two-impulse cost model, not in the regularized `(w, Y)`
substitution used in `01_lower_endpoint` and `02_lower_bridge`: away
from `w=0`, the raw formulas for `S1`, `S2`, `C`, and their derivatives
are already well-conditioned for interval evaluation, so no
substitution is needed here.

## 2. Locating the sign-change point: a Krawczyk-certified interior zero

A point `K_1 = (p_c, phi_c)` on the stationary branch (`Cp(p_c,phi_c)=0`)
is located by ordinary high-precision Newton iteration, then certified
as a genuine, isolated zero of the map `(p,phi) -> (Cp, H)` by a
Krawczyk interval test: with `J0` the Jacobian of this map at the
Newton-refined point and `Y := J0^{-1}` the preconditioner, the interval
evaluation of the Krawczyk operator

```
K(box) = (z0 - Y*F(z0)) + (I - Y*J(box)) * (box - z0)
```

on a box of half-width `eps = 1e-7` around `(p_c, phi_c)` is checked to
lie strictly inside that same box, certifying existence and uniqueness
of the zero within it. The preconditioner's own non-singularity is
certified by interval arithmetic (`det(Y)` strictly separated from 0),
not only by the floating-point diagnostic `det(J0)`.

## 3. Mesh construction and per-tile certification

Starting from the two edges of the certified `K_1` box, the branch is
followed outward in both directions (DOWN toward `w=0.15`, UP toward
`w=1.144701`) by an adaptive walk. Each tile `[w_lo, w_hi]` is
certified by:

  1. An interval-Newton contraction locating and certifying a unique
     root `phi(p)` of `Cp(p,.)=0` on the tile's `p`-box, via a
     mean-value enclosure of the numerator (avoiding the interval-
     dependency blow-up of the raw, deeply nested expression for `Cp`).
  2. An explicit check that `S1>0` and `S2>0` on the full tile box
     (well-definedness of the square roots in `C`).
  3. A first-order (mean-value) Taylor enclosure of `H` at the tile's
     center, with a certified interval remainder bound; the tile is
     accepted only if this remainder is strictly smaller in magnitude
     than the certified value of `H` at the center, which both bounds
     the tile and certifies the sign of `H` throughout it.
  4. A tile-to-tile continuity check: each tile already certifies
     existence and uniqueness of the branch root over its entire
     `p`-block, including the point shared with its neighbor; a
     nonempty intersection of the two tiles' certified `phi`-enclosures
     at that shared point therefore forces, by uniqueness, that both
     tiles certify the same root, hence the same branch.

If a tile's criteria are not met, its width is rescaled (never
relaxed) and recertified.

## 4. Results

- **DOWN walk** (`K_1` to `w=0.15`): 2052 tiles, constant sign `H<0`,
  2051 tile-to-tile boundary checks, all passed.
- **UP walk** (`K_1` to `w=1.144701`): 3090 tiles, constant sign `H>0`,
  3089 tile-to-tile boundary checks, all passed.
- **S1, S2 positivity**: checked explicitly on every one of the 5142
  tiles; 5460 checks passed (some tiles are checked more than once
  during adaptive rescaling), 0 failed.
- **Krawczyk preconditioner**: `|det(J0)| ~ 2.07e4` (floating-point
  diagnostic); `det(Y)` interval-certified strictly separated from 0
  (rigorous non-singularity certificate for the preconditioner
  actually used).

Since the DOWN walk's sign set is `{-1}` and the UP walk's is `{+1}`,
and the two walks originate from a single Krawczyk-certified isolated
root, `H` changes sign exactly once on `[0.15, 1.144701]`, at the
certified point `K_1` (`w_c approx 0.306980813901876...`).

## 5. Method: interval arithmetic throughout

Every certified quantity in this region -- the branch root enclosures,
the `S1`/`S2` positivity checks, the sign of `H` on each tile, the
Krawczyk non-singularity certificate -- is obtained by interval
arithmetic (`mpmath mp.iv`, working precision 50 significant digits,
with extractions to plain-precision values performed only at a working
precision of 60 digits, strictly above the interval precision, so that
every such extraction is exact). No step in this region relies on an
exact symbolic identity (unlike the endpoint identities `P(0)=0`,
`P'(0)=0` used in `01_lower_endpoint`) or on a floating-point stress
test presented as a proof.

## 6. Junction with 02_lower_bridge

`02_lower_bridge`'s certified domain is `[0.00015, 0.15]` (`H<0`
throughout, closed at both ends). This region's DOWN walk starts at the
`K_1` edge and is certified down to exactly `w=0.15`, with `H<0`
throughout. The two regions meet at `w=0.15` under the same definition
of `H(p,phi)` (native variables, no substitution needed on either side
of this point), with `H<0` certified by each region's own method on
its own side of that shared point.

## 7. Conclusion

`H(p,phi)` is strictly negative on the stationary branch for
`w in [0.15, w_c)`, strictly positive for `w in (w_c, 1.144701]`, and
changes sign exactly once, at the Krawczyk-certified point `K_1`
(`w_c`). Combined with `02_lower_bridge`, `H<0` is therefore certified
on the entire segment `[0.00015, w_c)`, with the sign change itself
certified rigorously rather than merely observed numerically.
