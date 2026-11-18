# Lower-Endpoint Sign Conclusion

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Scientific status:** PROVEN — analytic argument (Taylor-Lagrange) combined
with an interval-certified second-derivative bound.
**Scope:** `theorem_4`, R=2, region `01_lower_endpoint`.

## 1. Definitions

Let `s := sqrt(p)`, and parameterize the stationary branch near the lower
endpoint by

```
x(s) = -1 + s^2 * Y(s).
```

Define

```
P(s) = s^6 * H(s^2, x(s)),
```

where

```
H(p,x) = ((1-x^2)*C_xx - x*C_x)*C_pp - (1-x^2)*C_px^2.
```

**Notation note.** The historical symbolic scripts referenced below use
`w` for this same quantity: `w = s = sqrt(p)`. This `w` is local to the
lower-endpoint construction and is consistent with the project-wide
convention `w := sqrt(p)` already documented in
`docs/mathematical_conventions.md` (shared with `inflection_uniqueness_R2`,
`universal_convexity_z`, `curvature_sign_change_general_R`). It is a
**different** quantity from the parameter `w` used in
`theorem_4/05_upper_bridge`'s certificate, where `w := R^(-1/2)`. This
document uses `s` throughout precisely to avoid that collision; `w` is
used only when quoting a script's own variable name.

## 2. Exact endpoint identities

At the endpoint,

```
s = 0,   Y(0) = 1/2.
```

Exact symbolic substitution gives

```
P(0) = 0.
```

Likewise, by direct substitution and symbolic differentiation (not by
first differentiating formally in four independent arguments and
recomposing by the chain rule — that second, independent method is used
separately for the second derivative, see Section 3),

```
P_s(0, 1/2) + P_Y(0, 1/2) * Y'(0) = 0,
```

i.e. `P'(0) = 0` along the branch. Both identities are exact symbolic
results (`sympy`, confirmed to simplify to exactly `0`), not floating-point
approximations.

## 3. Strict second-derivative bound

The interval certificate `certify_Ppp_v2.py` establishes

```
P''(s) < 0
```

throughout the certified rectangle

```
s in [0, 0.00015],   Y in [0.49987175, 0.500038475].
```

Stationary-branch inclusion in this rectangle is established separately,
within the same script, by the certified endpoint signs of `G` and the
strict inequality `G_Y < 0` throughout the rectangle:

```
G(s, Y_lo) > 0 for all s in the box,
G(s, Y_hi) < 0 for all s in the box,
G_Y(s, Y) < 0 strictly for all (s, Y) in the box.
```

By continuity and strict monotonicity of `G` in `Y`, this gives existence
and uniqueness of the root `Y(s)` of `G(s, .) = 0` inside `[Y_lo, Y_hi]`
for every `s` in the box — i.e. the true physical branch is contained in
the certified rectangle, not merely assumed to be. Consequently, the
second-derivative bound above applies to the entire stationary branch
segment under consideration, not only to an arbitrarily chosen box.

The value of `P''(0)/2` obtained independently along the branch (Section
5, `crosscheck_Ppp_methods.py`) is

```
P''(0)/2 = -9/32 - sqrt(2)/8 = -(9 + 4*sqrt(2))/32,
```

consistent with, and strictly inside, the certified sign bound above.

## 4. Taylor-Lagrange argument

For every `s in (0, 0.00015]`, Taylor-Lagrange's theorem gives a point
`xi_s in (0, s)` such that

```
P(s) = P(0) + s*P'(0) + (s^2/2)*P''(xi_s).
```

Using the exact endpoint identities of Section 2,

```
P(s) = (s^2/2) * P''(xi_s).
```

Since `xi_s` lies inside the certified rectangle of Section 3,
`P''(xi_s) < 0`, hence

```
P(s) < 0.
```

Since `s^6 > 0`,

```
H(s^2, x(s)) = P(s) / s^6 < 0.
```

Therefore,

```
H < 0 on the physical stationary branch, for 0 < s <= 0.00015.
```

## 5. Supporting artifacts

- `../scripts/build_P.py`, `../scripts/build_P2.py` — construct `P`
  (`F0 := s^4*H` after substitution `p=s^2, x=-1+s^2*Y`, then
  `P := s^2*F0` in cancelled form, giving `P = s^6*H`) and evaluate `P`
  exactly at the endpoint `(s,Y)=(0,1/2)`: `P(0) = 0`.
- `../scripts/verify_Pprime_symbolic.py` — direct symbolic verification
  of `P'(0) = 0` along the branch, by substituting `u1=sqrt(S1(s,Y))`,
  `u2=sqrt(S2(s,Y))` before differentiating (method A).
- `../scripts/build_Ppp_chainrule.py` — independent computation of the
  same second derivative via formal differentiation of `P` in its four
  arguments `(s,Y,u1,u2)` followed by the chain rule (method B).
- `../scripts/crosscheck_Ppp_methods.py` — confronts methods A and B
  term by term at `(s,Y,u1,u2) = (0, 1/2, sqrt(2), 1/sqrt(2))`: all four
  second partials agree exactly, and the recomposed `P''(0)/2` from
  method B matches `-(9+4*sqrt(2))/32` exactly.
- `../scripts/certify_Ppp_v2.py` — interval certification of the strict
  second-derivative bound and of stationary-branch inclusion in the
  rectangle of Section 3.
- `../scripts/F0_full_pipeline.py` — builds `H(p,x)` from the original
  two-impulse cost formulas and `F0 := s^4*H`; its own numeric evaluation
  of `F0` along the true branch `Y(s)` (not at a fixed `Y`) converges to
  `-(9+4*sqrt(2))/32` as `s -> 0`, consistent with, and independent
  corroboration of, the exact identity above.

The interval certificate and the symbolic identities are distinct
components of this argument: the symbolic identities establish `P(0)=0`
and `P'(0)=0` exactly, with no floating-point step; the interval
certificate establishes `P''<0` on a rectangle together with the
branch-inclusion property of that rectangle. Any numerical stress test
present elsewhere in this project is corroborative and is never used as
a substitute for either component.
