# Assumption H — Complete Elementary Proof

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Scientific status:** PROVEN — ANALYTICALLY, without restriction.
**Reproducibility status:** REPRODUCED (two independent scripts, see below).

## Statement

For every `R > 1` and every `phi in (0,pi)`, with `p` the semi-latus
rectum and `C(p,phi,R) = sqrt(S1(p,phi,R)) + sqrt(S2(p,phi,R))` the
two-impulse transfer cost (`S1, S2` as defined in the manuscript, Sec.
2–3):

- **(H1) Uniqueness.** The minimizer `p*(phi,R)` of `p -> C(p,phi,R)`
  is unique.
- **(H2) Non-degeneracy.** At that minimizer, `C_pp(p*(phi,R),phi,R) > 0`.

This proof in fact establishes something stronger than (H1)+(H2) taken
separately: `C_pp > 0` at *every* stationary point of `C(.,phi,R)`, not
only at an assumed minimizer. Combined with the boundary behaviour as
`p -> 0+` and `p -> infinity` (established in
`lemma1_closure_stationary_point.py`, which depends on this result),
this gives existence and uniqueness of the minimizer over the full
physical domain `p > 0`, for every `R > 1` and `phi in (0,pi)`.

## Notation

With `z := sqrt(p)`, `w := R^(-1/2) in (0,1)`, `q := cot(phi/2)`, and
`u := q/z > 0`:

```
a := 1 - w^2        b := 1 + w^2        A0 := 1 + w^4        m := 1 - w^3
```

`S1`, `S2` are the squared impulses at the two burns (manuscript Sec.
2–3); `C_zz` denotes the second derivative of `C` with respect to `z`
at fixed `(w,u)`; `C_pp = C_zz / (4p)` at a stationary point (direct
differentiation, since `dC/dp = (dC/dz)/(2z)` and the first-derivative
term vanishes at a stationary point).

## The proof, in outline

1. **Positive Laurent decomposition.** `S1+S2` and `K(z)/2 := S1'' =
   S2''` (both second derivatives equal: an equality established
   independently of this proof) both admit the form
   `d1*u^2 + d0 + d(-1)/u^2` and `c1*u^2 + c0 + c(-1)/u^2` respectively,
   with **every coefficient `c_i, d_i` strictly positive**. Verified
   exactly from the manuscript's original `S1`, `S2` formulas (Sec.
   2–3), not merely from an already-substituted intermediate form —
   see `assumption_H_independent_verification.py`, Step 1.

2. **Three-term Cauchy-Schwarz** gives
   `(K/2)(S1+S2) >= L^2`, with
   `L = sqrt(c1*d(-1)) + sqrt(c0*d0) + sqrt(c(-1)*d1) =: T1 + T0 + T3`.

3. **The two "tangential" terms (`T1+T3`)** are bounded by `a*sqrt(b)`
   via an elementary one-variable lemma,
   `sqrt(X^2+12) + |X-2| >= 2*sqrt(2X)`,
   proved by two exact factorizations, one per sign branch of `X-2`,
   each manifestly non-negative on its own domain.

4. **The "radial" term (`T0`)** is bounded by `w(1-w)/sqrt(2)` via an
   exact completion of the square.

5. **Final comparison**, purely algebraic:
   `a*sqrt(b) + w(1-w)/sqrt(2) > m` (with `m = 1-w^3`), established by
   an exact rationalization identity plus a chain of elementary
   inequalities (`1+w+w^2+(1+w)*sqrt(b) > 1+w+w^2 >= 2w > sqrt(2)*w`).

6. Hence `L > m`, so `(K/2)(S1+S2) > m^2`: **this is `Phi > 0`, for
   every `0 < w < 1`, `u > 0`, `z > 0` — no remaining hypothesis.**

7. Via `(U+V)^2 > S1+S2` (AM-GM, with `U=sqrt(S1)`, `V=sqrt(S2)`), this
   gives `C_zz > 0` **at every stationary point** — not only at the
   presumed minimizer: (H2) "at the minimizer" alone is not sufficient
   to conclude (H1), but `C_zz > 0` *everywhere `C_z = 0`* is.

8. Boundary behaviour (`C(z) -> +infinity` and `C_z -> -infinity` /
   `+const` at the two ends `z -> 0+, infinity`) plus a topological
   argument (two increasing crossings of `C_z = 0` would force a
   decreasing crossing between them, hence a point with `C_zz <= 0`,
   a contradiction) `=>` **the stationary point is unique**, and it is
   the global minimizer.

## Verification performed

Two independent scripts are provided, both re-executed in this session
with no error:

- `assumption_H_proof_certificate.py` — the external collaborator's
  own script, working directly in the substituted variables
  `(w,z,u)`. All `assert` statements pass.
- `assumption_H_independent_verification.py` — a from-scratch
  reconstruction of every algebraic step, starting from the
  manuscript's original `S1`, `S2` formulas in `p`, `R`, `phi`
  (Sec. 2–3) rather than from the collaborator's already-substituted
  notation. Every identity used in the proof (Steps 1–6 above) is
  recomputed independently: all differences are confirmed to be
  exactly zero, under the standing hypothesis `0 < w < 1` (guaranteed
  by `R > 1`).

Neither script uses interval arithmetic, SOS solvers, or a Sturm
sequence: every step is an exact symbolic identity or an elementary
one-variable inequality, each proved by an explicit factorization.

## Consequence for the manuscript

Assumption H (as stated above) is not a standing hypothesis anywhere
this proof applies: it is an unconditional result, for every `R > 1`
and every `phi in (0,pi)`. In particular, it is the result invoked by
`lemma1_closure_stationary_point.py` in this same folder, and, through
that chain, by `05_upper_bridge`'s certificate in `theorem_4/` (see
`assumption_H_certificate.json` for the exact dependency record).

---
*Frederic G. Speyser & Kseniya Vyatkina, with a contribution from an
external collaborator, independently verified.*
