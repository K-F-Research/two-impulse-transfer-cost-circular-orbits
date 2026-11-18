# Theorem 3 -- Universal Strict Convexity in the Reciprocal-Angle Variable

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

## What this folder contains

`theorem3_universal_convexity_proof.md` proves that the optimized
two-impulse transfer cost, re-parametrized in z := cot^2(phi/2) (phi being
the angular deviation from the Hohmann angle), is **strictly convex on the
entire physical domain, for every radius ratio R>1** -- not only near the
Hohmann angle, and not only for a specific ratio.

This is logically independent of, and stronger in scope than, the local
curvature result established for the single ratio R=2 in
`circular/theorem_4/`: that result concerns curvature in the
angle phi itself, which does change sign at R=2, while this one concerns
curvature in z, which does not change sign for any R>1. The two results
describe two different functions of the same physical transfer cost and do
not contradict one another.

The proof is an exact algebraic argument: it reduces convexity to the sign
of an explicit polynomial on the physical branch of a defining quartic,
closes most of the domain by an exact modular reduction and a
rational-arithmetic positivity certificate, and closes the remaining case
by isolating the real roots of a degree-1807 resultant with certified
ball-arithmetic (an exact method, not a floating-point heuristic) combined
with a topological continuity argument.

`verify_convexity_independent.py` is an independent numerical verification
verification script (see below).

## Status

- **Scientific status:** `PROVEN -- COMPUTATIONALLY CERTIFIED`. Every step
  of the proof rests on an exact computation (resultants, gcd, certified
  root isolation, a rational-arithmetic positivity certificate) or a
  standard topological continuity argument -- no step depends on a
  numerical approximation as its underlying justification.
- **Reproducibility status:** `PARTIALLY REPRODUCED`. See below.

## Reproducibility note

The proof document above is complete and self-contained: every step is an
exact algebraic argument (resultants, modular reduction, a rational-arithmetic
positivity certificate, certified root isolation, and a standard topological
continuity argument), and nothing in its logical chain is missing.

The raw polynomial certificate data behind Sec. 4's resultant chain is
included in this repository (`certificate_data/`) and has been used for a
genuine independent re-verification, not just read and cited: the four
exact identities `A*G==RH`, `B*G==RJ`, `gcd(A,B)==1`, and `gcd(RH,RJ)==G`
have been recomputed from the stored objects -- the last two via a
freshly-run `gcd` computation, not by re-reading a previously reported
value. All four passed. `gcd(A,B)==1` is exactly the coprimality fact the
proof document's Lemma 4.2 topological-continuity argument relies on. See
`certificate_data/README.md` for the full detail and for the one piece not
re-executed: the expensive upstream resultant computation that originally
produced `RJ` from the two source polynomials.

Two further independent checks, written from the proof document's own
definitions, directly confirm the theorem's content:

1. Symbolic and dense numerical re-confirmation that the leading
   coefficient of the defining quartic, F4(t,w), is strictly positive on
   the physical domain (2000 random sample points, zero exceptions).
2. A direct, independent numerical confirmation of the theorem's actual
   conclusion: the optimized transfer cost was minimized over the
   semi-latus rectum at 8 radius ratios and 400 angles each, re-parametrized
   to z, and checked for strict convexity via second divided differences
   on the resulting non-uniform grid (3184 interior checks, zero sign
   violations).

The status below remains `PARTIALLY REPRODUCED` rather than `REPRODUCED`,
because the single most expensive step of the original chain (the
resultant producing `RJ`) was not rerun from scratch -- not because the
certificate data was unavailable or unverified.

## Notation used in this folder

`R := r2/r1`, `t := sqrt(R)`, `phi := pi - Delta_theta` (consistent with
`docs/mathematical_conventions.md`). The variable `z := cot^2(phi/2)` is
specific to this folder and to `circular/theorem_4/`; see
`docs/mathematical_conventions.md` for the flagged collision with the
unrelated `z := sqrt(p)` used in `circular/existence_uniqueness_minimizer/`.

## Reproduction

```
python3 verify_convexity_independent.py
```

Requires `numpy`, `scipy`, and `sympy`. Expected output: zero F4
positivity violations over 2000 sample points, and zero convexity
violations over 3184 checked points.
