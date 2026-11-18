# Theorem 2 -- Global Monotonicity of the Optimized Transfer Cost

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Numbering note.** This folder's name and the "Theorem 2" label above
follow this repository's own internal numbering. In the published
article, this same result is **Theorem 1**; the article's Theorem 2
(strict local minimum at Hohmann) is proved in
`circular/hohmann_local_minimum/theorem1_hohmann_strict_local_minimum.md`,
filed there under the same, swapped internal numbering. See
`docs/mathematical_conventions.md` for the full correspondence before
citing either theorem by number.

## What this folder contains

`theorem2_global_monotonicity_proof.md` is the complete, hand-reconstructible
proof that, for every fixed radius ratio R>1, the optimized two-impulse
transfer cost V_R(phi) := min_p F(p,phi,R) is **strictly increasing** in the
angular deviation phi from the Hohmann angle, over the entire open interval
(0,pi) -- not merely in a local neighborhood of the Hohmann angle (that
local result is Theorem 1 / `circular/hohmann_local_minimum/`).

The proof does not assume uniqueness of the minimizer p*(phi): it
establishes that every minimizer lies in a localization zone where the
partial derivative of the cost in phi is strictly positive (Sections 6-8),
then applies a convexity-free form of Danskin's theorem (Shapiro, Dentcheva
& Ruszczynski 2009, Theorem 7.21) to control the directional derivative of
the value function uniformly over any compact angular interval
(Sections 9-11), and concludes via an elementary real-analysis lemma (the
standard "sup technique", written out in full rather than merely cited).

`verify_monotonicity_independent.py` is an independent verification script
(see "Reproducibility note" below). `independent_numerical_checks/`
contains ten further independent numerical checks described there.

## Status

- **Scientific status:** `PROVEN -- ANALYTICALLY`. The proof is complete and
  self-contained; the two central algebraic identities have been
  independently recomputed from scratch, with no gap found.
- **Reproducibility status:** `REPRODUCED`. See the note below.

## Reproducibility note

An independent verification script, written directly from the proof
document's own Section 1 definitions (`verify_monotonicity_independent.py`):

1. Symbolically re-derives, with SymPy, the proof's two central closed-form
   identities -- `D'(p) = 2(1+c)/(p^3(1-c))` (Lemma 6.1) and the
   factorization `g(p) = A(p-p1*)(p-p2*)` (Lemma 8.3) -- directly from the
   Section 1 definitions, independently of the proof document's own
   derivation. Result: zero symbolic residual for both.
2. Numerically confirms global monotonicity itself, by directly minimizing
   `F(p,phi,R)` over `p` at 8 radius ratios (R = 1.2 to 50) and 60 angles
   each (0.5 deg to 179.5 deg), checking all 472 consecutive-angle pairs.
   Result: zero violations.

`independent_numerical_checks/` adds ten further independent numerical
checks: localization checks by Newton and bisection root-finding, a
bracketing confirmation, a reformulation via the classical Lagrange f,g
coefficients, double- and high-precision stress tests across extreme R and
boundary angles, a large-scale random sampling check, and an R<1 /
magnitude sanity check. See that folder's own `README.md` for details and
reproduction instructions.

## Notation used in this folder

Confirms the conventions in `docs/mathematical_conventions.md`:
`R := r2/r1` (Sec. 1, `r1=1`, `r2=R>1`) and `phi := pi - Delta_theta`
(Sec. 1, same convention as Theorem 1, Hohmann at `phi=0`).

## Reproduction

```
python3 verify_monotonicity_independent.py
```

Requires `numpy`, `scipy`, and `sympy`. Expected output: both symbolic
residuals equal to 0, and 0 monotonicity violations out of 472 checked
pairs.

See `independent_numerical_checks/README.md` for the additional numerical
validation scripts and how to reproduce them.
