# Theorem 1 -- Strict Local Minimum of the Hohmann Transfer

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Numbering note.** This folder's name and the "Theorem 1" label above
follow this repository's own internal numbering. In the published
article, this same result is **Theorem 2**; the article's Theorem 1
(global monotonicity) is proved in
`circular/global_monotonicity/theorem2_global_monotonicity_proof.md`,
filed there under the same, swapped internal numbering. See
`docs/mathematical_conventions.md` for the full correspondence before
citing either theorem by number.

## What this folder contains

`theorem1_hohmann_strict_local_minimum.md` is the complete, hand-written
proof that the Hohmann transfer is not only the global cost minimum among
two-impulse transfers between two coplanar circular orbits (Prussing, 1992),
but also a **strict local minimum** in the transfer angle: the cost strictly
increases as soon as the transfer angle departs, even slightly, from the
Hohmann value. The core proof (Sections 1-9) is a pure analytic argument,
not a computation; `verify_proposition_4_3_subquadratic_scales.py` (see
below) is a later, separate numerical confirmation of one specific
generalization used within it (Proposition 4.3).

The proof is organized as five numbered lemmas plus one local
characterization (named "Lemma 3" for ease of reference), concluding with
the Theorem in Section 9: a closed-form second-order coefficient
$K(r_1,r_2)$ for the cost expansion around the Hohmann angle, shown positive
for every ratio $R=r_2/r_1 \neq 1$ via an explicit rationalizing
substitution and a difference-of-squares factorization into three manifestly
positive factors.

`verify_proposition_4_3_subquadratic_scales.py` independently confirms
Proposition 4.3's generalized form: that for every exponent
$\gamma\in[1,2)$, a path $p(\varphi)=p_H+\beta\varphi^\gamma$ is
non-competitive with the Hohmann-adjacent family, because its cost excess
decays like $\varphi^{2\gamma-2}$ -- strictly slower than $\varphi^2$. It
evaluates the full, unexpanded two-impulse cost function directly (no
series expansion) at 60-digit precision for $\gamma=1,\,5/4,\,3/2,\,7/4$
and a representative ratio $R=3$, and measures the empirical log-log decay
exponent as $\varphi\to0$.

`surrogate_benchmark/benchmark_quadratic_surrogate.py` (in its own
subfolder, with its own README) is a separate numerical benchmark
supporting Section 7.2's discussion of the curvature coefficient as a
cheap surrogate cost: it compares the quadratic surrogate
$V_R(0)+K(r_1,r_2)\varphi^2$ against the exact two-impulse cost (minimized
over $p$) on the constrained problem $\min V_R(\varphi)$ s.t.
$\varphi\ge\varphi_{\min}$, at $R=1.5$, for $\varphi_{\min}$ from $5^\circ$
to $30^\circ$. This constraint is necessary to make the comparison
non-trivial: by Theorem 1, the unconstrained problem always has its
minimum at $\varphi=0$.

## Status

- **Scientific status:** `PROVEN -- ANALYTICALLY`. The proof is complete and
  self-contained.
- **Reproducibility status:** `PARTIALLY APPLICABLE`. The core proof
  (Sections 1-9) is a hand proof with no associated script; every algebraic
  step in Sections 4bis through 9 and both appendices has been
  independently re-read against the stated conclusions, with no error
  found. Proposition 4.3's generalized exponent claim is additionally
  `REPRODUCED` by `verify_proposition_4_3_subquadratic_scales.py`: running
  it gives empirical exponents converging to $0,\,0.5,\,1,\,1.5$ for
  $\gamma=1,\,1.25,\,1.5,\,1.75$ respectively, matching the predicted
  $2\gamma-2$ in all four cases. `surrogate_benchmark/benchmark_quadratic_surrogate.py`
  is likewise `REPRODUCED`: it gives a surrogate relative error from
  $0.0003\%$ ($\varphi_{\min}=5^\circ$) to $0.24\%$ ($\varphi_{\min}=30^\circ$)
  at $R=1.5$, and a per-evaluation speedup on the order of several hundred
  times over the exact (numerically minimized) cost, on this reference
  implementation.

## Notation used in this folder

This folder confirms, rather than merely carries over, two of the entries in
`docs/mathematical_conventions.md`:

- `phi := pi - Delta_theta` is explicitly defined in Section 4 of the proof
  itself (not assumed from elsewhere), with the Hohmann transfer at
  `phi = 0` (i.e. `Delta_theta = pi`).
- `R := r2/r1` is explicitly used in Section 8 (`t := sqrt(R)`, with
  `r1 = 1`, `r2 = R > 1` after the Lemma 4/5 reduction).

Both conventions are now confirmed for this folder in
`docs/mathematical_conventions.md`, rather than left as an unconfirmed
carry-over note.

## Reproduction

The core proof has nothing to execute: to check it, read
`theorem1_hohmann_strict_local_minimum.md` directly; every step is intended
to be verifiable by hand, as stated throughout the document.
`verify_proposition_4_3_subquadratic_scales.py` requires only `mpmath` and
runs in well under a minute.
`surrogate_benchmark/benchmark_quadratic_surrogate.py` requires `numpy`
and `scipy` and runs in a few seconds; see its own subfolder README for
details.
