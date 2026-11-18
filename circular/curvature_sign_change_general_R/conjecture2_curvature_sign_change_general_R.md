# Conjecture — Single Curvature Sign Change of the Optimized Transfer Cost for General R

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

## Statement

Let `V_R(phi) := min_p C(p,phi)` be the optimized two-impulse transfer cost for
radius ratio `R>1`, as in `circular/global_monotonicity/` and
`circular/theorem_4/`, and let `w := sqrt(p)` on the physical
stationary branch `C_p(p,phi) = 0`. Define `P(w,R)` as the explicit algebraic
quantity (a polynomial combination of the problem's primitive radicands,
denoted `P` and reduced from `Delta_H` by the exact identity `Delta_H = P / S1^3`)
whose sign on the branch coincides with the sign of `V_R''(phi)`, exactly as in
`circular/theorem_4/`, Section 1, but with `R` now left free.

**Conjecture.** *For every `R>1`, `P(w,R)` -- equivalently `V_R''(phi)` -- changes
sign exactly once on the physical interval `w in (0, w_H(R))` (equivalently
`phi in (0,pi)`).*

For `R=2` this is exactly Theorem 4 (`circular/theorem_4/`), which is fully
certified: six certified regions, joined at four certified junctions, cover
the entire physical interval with no gap (see
`circular/theorem_4/06_global_audit/global_coverage_summary.md`). For
general `R`, the statement is at the stage of a conjecture: strongly
corroborated numerically, but not established by either an exact algebraic
argument or a certified computation.

## What is established

- **The reduction itself is unconditional and does not depend on this
  conjecture.** The identity `Delta_H = P/S1^3` (hence the reduction of the
  curvature-sign question to the sign of `P`) holds for every `R>1`, exactly
  as it does for `R=2`; it was proven without reference to any particular
  value of `R` and without assuming the result below.
- **The existence and non-degeneracy of the minimizer `p*(phi,R)`** on which
  the whole construction (`V_R`, the branch `w(phi,R)`, `P`) rests is itself
  now an unconditional theorem for every `R>1` and every `phi in (0,pi)` (an
  auxiliary technical hypothesis used in earlier drafts of this material, and
  historically labeled as an assumption, was closed by an exact three-term
  Cauchy-Schwarz argument combined with a one-variable lemma; it no longer
  restricts the scope of this conjecture or of any other result in this
  repository).
- **Dense numerical corroboration.** Direct computation of `P(w,R)` along the
  physical branch, for a wide range of `R>1`, shows a single sign change in
  every case tested, consistent with the conjecture; no counterexample --
  no second sign change, and no value of `R` with zero or more than one sign
  change -- has been found at any point during this project's work.
- **R=2 is a fully worked model case**, not merely a data point: its proof
  (Krawczyk certificate at the known critical point, a certified interval
  mesh over a large explicit sub-interval, and two exact endpoint limits;
  see `circular/theorem_4/`) shows concretely what closing
  the conjecture for a given `R` requires.

## Why this is not a mechanical extension of the R=2 proof

The R=2 proof is not a generic argument instantiated at `R=2`: several of its
key ingredients are numerical constants or certified objects computed
specifically for that ratio --

- the exact endpoint limit values (Section 5 of the Theorem 4 proof document)
  are closed-form numbers obtained from a Laurent-Puiseux expansion whose
  coefficients depend on `R` and were only carried out at `R=2`;
- the certified interval mesh of Section 4 there (covering a large explicit
  sub-interval of `w`) is a numerical certification tied to the specific
  endpoints `w_H(2)` and the specific critical point `w_c(2)`;
- the Krawczyk box of Section 3 there is centered on the specific numerical
  root for `R=2`.

None of these carries over to another value of `R` without being recomputed,
and a literal repetition of the same computation at a dense grid of `R`
values would only ever produce evidence for finitely many ratios, never a
proof for all `R>1` simultaneously. Closing the conjecture in general requires
either an argument that is uniform in `R` from the start (most plausibly via
the same `t := sqrt(R)`, `w := sqrt(p)`, `z := cot^2(phi/2)` reduction used to
close Theorem 3's universal convexity result, which does treat all `R>1` at
once), or a demonstrably `R`-uniform version of the certified-computation
route used for `R=2`.

## Recommended approach, not yet carried out

The most promising path identified so far reduces the conjecture to a single
transversality lemma in the `(t,w,z)` variables already used throughout this
repository's `circular/` results: showing that the relevant defining
quantity crosses zero transversally (nonzero derivative) at its unique
physical zero, uniformly in `t>1`. This would mirror, in structure, the
exact-algebra closure already achieved for Theorem 3's universal convexity
result, which successfully handled the analogous "uniform in `R`" difficulty
for a related but logically independent statement. This path has been
identified but not executed: no uniform transversality argument has yet been
written out or verified for this conjecture.

## Status

- **Scientific status:** `OPEN CONJECTURE`. Proven for `R=2`
  (`circular/theorem_4/`, fully certified with no gap in domain coverage);
  not established for general `R>1` by any exact or certified argument.
- **Reproducibility status:** `NOT APPLICABLE` -- an open conjecture has no
  claimed computation to reproduce. The numerical corroboration described
  above is exploratory evidence, not a certified result, and is not claimed
  as such.

## Notation

`R`, `phi`, `w := sqrt(p)`, `t := sqrt(R)`, `z := cot^2(phi/2)` -- all
consistent with `docs/mathematical_conventions.md`. `P` here denotes the same
quantity (up to the stated exact factor `S1^3`) as `H`'s numerator in
`circular/theorem_4/`, generalized to free `R`; this is a
notational continuation, not a new collision.
