# Proof Dependencies

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

This document maps, for each result in `circular/`, the hypotheses it
relies on and the other results it depends on. See
`mathematical_conventions.md` for the article-vs-folder numbering
correspondence before reading theorem numbers below, and
`script_article_correspondence.md` for script-level notation.

## Independent results (no dependency on another `circular/` result)

- **Article Theorem 1 -- global monotonicity**
  (`circular/global_monotonicity/`, filed internally as "Theorem 2"). Proved
  directly from the cost function's own definition; explicitly does **not**
  assume uniqueness of the minimizer `p*(phi)` (see the proof document's
  Section 3 and this folder's own README).
- **Article Theorem 2 -- strict local minimum at Hohmann**
  (`circular/hohmann_local_minimum/`, filed internally as "Theorem 1"). A
  self-contained analytic argument (five lemmas plus one local
  characterization), independent of every other folder in `circular/`.
- **Article Theorem 3 -- universal convexity in `z`**
  (`circular/universal_convexity_z/`). An independent algebraic argument
  reducing convexity of `V(z)` to the sign of an explicit polynomial on the
  physical branch of a defining quartic; does not rely on Theorem 1, 2, or 4.
  Explicitly does not contradict Theorem 4 (see below): the two results
  concern the sign of the second derivative of the same cost function under
  two different choices of horizontal variable (`z` vs. `phi`).

## Supplementary result: existence and uniqueness of the minimizer

`circular/existence_uniqueness_minimizer/` proves that, for every `R>1` and
every `phi` in `(0,pi)`, the cost function admits a unique global minimizer
`p*(phi,R)`. It is **not** a numbered result of the published article: the
article's own Theorem 1 (global monotonicity) explicitly does not require
it, by a stronger, more direct argument (see above). This folder's result
is kept because it is true, proved, and still useful elsewhere in this
repository:

- It is the stated justification for the numerical solver (`otca_c.core`)
  returning the actual global minimum via root-finding, rather than an
  arbitrary local one.
- It is a direct dependency of **article Theorem 4**'s region
  `circular/theorem_4/05_upper_bridge/`: that region's certificate relies on
  global uniqueness of the stationary root at fixed `phi` to establish that
  the certified Hohmann tile family represents the entire physical branch
  segment it covers (see `05_upper_bridge/upper_bridge_certificate.json`,
  field `branch_identification_note.closure`, and
  `existence_uniqueness_minimizer/README.md`).

## Article Theorem 4 -- certified curvature sign change (R=2)

`circular/theorem_4/` proves that the curvature determinant `H(p,phi)`
changes sign exactly once along the physical stationary branch `p*(phi)`
at the single ratio `R=2`. It depends on:

- the existence of the physical branch `p*(phi)` and, specifically in the
  region `05_upper_bridge/`, on its **uniqueness** (see
  `existence_uniqueness_minimizer/`, above);
- no dependency on Theorem 1, 2, or 3.

Internally, the proof is a region-by-region decomposition of the physical
domain into six independently certified regions (`01_lower_endpoint/`
through `06_global_audit/`), joined at their shared boundaries by explicit
junction certificates; see `circular/theorem_4/README.md` and
`THEOREM4_PROOF_MAP.md` for the exact decomposition and junction structure.

## Conjecture 2 -- open generalization of Theorem 4 to general R

`circular/curvature_sign_change_general_R/` states, as this repository's
own formalization of a one-sentence remark in the article's Section 7.3,
the open generalization of Theorem 4 from the single ratio `R=2` to every
`R>1`. It is not a numbered conjecture of the published article. It
depends on:

- **Theorem 4** (`circular/theorem_4/`) as the fully worked `R=2` model
  case that the open generalization targets;
- the reduction framework of **Theorem 3** (`circular/universal_convexity_z/`):
  this folder explicitly reuses the full `t, w, z` reduction from that
  folder, together with the `w` variable of `theorem_4`;
- the exact identity `Delta_H = P/S1^3`, which is itself unconditional and
  already proven, and on which the open reduction rests.

This folder is an honest record of an open problem, not a proof: nothing in
it is claimed as proven or certified (see its own README for the full
scope statement).
