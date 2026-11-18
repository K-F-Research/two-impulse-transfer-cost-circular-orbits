# Open Question -- Single Curvature Sign Change for General R

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Not a numbered conjecture of the published manuscript.** The published
manuscript (Reference [8]) does not state this as a numbered conjecture.
Section 7.3 ("Limitations and future work") mentions it only in a single
sentence: "Extending Theorem 4 to general R > 1 and closing its remaining
endpoint gap are left for future work." This folder is this repository's own
precise formalization of that one-sentence remark, not a restatement of an
article-numbered result.

## What this folder contains

`conjecture2_curvature_sign_change_general_R.md` states, precisely, the open
generalization of Theorem 4 (`circular/theorem_4/`) from the
single ratio `R=2` to every `R>1`: that the curvature of the optimized
transfer cost `V_R(phi)` changes sign exactly once on `(0,pi)`, for every
`R`, not only `R=2`.

This folder is an honest record of an open problem, not a proof. It exists
in the public repository because the reduction it rests on (the exact
identity `Delta_H = P/S1^3`) is itself unconditional and already proven, and
because the `R=2` case is a fully worked model of what a complete proof
would require -- so this formalized target is well-posed and precisely
stated, not a vague hope, even though the article itself does not number it.

## Status

- **Scientific status:** `OPEN CONJECTURE` (repository's own formalization of
  an article remark, not an article-numbered conjecture).
- **Reproducibility status:** `NOT APPLICABLE` (nothing is claimed as proven
  or certified here; see the proof document for the numerical evidence that
  does exist and its explicitly limited, exploratory status).

## Notation used in this folder

`R`, `phi := pi - Delta_theta`, `w := sqrt(p)`, `t := sqrt(R)`, `z :=
cot^2(phi/2)`, consistent with `docs/mathematical_conventions.md`.
