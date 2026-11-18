# Global Existence and Uniqueness of the Minimizer p*(phi)

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Scientific status:** SUPPLEMENTARY (not part of the published manuscript,
Reference [8]) -- scientifically PROVEN -- ANALYTICALLY on its own terms.
**Reproducibility status:** REPRODUCED (technical audit, 2026-10-25).

**Not an article lemma.** This folder's result does **not** correspond to a
numbered lemma in the published manuscript. The article's own Theorem 1
(global monotonicity) is proved by a stronger, more direct argument that
explicitly does **not** require uniqueness of the minimizer -- see the
article's abstract ("Theorem 1 proves that V_R is strictly increasing on
(0,pi) for every R > 1, without requiring uniqueness of the minimizer") and
Section 3 ("We do not assume uniqueness of the minimizer. What is
established below is in fact stronger than uniqueness would require..."). An
earlier line of work on this project pursued uniqueness as "Assumption H"
before that stronger, uniqueness-free argument was found for Theorem 1; this
folder is the completed proof of that earlier, no-longer-necessary
assumption, kept because it is true, proved, and still useful (see below),
not because the published Theorem 1 depends on it.

**Why it is kept in this repository.** The package's numerical solver
(`otca_c.core`) finds p* by root-finding: it looks for a point where the
cost function's derivative with respect to p vanishes. That procedure is
only guaranteed to return *the* global minimum -- not an arbitrary local one
-- if the global minimizer is unique. This folder's result is what justifies
that guarantee for the software, independently of how the published
Theorem 1 itself is proved.

**Also a dependency of `theorem_4/05_upper_bridge`.** The region's
certificate for `theorem_4` at R=2 relies on global uniqueness of the
stationary root at fixed `phi`, to establish that the certified Hohmann
tile family represents the entire physical branch segment it covers (see
`theorem_4/05_upper_bridge/upper_bridge_certificate.json`,
field `branch_identification_note.closure`). The present folder is that
dependency's source.

## Statement

For every R>1 and every phi in (0,pi), the cost function C(p,phi,R) admits a
unique global minimizer p*(phi,R) over p>0.

## Contents

- `lemma1_closure_stationary_point.py` -- closes the last analytic step of the
  proof (boundary behaviour as s -> 0+ and s -> +infinity, s := sqrt(p)),
  combined with the non-degeneracy result (C_ss > 0 at every stationary
  point, historically called "Assumption H," proved unconditionally for
  every R>1 and every phi in (0,pi); see `assumption_H_proof.md`). Together
  these give existence and uniqueness without reservation.
- `lemma1_stress_test_2000_trials.py` -- independent numerical corroboration
  over 2000 random parameter pairs; not a proof, a cross-check.
- `assumption_H_proof.md` -- the complete elementary proof of Assumption H
  (`C_pp > 0` at every stationary point, for every R>1 and phi in (0,pi)),
  the non-degeneracy result on which `lemma1_closure_stationary_point.py`
  depends. Not shipped elsewhere in this repository prior to this folder.
- `assumption_H_proof_certificate.py` -- the algebraic certificate for
  Assumption H, working directly in the substituted variables (w,z,u).
- `assumption_H_independent_verification.py` -- a from-scratch
  reconstruction of every step of that certificate, starting from the
  manuscript's original S1, S2 formulas rather than from the
  already-substituted variables.
- `assumption_H_certificate.json` -- certification record: dependency
  chain, scripts, hashes, and verification results.

## Notation

See `docs/mathematical_conventions.md` for the full symbol table. This folder
uses `s := sqrt(p)`, chosen specifically to avoid a collision with the
unrelated `z := cot^2(phi/2)` substitution used in
`circular/universal_convexity_z/`, Theorem 3. `assumption_H_proof.md` and its
two scripts use `z := sqrt(p)` as well (consistent with `s` above; the
collaborator's original script used the symbol `z` for this quantity, kept
here as `z` rather than renamed to `s`, since both names denote the same
object and the certificate records the correspondence explicitly). All
scripts in this folder have been verified to produce consistent results
under this convention.

## How to reproduce

```
python3 lemma1_closure_stationary_point.py
python3 lemma1_stress_test_2000_trials.py
python3 assumption_H_proof_certificate.py
python3 assumption_H_independent_verification.py
```

No external inputs; all four scripts are self-contained. Dependencies:
`sympy`, and `mpmath` for the stress test.
