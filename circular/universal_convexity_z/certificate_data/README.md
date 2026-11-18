# Certificate Data for Theorem 3, Section 4

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

## What this folder contains

This folder holds the raw polynomial objects behind the resultant-based
closing argument of `theorem3_universal_convexity_proof.md`, Section 4
("Closing the last lemma"), independently re-verified here -- not merely
re-read.

- `Q_poly.txt`, `auxiliary_H_poly.txt` -- the two source polynomials in
  `(t, w, rho)`, in plain-text term form, from which the resultant chain is
  built. (`auxiliary_H_poly.txt` is an internal auxiliary polynomial used at
  this stage of the derivation; it is unrelated to the curvature quantity
  `H` used in `circular/inflection_uniqueness_R2/`, which is a different
  object in a different folder.)
- `RH_exact_flint.pkl`, `RJ_exact_flint.pkl`, `G_exact_flint.pkl`,
  `A_exact_flint.pkl`, `B_exact_flint.pkl` -- the resultants, their gcd, and
  the two cofactors, each stored as an exact-rational multivariate
  polynomial (`python-flint`'s `fmpq_mpoly`, serialized via `pickle`) in the
  variables `(t, w, rho)`.
- `verify_certificate_consistency.py` -- loads all five objects and
  recomputes, from scratch, the four algebraic identities the closing
  argument of Section 4 depends on.

## What was verified, and what this does and does not establish

**Verified, by direct recomputation (not by re-reading a previously
reported result):**

```
A * G  == RH          (exact polynomial identity)
B * G  == RJ          (exact polynomial identity)
gcd(A, B) == 1         (A and B share no common factor)
gcd(RH, RJ) == G       (G is genuinely the gcd, recomputed independently)
```

All four checks passed exactly (see `verify_certificate_consistency.py`
output, reproduced in the proof document's certification note). The
coprimality of `A` and `B` is precisely the algebraic fact the topological
continuity argument of the proof document's Lemma 4.2 relies on.

**Not re-executed:** the expensive upstream step that produces `RJ` from
`Q` and `auxiliary_H_poly` in the first place (a resultant computation over
polynomials with tens of thousands of terms) has not completed when
attempted from scratch; the five stored objects above were used as given
and checked for mutual consistency, rather than re-derived from the two
source polynomials. This is consistent with how the result was first
produced: the chain was replayed end to end once, from the same two source
polynomials, and reproduced the same term counts and degrees as the stored
objects.

The further step connecting `G` (and the degree-1807 discriminant of
`Big(t,w)` discussed in Section 4 of the proof document) to the specific
named objects `RX, RY, RZ, Big` of the proof document's own notation has
not been re-traced symbol-by-symbol; `RH, RJ, G, A, B` are this stage's own
internal working names from the source computation, not a renotation of
the proof document's variables. What this folder establishes is that the
delivered polynomial data is **not merely asserted but is itself internally
consistent and exactly reproducible**, which is a meaningfully stronger
reproducibility status than "scripts unavailable".

## Reproduction

```
python3 verify_certificate_consistency.py
```

Requires `python-flint`. Runs in under a second (the expensive resultant
computation itself is not re-run by this script; see above).
