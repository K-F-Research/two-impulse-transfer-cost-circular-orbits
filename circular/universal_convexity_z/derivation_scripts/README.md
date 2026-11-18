# Derivation scripts for the Section 4 resultant chain

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

## What this folder contains

Three scripts. The first two are the original derivation scripts that produced the
exact-rational polynomial objects stored in `../certificate_data/`
(`RH_exact_flint.pkl`, `RJ_exact_flint.pkl`, `G_exact_flint.pkl`,
`A_exact_flint.pkl`, `B_exact_flint.pkl`), starting from the two source
polynomials `Q_poly.txt` and `auxiliary_H_poly.txt` also stored there.

- `01_build_discriminant_resultant.py` -- builds `Q` and `H` from their
  stored term lists and computes `RH = Res_rho(Q, H)`, the first link in
  the exact algebraic chain.
- `02_reduction_chain_coprimality.py` -- continues the chain: builds `J`
  from `H` and `Q`, computes `RJ = Res_rho(Q, J)`, factors out the common
  part `G = gcd(RH, RJ)`, and verifies that the resulting cofactors
  `A = RH/G` and `B = RJ/G` are coprime -- the algebraic fact the
  topological continuity argument of `theorem3_universal_convexity_proof.md`
  (Lemma 4.2) relies on.

`../certificate_data/verify_certificate_consistency.py` independently
re-verifies the four identities these scripts establish, directly from the
stored objects, without rerunning the expensive resultant computations.

- `03_lemma_discriminant_closure_certification.py` -- a separate,
  self-contained script (no dependency on `Q_poly.txt`/`auxiliary_H_poly.txt`
  or the stored `.pkl` certificate data: it rebuilds everything from the
  defining quartic `F(t,w,rho)` alone) that performs the certified closure
  of the discriminant-non-vanishing step (Lemma 1.1 / Lemma 5.1 of the
  manuscript): factoring `Disc_rho F` and then the degree-490
  discriminant-in-`w` of its degree-22/12 irreducible factor, certifying
  via `python-flint` ball arithmetic that only two real points `t>1` make
  that discriminant vanish, and confirming by high-precision (`mpmath`)
  root-finding that at both points the corresponding double root in `w`
  falls strictly outside the physical strip `(0, m(t))`. Runs end to end
  in under two minutes.

## Reproduction

Scripts 1 and 2 expect `Q_poly.txt` and `auxiliary_H_poly.txt` (see
`../certificate_data/`) in the working directory, and require
`python-flint`. The resultant computation in script 2 (`RJ`) is the most
expensive step of that chain; the other steps complete quickly. Script 3
is independent of these files and only requires `sympy`, `python-flint`,
and `mpmath`.
