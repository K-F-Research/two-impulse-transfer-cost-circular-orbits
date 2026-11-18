# Symbol Glossary

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

This document is the single source of truth for notation across the paper,
the `circular/` proof folders, and the software package. A symbol is
listed here once, with its canonical definition; a folder-local README
states only that it follows this glossary, rather than repeating the
definition.

## Article-level notation

| Symbol | Definition | Used in |
|---|---|---|
| `R` | ratio of the two circular-orbit radii, `R := r2/r1 > 1` | `circular/` (all folders) |
| `phi` | angular deviation from the Hohmann transfer angle, `phi := pi - Delta_theta` (Hohmann at `phi = 0`) | `circular/` (all folders) |
| `p` | semi-latus rectum of the transfer orbit | `circular/` (all folders) |
| `t` | `t := sqrt(R)` | `circular/hohmann_local_minimum`, `global_monotonicity`, `universal_convexity_z`, `curvature_sign_change_general_R` |
| `w` | `w := sqrt(p)` | `circular/theorem_4`, `universal_convexity_z`, `curvature_sign_change_general_R` |
| `s` | `s := sqrt(p)`, used only as the proof substitution in Lemma 1 | `circular/existence_uniqueness_minimizer` |
| `z` | `z := cot^2(phi/2)`, the reparametrized variable of Theorem 3 | `circular/universal_convexity_z`, `circular/curvature_sign_change_general_R` |
| `x` | `x := cos(phi)` | `circular/theorem_4` only |
| `rho` | ratio of the two impulse-cost radicands' square roots, root of the governing quartic `F(t,w,rho)=0`; the physical branch `rho*(t,w)` satisfies `0<rho*<1` and `N_rho(rho*)>0` | `circular/universal_convexity_z` |
| `V_R(phi)` | optimized two-impulse transfer cost (value function), `V_R(phi) := min_p F(p,phi,R)` | all `circular/` theorem statements |
| `V(z)` | `V_R(phi)` expressed as a function of `z` instead of `phi` | `circular/universal_convexity_z` |
| `K(r1,r2)` | local curvature coefficient at Hohmann, `K(r1,r2) := (1/2) V_R''(0)` | article Theorem 2; **see factor-of-2 note below** |
| `K_circ(R)` | `K_circ(R) := V_R''(0)` (no 1/2 factor) | article Section 4 opening; **see factor-of-2 note below** |

### Factor-of-2 convention between `K(r1,r2)` and `K_circ(R)` (resolved)

`K_circ(R) := V_R''(0)` and `K(r1,r2) := (1/2) V_R''(0)` are two names for
two different normalizations of the same quantity: `K(r1,r2) = (1/2)
K_circ(R)`. The article fixes the Taylor-expansion statement as
`V_R(phi) = V_R(0) + (1/2) K_circ(R) phi^2 + o(phi^2)`, so that `K_circ(R)`
is the literal second derivative and `(1/2) K_circ(R)` is the quadratic
Taylor coefficient. `K(r1,r2)`, used in the abstract and Theorem 2's
statement, is that same Taylor coefficient written under its own name. Any
future text or script computing a curvature value under one name must state
explicitly which of the two it reports, and the conversion factor is exactly
2, not something to re-derive each time.

## Flagged notation collisions

### `H`: two unrelated objects sharing a letter

- `circular/theorem_4/`: `H(p,phi) := C_phiphi*C_pp - C_pphi^2` is the
  curvature (Hessian-determinant) quantity whose sign change along the
  physical stationary branch is the subject of Theorem 4; see
  `circular/theorem_4/README.md` and
  `circular/theorem_4/06_global_audit/global_coverage_summary.md` for the
  single, consistent definition used across all five regions.
- `circular/universal_convexity_z/certificate_data/` and
  `derivation_scripts/`: `H` (stored as `auxiliary_H_poly.txt`) is an
  internal auxiliary polynomial in `(t, w, rho)` used only as a
  resultant-chain intermediate in the Section 4 closing argument. It has no
  relationship to the curvature `H(w)` of `theorem_4`.

These two uses of `H` must never be conflated; each folder's own README
already names the other folder's `H` explicitly to prevent confusion
(`certificate_data/README.md` does this one-directionally -- this glossary
is now the authoritative, bidirectional record).

### `z`: resolved, not a live collision

An earlier draft of `circular/existence_uniqueness_minimizer/lemma1_closure_stationary_point.py`
used `z` for its own proof substitution `sqrt(p)`, which collided with the
unrelated `z := cot^2(phi/2)` of Theorem 3. The script-internal variable was
renamed to `s` to remove the collision; `z` is now reserved, throughout this
repository, for the Theorem 3 meaning only. No further action needed.

### `w` vs `s`: same meaning, different folders, not a collision

`w := sqrt(p)` (`theorem_4`, `universal_convexity_z`,
`curvature_sign_change_general_R`) and `s := sqrt(p)`
(`existence_uniqueness_minimizer`) denote the same quantity under two
different names, used in folders that do not reference each other's
internal variables. Not an error, but worth knowing: a reader implementing
both proofs side by side should expect `w` and `s` to be interchangeable.

## Flagged cross-document numbering inconsistency

The published article numbers the
four central theorems as: Theorem 1 = global monotonicity, Theorem 2 =
local curvature at Hohmann, Theorem 3 = universal convexity in `z`, Theorem
4 = certified curvature sign change at R=2. The proof-folder filenames use a
**different, swapped** internal numbering for the first two results:

- `circular/hohmann_local_minimum/theorem1_hohmann_strict_local_minimum.md`
  is filed as "Theorem 1" but its content is the article's **Theorem 2**
  (local curvature / strict local minimum at Hohmann).
- `circular/global_monotonicity/theorem2_global_monotonicity_proof.md` is
  filed as "Theorem 2" but its content is the article's **Theorem 1**
  (global monotonicity).

Theorem 3 and Theorem 4 are numbered consistently between the article and
their respective folders (`universal_convexity_z`, `theorem_4`).
This is a real inconsistency, not a typo in this glossary: the two proof
documents' own titles and internal theorem labels use the swapped
numbering throughout their text, so rather than renaming the two
folders/files, an explicit "Numbering note" cross-referencing this table
has been added to the top of both proof documents
(`theorem1_hohmann_strict_local_minimum.md`,
`theorem2_global_monotonicity_proof.md`) and both folder READMEs
(`circular/hohmann_local_minimum/README.md`,
`circular/global_monotonicity/README.md`). Any citation from the article
to "the Theorem 1 proof folder" or "the Theorem 2 proof folder" should
still be checked against this table, not assumed.

## Units

Radii, velocities, and angles: each folder's own README states which
convention its own scripts use (normalized radii throughout the circular
proof folders; both degrees and radians appear in
different numerical scripts and are stated locally at each use).
