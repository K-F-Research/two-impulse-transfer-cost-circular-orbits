# Theorem 4 (R=2) — Region 03_central_mesh

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Scope:** `w in [0.15, 1.144701]` (`w := sqrt(p)`), continuing
`02_lower_bridge` (`H<0` on `[0.00015, 0.15]`). Establishes that the
curvature determinant `H(p,phi)` changes sign **exactly once** on this
segment, at a Krawczyk-certified point `K_1` (`w_c`), with `H<0` for
`w<w_c` and `H>0` for `w>w_c`.

## Contents

- `central_mesh_proof.md` — the complete argument: the Krawczyk-certified
  sign-change point, the per-tile certification method (interval-Newton
  branch root, `S1,S2>0`, mean-value sign enclosure of `H`, tile-to-tile
  continuity), the results of the two walks, and the junction with
  `02_lower_bridge`.
- `mesh_certificate_R2.py` — the certificate-producing script: builds the
  adaptive mesh in both directions from the Krawczyk-certified `K_1` box
  and certifies the sign of `H` on every tile. Its output is recorded in
  `central_mesh_certificate.json`.
- `mesh_certificate_R2_log.txt` — the complete execution trace: the
  Krawczyk certificate, both walks' tile counts and sign sets, the audit
  counters, and the final summary.
- `central_mesh_certificate.json` — certification record: domain, claims,
  method, artifact hashes, dependencies.

## Notation

This region works directly in the native variables `(p, phi)` of the
cost model; `w := sqrt(p)` is used only to describe the domain and
report results, consistently with the rest of `theorem_4`.

## How to reproduce

```
python3 mesh_certificate_R2.py
```

The script takes several minutes (interval certification of 5142
tiles with adaptive rescaling). It checkpoints each half-walk (DOWN,
then UP) to disk as it completes; if interrupted, rerunning it resumes
from the last completed half-walk rather than recomputing it.
Dependencies: `sympy`, `mpmath`.
