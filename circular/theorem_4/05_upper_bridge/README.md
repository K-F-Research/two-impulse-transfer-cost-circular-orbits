# Theorem 4 (R=2) — Region 05_upper_bridge

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Scope:** `t_archive <= t <= 0.19055605997847746...` (`t := sin(phi/2)^2`),
continuing `04_hohmann_endpoint` up to the junction with `03_central_mesh`
at `w = 1.144701` (`w := sqrt(p)`). Establishes `H > 0` (strict positivity
of the curvature determinant) throughout this segment, via the rescaled
quantity `t*H`, and certifies the branch-identity junctions with both
neighboring regions.

## Contents

- `upper_bridge_certificate.json` — certification record: the 98-tile
  adaptive mesh in `t`, strict positivity of `S1`, `S2`, `t*H`, and
  `C_pp`, and the two junction certificates (`J_HOHMANN_LOCAL_UPPER_BRIDGE`
  with `04_hohmann_endpoint`, `J_UPPER_BRIDGE_CENTRAL_MESH` with
  `03_central_mesh`), with artifact hashes and provenance.
- `hohmann_full_mesh.py` — the certificate-generating script: builds the
  98-tile adaptive mesh over the region's `t`-domain and certifies the
  sign of `t*H` on every tile. Its archived output is
  `hohmann_full_mesh_archive.json`.
- `hohmann_full_mesh_archive.json` — the complete mesh archive consumed
  by `upper_bridge_certificate.json` and by `06_global_audit`'s
  independent re-certification.
- `hohmann_central_join.py` — certifies the junction with
  `03_central_mesh`: non-empty overlap of the two regions' `p`- and
  `phi`-boxes at the shared boundary, strict sign of `C_pphi` on the
  overlap (uniqueness of `phi` at fixed `p`), and existence of a
  stationary root on both sides throughout the overlap. Output recorded
  in `hohmann_central_join_log_en.txt`.
- `hohmann_central_join_log_en.txt` — the complete execution trace of
  the junction check above.

The junction with `04_hohmann_endpoint` at `t = t_archive` is certified by
`hohmann_raccord_t1e4.py`, filed in `04_hohmann_endpoint/` since it also
re-verifies quantities local to that region; see
`04_hohmann_endpoint/README.md`.

## Notation

`t := sin(phi/2)^2`, consistent with `04_hohmann_endpoint`; `w := sqrt(p)`
is used only to describe the junction with `03_central_mesh`. See
`upper_bridge_certificate.json`'s `definitions` block for the exact maps
linking `t` to `(p, phi)`.

## How to reproduce

```
python3 hohmann_full_mesh.py
python3 hohmann_central_join.py
```

`hohmann_full_mesh.py` takes several minutes (interval certification of
98 tiles). `hohmann_central_join.py` depends on the archived mesh produced
by `hohmann_full_mesh.py` and on `03_central_mesh`'s own certified
boundary tile; run `hohmann_full_mesh.py` and
`03_central_mesh/mesh_certificate_R2.py` first. Dependencies: `sympy`,
`mpmath`.
