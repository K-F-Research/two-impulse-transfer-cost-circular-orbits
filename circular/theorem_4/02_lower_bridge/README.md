# Theorem 4 (R=2) — Region 02_lower_bridge

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Scope:** the segment `0.00015 <= s <= 0.15` (`s := sqrt(p)`), bridging
the lower-endpoint region (`01_lower_endpoint`, `0 < s <= 0.00015`) up
to `s = 0.15`. Establishes `H < 0` (strict negativity of the curvature
determinant) throughout this segment by a certified adaptive mesh.

## Contents

- `lower_bridge_proof.md` — the complete argument: the direct
  interval (mean-value) bound on `P` away from `s=0`, the mesh
  construction and coverage check, the tile-to-tile branch-stitching
  argument, the two independent audit levels, and the junction with
  `01_lower_endpoint` at `s=0.00015`.
- `scripts/`:
  - `full_mesh_bridge.py` — builds and certifies the adaptive mesh over
    `[0.00015, 0.15]`; writes the archive consumed by the two audit
    scripts below.
  - `pilot_mesh_audit.py` — independent structural audit (level 1):
    re-checks every criterion from the archived bounds, without any
    new symbolic or interval computation.
  - `pilot_mesh_audit_level2.py` — independent numerical audit (level
    2): reconstructs the symbolic chain separately and re-evaluates
    every tile and stitching by interval arithmetic from the archived
    `W`/`Y0`/`Y_B` boxes alone, cross-checking against the archived
    values.
- `lower_bridge_archive.json` — the complete mesh archive: 464
  tiles, 463 tile-to-tile stitchings, parameters, and per-tile
  certified bounds.
- `lower_bridge_certificate.json` — certification record:
  domain, claims, method distinctions, artifact hashes, dependencies.

## Notation

`s := sqrt(p)` throughout this region's documents. The historical
scripts use `w` for this same quantity (`w = s`), consistent with the
project-wide convention `w := sqrt(p)` in `docs/mathematical_conventions.md`.
This is a **different** quantity from the `w := R^(-1/2)` used in
`theorem_4/05_upper_bridge`'s certificate. See
`lower_bridge_proof.md`, Section 1.

## How to reproduce

```
python3 scripts/full_mesh_bridge.py
python3 scripts/pilot_mesh_audit.py /tmp/full_mesh_archive.json
python3 scripts/pilot_mesh_audit_level2.py /tmp/full_mesh_archive.json
```

`full_mesh_bridge.py` depends on the pickled symbolic form of `P`
(`P_built.pkl`), produced by `01_lower_endpoint/scripts/build_P.py`;
run that script first. `full_mesh_bridge.py` takes a few minutes
(interval certification of 464 tiles with automatic subdivision); the
two audit scripts take well under a minute each. Dependencies: `sympy`,
`mpmath`.
