# Theorem 4 (R=2) — Region 06_global_audit

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Scope:** cross-cutting checks that span the whole of `theorem_4`
rather than any single region: consistency of sign conventions and
variable-transformation factors across all six regions, a symbolic
re-verification of one derivation step flagged as numerically heavy,
and an independent re-certification of one region's archived interval
data from its own frozen inputs.

## Contents

- `scripts/audit_sign_conventions_R2.py` — checks that every factor
  used to relate a signed quantity to `H` itself (`w^6` in `P := w^6*H`;
  `t` in `t*H`; `t(1-t)` in the `C_phiphi` identity; the Jacobian
  `dt/dphi`; the map `p = w^2`; the Krawczyk preconditioner `Y`) has a
  constant, known sign on the relevant physical domain, so that no sign
  inversion can be introduced silently by a rescaling or a coordinate
  change between regions. 25 checks, all passed (see
  `archive/audit_sign_conventions_R2_log.txt`).
- `scripts/audit_transformations_R2.py` — independently re-derives, in
  purely rational form (the two square roots `sqrt(S1)`, `sqrt(S2)` are
  kept as opaque symbols rather than expanded, to keep the symbolic
  comparisons tractable), the variable-transformation identities used
  across the endpoint regularizations: `x = 1-2t`, the `x`-to-`t` map
  of the Hessian, the Hohmann regularization `Ghat`, the assembly of
  `t*H`, and the closed form of the Hohmann-endpoint limit `L(a)`. 14
  checks, all passed (see `archive/audit_transformations_R2_log.txt`).
- `scripts/audit_section4_symbolic_light.py` — re-verifies, by a
  memory-light method (expansion of the isolated numerator only, never
  a full `cancel()`/`factor()` of the combined expression), that one
  numerator required to vanish identically in the Hohmann-endpoint
  derivation does in fact vanish identically. Result: confirmed (see
  `archive/audit_section4_symbolic_light_log.txt`).
- `scripts/audit_hohmann_archive.py` — re-certifies the archived
  interval data of `05_upper_bridge` (`hohmann_full_mesh_archive.json`)
  directly from its own frozen, hash-checked inputs: it does not reuse
  any status field recorded by the generating script, and it
  recomputes the four certification criteria (radicand positivity,
  strict positivity of the two rescaled quantities, tile-to-tile
  junction consistency) and the domain-coverage topology from the
  archived interval bounds alone, after applying an explicit outward
  safety margin to each bound to compensate for the non-directed
  rounding of its decimal serialization. It runs the full
  re-certification at two independent working precisions (50 and 100
  decimal digits) and confirms that both precisions agree on every
  conclusion. Result: both precisions fully re-certify the archive, 0
  gaps, 0 junction failures (see
  `archive/audit_hohmann_archive_log.txt`).

## What this region does not establish

These four checks are consistency and re-certification audits of
specific, already-certified pieces of the proof; they are not a
substitute for, and do not themselves perform, an independent review
of the proof by a party other than its authors. See
`global_coverage_summary.md` in this same directory for the explicit
statement of what remains to be done before the proof can be
considered independently audited in that sense.

## How to reproduce

```
python3 scripts/audit_sign_conventions_R2.py
python3 scripts/audit_transformations_R2.py
python3 scripts/audit_section4_symbolic_light.py
python3 scripts/audit_hohmann_archive.py
```

`audit_hohmann_archive.py` depends on the archived interval data and
frozen symbolic expressions produced in `05_upper_bridge/`; the paths
it reads are set at the top of the script and should be adjusted to
the local checkout if they do not already point to that region's
archived files. Dependencies: `sympy`, `mpmath`.

## Copyright

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
