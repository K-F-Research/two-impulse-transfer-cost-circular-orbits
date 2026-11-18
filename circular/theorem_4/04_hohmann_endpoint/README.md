# Theorem 4 (R=2) — Region 04_hohmann_endpoint

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Scope:** the segment `0 < t <= t_archive` (`t := sin(phi/2)^2`,
`t_archive` the archive-recorded upper boundary, numerically ~1e-4),
the endpoint of the domain covered by `theorem_4` at R=2 nearest the
Hohmann point `p = 4/3`. Establishes `H > 0` (strict positivity of the
curvature determinant) throughout this segment, via the rescaled
quantity `t*H`.

## Contents

- `certificate/hohmann_endpoint_certificate.json` — certification
  record: parameterization, branch existence/uniqueness, radicand
  positivity, the sign of `t*H`, the standalone `C_pp` bound, the
  endpoint-limit consistency check, the junction with
  `05_upper_bridge`, and artifact hashes.
- `hohmann_certificate_v3.py` — the certificate-generating script:
  builds the branch map, certifies existence/uniqueness of the
  stationary root, certifies `S1, S2 > 0`, and certifies `t*H > 0` by
  subdivision of the `a`-interval.
- `hohmann_raccord_t1e4.py` — certifies the junction with
  `05_upper_bridge` at `t = t_archive` (shared with that region; see
  both regions' certificates, junction `J_HOHMANN_UPPER_BRIDGE` /
  `J_HOHMANN_LOCAL_UPPER_BRIDGE`).
- `hohmann_w_bound_exact.py` — independent exact-arithmetic check,
  separate from the `t*H` sign certificate above: certifies, by
  interval arithmetic plus an exact-rational (`fractions.Fraction`)
  equality check on the parsed boundary string, that the lower bound
  `w_-` of the upper-bridge region in the `w = sqrt(p)` variable, taken
  at the junction `t = t_archive`, is strictly below `w_H = sqrt(4/3)`
  (the Hohmann point itself, a limit of the branch as `t -> 0`, never
  attained for `t > 0`). This supports the parameterization note in
  `certificate/hohmann_endpoint_certificate.json`
  (`parameterization.domain.upper_note`) that `t_archive` is used
  verbatim rather than the bare literal `1e-4`; it is not itself one of
  the sign-certificate scripts above and carries its own hash below.
- `hohmann_independent_audit.py` — an independent re-derivation of the
  regularized Hohmann-endpoint quantities (`Ghat(t,a)`, `t*H(t,a)`,
  `S1_hat`, `S2_hat`) built from the primitive cost function along its
  own symbolic path, without importing or reusing any object from
  `hohmann_certificate_v3.py`.
- `numeric_check_tHhat.py` — a high-precision (60-digit) numerical
  cross-check of the regularized quantity `t*H` against a second,
  independent evaluation route (differentiating `S1`, `S2` directly at
  `p = p_H + a*t` rather than reusing the already-assembled
  `S1_hat`/`S2_hat`/`t*C_pphi_hat`/`t*C_pt_hat`/`t*C_tt_hat`
  expressions); agreement is confirmed to a relative error of at most
  `2e-55` across 40 sample points. Output recorded in
  `numeric_check_tHhat_log.txt`.

## Notation

`t := sin(phi/2)^2` is this region's own parameter, distinct from
`s := sqrt(p)` used in `01_lower_endpoint`/`02_lower_bridge` and from
`w := sqrt(p)` used to describe domains in `03_central_mesh`; see
`certificate/hohmann_endpoint_certificate.json`'s `definitions` block
for the exact maps `p(t,a) = 4/3 + a*t` and `phi(t) = 2*arcsin(sqrt(t))`
linking them.

## How to reproduce

```
python3 hohmann_certificate_v3.py
python3 hohmann_raccord_t1e4.py
python3 hohmann_w_bound_exact.py
python3 hohmann_independent_audit.py
python3 numeric_check_tHhat.py
```

`hohmann_raccord_t1e4.py` and `hohmann_w_bound_exact.py` depend on the
frozen symbolic expressions pickle produced by `hohmann_certificate_v3.py`;
run it first. `hohmann_independent_audit.py` and `numeric_check_tHhat.py`
are self-contained and depend on no artifact produced by the other
scripts in this directory. Dependencies: `sympy`, `mpmath`.
