# Theorem 1 -- Independent Numerical Checks

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Ten independent numerical checks supporting the global monotonicity proof
(`../theorem2_global_monotonicity_proof.md`; the main independent
verification script, `../verify_monotonicity_independent.py`, is at the
parent level). Each script below uses a distinct method, precision level,
or input range, so together they corroborate the theorem from several
independent angles rather than repeating the same computation.

## Contents

- `bracket_strategy.py` -- confirms, at arbitrary precision, that the sign
  of dF/dp at the two closed-form points p1*, p2* alone localizes the true
  minimizer p* without solving for it explicitly (200 configurations,
  consistent with Lemma 7.3 of the proof document).
- `independent_fg_check.py` -- re-derives the two-impulse Delta-v total via
  the classical Lagrange f,g coefficients (vector Lambert formulation)
  instead of the vis-viva decomposition used throughout the proof, ruling
  out a shared transcription error between the two formulations.
- `newton_root_localization_check.py` -- independent high-precision
  localization check using Newton root-finding on dF/dp=0.
- `bisection_root_localization_check.py` -- independent high-precision
  localization check using bisection root-finding instead of Newton's
  method, over angles from 1 deg to 179 deg.
- `bisection_root_localization_check_near_pi.py` -- extends the bisection
  check to angles from 150 deg to 179.9 deg, where the minimizer shrinks
  quickly and requires a refined search bracket. Combined with
  `bisection_root_localization_check.py`, this gives the 243-configuration
  arbitrary-precision localization check cited in the proof document.
- `random_stress_test.py` -- large-scale random sampling (50,000 draws)
  of the localization and monotonicity properties over a wide range of R
  and phi, well beyond the hand-picked grids used elsewhere.
- `high_precision_recheck_flagged_cases.py` -- re-checks, at arbitrary
  precision, the handful of apparent anomalies flagged by
  `random_stress_test.py` in the extreme-R, near-boundary regime, showing
  they are double-precision resolution artifacts rather than violations of
  the theorem.
- `stress_test_double_precision.py` -- targeted search for pathological
  regimes (R -> 1+, R -> infinity, multiple local minima in p, phi near 0
  or pi) at double precision.
- `stress_test_high_precision.py` -- arbitrary-precision re-check of the
  boundary cases flagged by `stress_test_double_precision.py`.
- `R_less_than_1_and_magnitude_test.py` -- two further checks: that the
  localization mechanism is symmetric under R<1 (relabeling which body is
  "1" and which is "2"), and that V_R(phi2)-V_R(phi1) scales linearly with
  small angular gaps.

## Reproduction

Each script is self-contained:

```
python3 bracket_strategy.py
python3 independent_fg_check.py
python3 newton_root_localization_check.py
python3 bisection_root_localization_check.py
python3 bisection_root_localization_check_near_pi.py
python3 random_stress_test.py
python3 high_precision_recheck_flagged_cases.py
python3 stress_test_double_precision.py
python3 stress_test_high_precision.py
python3 R_less_than_1_and_magnitude_test.py
```

Requires `numpy`, `scipy`, and `mpmath`.
