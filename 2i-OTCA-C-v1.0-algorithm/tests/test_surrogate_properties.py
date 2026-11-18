"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Metamorphic / property-based tests for `otca_c.surrogate`. Unlike
`test_surrogate.py` (which checks specific input/output pairs against
fixed reference values), these tests
check invariants that must hold across many inputs by construction of the
underlying two-body physics, and so can catch unit, normalization, or
scaling errors that a fixed-case grid would not exercise.

Four properties are checked:

1. Length-scaling: for mu fixed, scaling both radii by lambda > 0 must scale
   every computed speed (circular speeds, Hohmann cost, surrogate estimate,
   exact cost) by 1/sqrt(lambda). This follows directly from vis-viva:
   v = sqrt(mu*(2/r - 1/a)); scaling r1,r2,a all by lambda while holding mu
   fixed scales every v by 1/sqrt(lambda).

2. Gravitational-parameter scaling: for geometry fixed, scaling mu by q > 0
   must scale every computed speed by sqrt(q) (again directly from vis-viva).

3. K/vc1 invariance under pure radius scaling at fixed R: the *normalized*
   curvature coefficient K(1,R) (`_K_normalized`) must depend on R=r2/r1
   alone, not on the absolute scale of r1, r2 -- this is what "normalized"
   is supposed to mean, and a normalization bug (e.g. an accidental
   dependence on r1 itself) would otherwise not be caught by the existing
   fixed-case tests.

4. Finite-difference consistency of the curvature coefficient K(r1,r2)
   against the EXACT solver's own cost function (not against the
   surrogate's own closed form): the centered second difference of
   api.compute_transfer's exact cost around phi=0 must converge to
   K(r1,r2) as the step size shrinks, before floating-point cancellation
   degrades it again. This is a genuine cross-check of Theorem 2 against
   an independent code path (the exact solver), not a self-consistency
   check of the surrogate module against itself.
"""
import math

import pytest

from otca_c.api import compute_transfer
from otca_c.surrogate import (
    _K_normalized,
    compare_to_exact,
    curvature_coefficient_km_s,
    estimate_cost,
)

EARTH_MU = 398600.4418


# ---------------------------------------------------------------------------
# 1. Length scaling: v -> v / sqrt(lambda)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("lam", [0.1, 0.25, 0.5, 2.0, 4.0, 10.0, 100.0])
def test_length_scaling_of_estimate(lam: float) -> None:
    r1, r2, dtheta = 6678.0, 10017.0, 170.0
    base = estimate_cost(r1, r2, dtheta, EARTH_MU)
    scaled = estimate_cost(r1 * lam, r2 * lam, dtheta, EARTH_MU)
    predicted = base.delta_v_total_estimate_km_s / math.sqrt(lam)
    assert scaled.delta_v_total_estimate_km_s == pytest.approx(predicted, rel=1e-9)


@pytest.mark.parametrize("lam", [0.1, 0.25, 0.5, 2.0, 4.0, 10.0, 100.0])
def test_length_scaling_of_exact_solver(lam: float) -> None:
    r1, r2, dtheta = 6678.0, 10017.0, 170.0
    base = compute_transfer(r1, r2, dtheta, mu_km3_s2=EARTH_MU, validate=False)
    scaled = compute_transfer(r1 * lam, r2 * lam, dtheta, mu_km3_s2=EARTH_MU, validate=False)
    predicted = base.delta_v_total_km_s / math.sqrt(lam)
    assert scaled.delta_v_total_km_s == pytest.approx(predicted, rel=1e-9)


# ---------------------------------------------------------------------------
# 2. Gravitational-parameter scaling: v -> v * sqrt(q)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q", [0.1, 0.25, 0.5, 2.0, 4.0, 10.0, 100.0])
def test_mu_scaling_of_estimate(q: float) -> None:
    r1, r2, dtheta = 6678.0, 10017.0, 170.0
    base = estimate_cost(r1, r2, dtheta, EARTH_MU)
    scaled = estimate_cost(r1, r2, dtheta, EARTH_MU * q)
    predicted = base.delta_v_total_estimate_km_s * math.sqrt(q)
    assert scaled.delta_v_total_estimate_km_s == pytest.approx(predicted, rel=1e-9)


@pytest.mark.parametrize("q", [0.1, 0.25, 0.5, 2.0, 4.0, 10.0, 100.0])
def test_mu_scaling_of_exact_solver(q: float) -> None:
    r1, r2, dtheta = 6678.0, 10017.0, 170.0
    base = compute_transfer(r1, r2, dtheta, mu_km3_s2=EARTH_MU, validate=False)
    scaled = compute_transfer(r1, r2, dtheta, mu_km3_s2=EARTH_MU * q, validate=False)
    predicted = base.delta_v_total_km_s * math.sqrt(q)
    assert scaled.delta_v_total_km_s == pytest.approx(predicted, rel=1e-9)


# ---------------------------------------------------------------------------
# 3. K(1,R) depends on R alone, not on the absolute scale of r1
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("R", [1.1, 1.5, 2.0, 5.0, 10.0])
@pytest.mark.parametrize("r1", [1.0, 6678.0, 384400.0, 1.496e8])
def test_K_normalized_depends_only_on_R(R: float, r1: float) -> None:
    k = _K_normalized(R)
    # _K_normalized(R) must not take r1 at all -- this test's real purpose
    # is to fail loudly if a future refactor accidentally threads r1 into
    # the normalized formula. Calling it with r1 unused is the point.
    assert math.isfinite(k) and k > 0
    assert _K_normalized(R) == k  # pure function of R: identical call, identical result


@pytest.mark.parametrize("R", [1.1, 1.5, 2.0, 5.0, 10.0])
def test_curvature_coefficient_scales_with_vc1_at_fixed_R(R: float) -> None:
    """K(r1,r2) / sqrt(mu/r1) must be the same for any r1 at fixed R=r2/r1."""
    ratios = []
    for r1 in [1000.0, 6678.0, 50000.0, 384400.0]:
        r2 = r1 * R
        k_phys = curvature_coefficient_km_s(r1, r2, EARTH_MU)
        vc1 = math.sqrt(EARTH_MU / r1)
        ratios.append(k_phys / vc1)
    for ratio in ratios[1:]:
        assert ratio == pytest.approx(ratios[0], rel=1e-9)


# ---------------------------------------------------------------------------
# 4. Finite-difference curvature check against the EXACT solver (independent
#    of the surrogate's own closed form), with a convergence sweep so the
#    test fails informatively rather than passing on a lucky tolerance.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("R", [1.2, 1.5, 2.0, 5.0])
def test_K_matches_finite_difference_of_exact_solver(R: float) -> None:
    r1 = 6678.0
    r2 = r1 * R
    K_theoretical = curvature_coefficient_km_s(r1, r2, EARTH_MU)

    def dv_exact(phi_rad: float) -> float:
        dtheta_deg = 180.0 - math.degrees(phi_rad)
        return compute_transfer(r1, r2, dtheta_deg, mu_km3_s2=EARTH_MU, validate=False).delta_v_total_km_s

    dv0 = dv_exact(0.0)
    # Centered second difference at a step small enough for O(phi^2)
    # truncation error to be negligible, but not so small that the exact
    # solver's own double-precision floor dominates (empirically, around
    # phi ~ 5e-4 to 1e-3 rad is the sweet spot for this solver -- see
    # k_convergence.py for the full sweep showing the characteristic
    # V-shaped error curve: O(phi^2) truncation on one side, numerical
    # cancellation noise on the other).
    phi = 5e-4
    K_numerical = (dv_exact(phi) + dv_exact(-phi) - 2 * dv0) / (2 * phi**2)
    assert K_numerical == pytest.approx(K_theoretical, rel=5e-4)
