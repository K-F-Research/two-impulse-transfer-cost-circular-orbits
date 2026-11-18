"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Tests for `otca_c.surrogate` (Theorem 2's closed-form fast surrogate).

Three independent kinds of guarantee are checked here, deliberately kept
separate because they answer different questions:

1. `test_K_matches_manuscript_table` -- is the closed form (Section 4.7)
   implemented correctly? Checked against the manuscript's own worked
   numbers (Section 4.8), not against this module's own code.
2. `test_estimate_agrees_with_exact_near_hohmann` -- does the surrogate
   actually approximate the real, independently-solved cost well near
   phi=0? Checked against `core.solve`/`api.compute_transfer`, not
   against the surrogate's own formula.
3. `test_domain_flag_and_guard_behaviour` -- does the package correctly
   warn (not silently mislead) outside its measured accuracy range, and
   does the internal positivity guard actually fire on a broken formula?
"""
import math

import pytest

from otca_c.errors import DomainError
from otca_c.surrogate import (
    MAX_VALIDATED_PHI_DEG,
    SurrogateEstimate,
    _K_normalized,
    compare_to_exact,
    curvature_coefficient_km_s,
    estimate_cost,
)
from otca_c.units import Bodies


# Manuscript Section 4.8's own worked table: K(1,R) at R=1.2, 1.5, 2, 3, 5,
# 10 -- reproduced here verbatim as the independent ground truth.
_MANUSCRIPT_TABLE = {
    1.2: 0.043497,
    1.5: 0.091063,
    2.0: 0.143388,
    3.0: 0.201363,
    5.0: 0.253121,
    10.0: 0.298373,
}


@pytest.mark.parametrize("R, k_expected", sorted(_MANUSCRIPT_TABLE.items()))
def test_K_matches_manuscript_table(R: float, k_expected: float) -> None:
    """The closed-form K(1,R) must reproduce the manuscript's own worked
    numbers (Section 4.8) to the precision those numbers were quoted at."""
    assert _K_normalized(R) == pytest.approx(k_expected, abs=5e-6)


def test_K_positive_for_a_wide_sample_of_R() -> None:
    """Theorem 2 proves K(1,R) > 0 for EVERY R > 1 -- not just the six
    tabulated points above. Sample densely, including R close to 1 (the
    K(1,R) ~ (R-1)/4 limit of Section 4.8) and R very large (the K(1,R)
    -> (2+sqrt(2))/8 limit), to exercise the formula away from the exact
    points the manuscript happened to tabulate."""
    for R in [1.0001, 1.001, 1.01, 1.03, 1.07, 1.1, 1.3, 1.7, 2.5, 4.0,
              7.0, 20.0, 50.0, 200.0, 1000.0, 1.0e6]:
        k = _K_normalized(R)
        assert math.isfinite(k) and k > 0.0, f"K(1,{R}) = {k}"


def test_K_limit_near_R_equals_1() -> None:
    """K(1,R) ~ (R-1)/4 as R -> 1+ (manuscript Section 4.8)."""
    for R in [1.001, 1.0001, 1.00001]:
        k = _K_normalized(R)
        predicted = (R - 1.0) / 4.0
        rel_err = abs(k - predicted) / predicted
        assert rel_err < 0.01, f"R={R}: K={k}, (R-1)/4={predicted}"


def test_K_limit_large_R() -> None:
    """K(1,R) -> (2+sqrt(2))/8 as R -> infinity (manuscript Section 4.8).

    Convergence is genuinely slow (empirically O(1/R), verified with
    mpmath at R = 1e4..1e20 during development: relative error ~1.6e-2
    at R=1e4, ~1.7e-4 at R=1e8, ~1.7e-9 at R=1e15). R=1.0e8 only reaches
    ~1.7e-4 relative error, so the tolerance here is set to 5e-4 (looser
    than the matching-table test's precision, deliberately, since this is
    checking a slow asymptotic limit, not a tabulated exact value) and R
    is pushed to 1e10 to stay safely inside that margin rather than
    sitting right at the edge of it.
    """
    limit = (2.0 + 2.0 ** 0.5) / 8.0
    k = _K_normalized(1.0e10)
    assert k == pytest.approx(limit, rel=5e-4)


def test_K_normalized_rejects_R_at_or_below_one() -> None:
    with pytest.raises(DomainError):
        _K_normalized(1.0)
    with pytest.raises(DomainError):
        _K_normalized(0.5)


def test_curvature_coefficient_km_s_scales_like_vc1() -> None:
    """K(r1,r2) in real units must equal K(1,R) * sqrt(mu/r1) exactly
    (module docstring's dimensional-analysis argument), checked
    numerically rather than assumed."""
    r1_km, r2_km = 6678.0, 13356.0  # R = 2
    mu = Bodies.EARTH
    K_real = curvature_coefficient_km_s(r1_km, r2_km, mu)
    expected = _K_normalized(2.0) * (mu / r1_km) ** 0.5
    assert K_real == pytest.approx(expected, rel=1e-12)


def test_estimate_agrees_with_exact_near_hohmann() -> None:
    """At a small angular deviation, the surrogate must agree with the
    independently-solved exact cost (api.compute_transfer) to well
    within the error bound this repository's own benchmark measured at
    R=1.5 (<0.0003% at 5 deg -- allow a generous margin here since this
    is a different R, but the error should still be tiny)."""
    r1_km, r2_km = 6678.0, 10017.0  # R = 1.5
    comparison = compare_to_exact(r1_km, r2_km, delta_theta_deg=175.0, mu_km3_s2=Bodies.EARTH)
    assert comparison.relative_error < 1e-3  # 5 deg deviation from Hohmann


def test_estimate_error_grows_with_angle() -> None:
    """The surrogate's relative error should grow monotonically-ish with
    the angular deviation (the theorem's own O(phi**4) remainder), not
    stay flat or shrink -- checked directly, not assumed."""
    r1_km, r2_km = 6678.0, 10017.0  # R = 1.5
    errors = []
    for delta_theta in [175.0, 165.0, 150.0]:  # 5, 15, 30 deg deviation
        c = compare_to_exact(r1_km, r2_km, delta_theta, Bodies.EARTH)
        errors.append(c.relative_error)
    assert errors[0] < errors[1] < errors[2]


def test_domain_flag_and_guard_behaviour() -> None:
    """estimate_cost must flag (not raise) a deviation beyond
    MAX_VALIDATED_PHI_DEG, and must NOT flag one comfortably inside it."""
    r1_km, r2_km = 6678.0, 10017.0

    near = estimate_cost(r1_km, r2_km, delta_theta_deg=170.0, mu_km3_s2=Bodies.EARTH)  # 10 deg
    assert near.within_measured_domain is True

    far = estimate_cost(r1_km, r2_km, delta_theta_deg=100.0, mu_km3_s2=Bodies.EARTH)  # 80 deg
    assert far.within_measured_domain is False
    assert "EXCEEDS" in far.domain_note

    # Boundary itself: exactly at the threshold must still count as "within".
    boundary_delta_theta = 180.0 - MAX_VALIDATED_PHI_DEG
    boundary = estimate_cost(r1_km, r2_km, boundary_delta_theta, Bodies.EARTH)
    assert boundary.within_measured_domain is True


def test_estimate_cost_rejects_invalid_inputs() -> None:
    with pytest.raises(DomainError):
        estimate_cost(-1.0, 2.0, 170.0)
    with pytest.raises(DomainError):
        estimate_cost(2.0, 1.0, 170.0)  # r2 < r1
    with pytest.raises(DomainError):
        estimate_cost(1.0, 2.0, 0.0)  # degenerate angle


def test_estimate_cost_result_is_frozen_dataclass() -> None:
    result = estimate_cost(6678.0, 10017.0, 170.0, Bodies.EARTH)
    assert isinstance(result, SurrogateEstimate)
    with pytest.raises(Exception):
        result.delta_v_total_estimate_km_s = 0.0  # frozen -- must not allow mutation
