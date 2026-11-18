"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Adversarial / numerical-stability tests for `otca_c.surrogate`, added to
widen coverage into regions most likely to expose a closed-form formula's
weaknesses: near-degenerate radius ratios, very large radius ratios, the
angular limit phi -> 0 (where the exact solver's own cancellation failure
mode is independently documented in the manuscript, Section 5.6.1), and
malformed/boundary-adjacent numerical input. These are deliberately
distinct from `test_surrogate.py` (fixed reference values),
`test_surrogate_properties.py` (scaling/metamorphic invariants), and
`test_surrogate_extra.py` (domain-flag and package-exposure edge cases):
this file's purpose is specifically to stress the implementation at
inputs an operational user is unlikely to choose deliberately but that a
formally verified component must not silently mishandle.

All assertions check a specific, pre-stated numerical property (finite,
correctly signed, converges in a stated direction) rather than merely
absence of an exception, per this repository's own
never-return-a-silently-wrong-value discipline (see `errors.py`).
"""
import math

import pytest

from otca_c import Bodies, DomainError
from otca_c.surrogate import (
    MAX_VALIDATED_PHI_DEG,
    _K_normalized,
    compare_to_exact,
    curvature_coefficient_km_s,
    estimate_cost,
)

EARTH_MU = Bodies.EARTH


# ---------------------------------------------------------------------------
# Near-degenerate R -> 1+: K(1,R) must stay finite, positive, and track the
# proven (R-1)/4 asymptote increasingly tightly as R approaches 1, down to
# the edge of double-precision resolution. Theorem 2 proves this limit;
# test_surrogate.py already checks it at R=1.001..1.00001, this pushes an
# order of magnitude further to where cancellation first threatens finite
# differences inside the closed form itself.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("R", [1.0 + 1e-6, 1.0 + 1e-8, 1.0 + 1e-10])
def test_K_near_unity_remains_finite_and_tracks_asymptote(R: float) -> None:
    k = _K_normalized(R)
    assert math.isfinite(k) and k > 0.0
    predicted = (R - 1.0) / 4.0
    rel_err = abs(k - predicted) / predicted
    # Looser tolerance than the main asymptote test: this probes closer to
    # R=1, where the closed form's own internal subtractions (a_K, c_K,
    # b_K**2/(4*a_K)) lose relative precision before the physical limit
    # itself is reached -- the point of this test is to confirm that loss
    # stays bounded and does not blow up or change sign, not to demand
    # the same tight tolerance used farther from the degenerate point.
    assert rel_err < 0.05, f"R={R}: K={k}, (R-1)/4={predicted}, rel_err={rel_err}"


def test_K_normalized_rejects_R_extremely_close_to_one_from_below() -> None:
    """R < 1 (even by a tiny margin) is outside Theorem 2's domain
    (r1=1 is defined as the smaller radius); must raise, not silently
    extrapolate."""
    with pytest.raises(DomainError):
        _K_normalized(1.0 - 1e-9)


# ---------------------------------------------------------------------------
# Very large R: K(1,R) must converge toward (2+sqrt(2))/8 from below (per
# the manuscript's own documented O(1/R) convergence) without overshoot,
# oscillation, or loss of finiteness, pushed to R=1e15 -- an order of
# magnitude beyond the existing test_K_limit_large_R (which stops at
# R=1e10) and close to where float64 R**0.5 itself starts losing digits.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("R", [1.0e6, 1.0e10, 1.0e15])
def test_K_large_R_converges_monotonically_without_overshoot(R: float) -> None:
    limit = (2.0 + 2.0**0.5) / 8.0
    k = _K_normalized(R)
    assert math.isfinite(k) and k > 0.0
    assert k < limit, (
        f"K(1,R={R!r})={k} meets or exceeds the proven supremum "
        f"{limit}; the manuscript's O(1/R) convergence is from below, "
        "so this would indicate a sign or algebra defect, not noise."
    )


def test_K_large_R_sequence_is_monotonically_increasing_toward_limit() -> None:
    """A direct, ordered check that K(1,R) increases as R grows (consistent
    with convergence from below to a finite supremum), rather than merely
    checking a few isolated points against the limit separately."""
    Rs = [2.0, 10.0, 1.0e2, 1.0e4, 1.0e6, 1.0e8, 1.0e10]
    ks = [_K_normalized(R) for R in Rs]
    for a, b in zip(ks, ks[1:]):
        assert b > a, f"K(1,R) sequence not monotonically increasing: {ks}"


# ---------------------------------------------------------------------------
# phi -> 0 (delta_theta_deg -> 180 deg from below): the exact solver's own
# finite-difference curvature check is documented (manuscript Section
# 5.6.1, and this repository's k_convergence.py) to suffer catastrophic
# cancellation extremely close to phi=0. The CLOSED-FORM surrogate must
# NOT inherit this failure mode -- it never differences two nearby exact
# solves, so its estimate at a vanishingly small phi should remain stable
# and should converge to the exact Hohmann value, unlike a naive
# finite-difference approach would.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("phi_deg", [1e-3, 1e-6, 1e-9])
def test_estimate_stable_at_vanishingly_small_phi(phi_deg: float) -> None:
    r1_km, r2_km = 6678.0, 10017.0  # R = 1.5
    delta_theta_deg = 180.0 - phi_deg
    result = estimate_cost(r1_km, r2_km, delta_theta_deg, EARTH_MU)
    assert math.isfinite(result.delta_v_total_estimate_km_s)
    assert result.within_measured_domain is True
    # At phi this small the quadratic term is negligible to double
    # precision; the estimate must not drift from the closed-form Hohmann
    # reference it carries internally.
    assert result.delta_v_total_estimate_km_s == pytest.approx(
        result.hohmann_delta_v_total_km_s, rel=1e-6
    )


@pytest.mark.parametrize("phi_deg", [1e-3, 1e-6])
def test_compare_to_exact_stable_at_vanishingly_small_phi(phi_deg: float) -> None:
    """Unlike the surrogate's own internal Hohmann reference (previous
    test), this cross-checks against the independently solved exact cost
    (api.compute_transfer) at the same vanishingly small phi -- the
    regime the manuscript documents as numerically hardest for the exact
    solver's own internals (Section 5.6.1), making this the sharpest
    available test of whether the surrogate stays trustworthy exactly
    where the alternative (a finite-difference curvature estimate) would
    not."""
    r1_km, r2_km = 6678.0, 10017.0
    delta_theta_deg = 180.0 - phi_deg
    comparison = compare_to_exact(r1_km, r2_km, delta_theta_deg, EARTH_MU)
    assert math.isfinite(comparison.relative_error)
    assert comparison.relative_error < 1e-4


# ---------------------------------------------------------------------------
# Extreme radius ratios combined with small phi: the full 2D adversarial
# corner (R near 1 AND phi near 0 simultaneously; R very large AND phi
# near the 30 deg validated boundary simultaneously), not covered by any
# single-axis test above or in the other three surrogate test files.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "r2_over_r1, phi_deg",
    [
        (1.0 + 1e-4, 1e-3),     # near-degenerate R, near-zero phi
        (1.0 + 1e-4, 29.0),     # near-degenerate R, near the validated boundary
        (1.0e4, 1e-3),          # huge R, near-zero phi
        (1.0e4, 29.0),          # huge R, near the validated boundary
    ],
)
def test_estimate_well_formed_at_adversarial_corners(r2_over_r1: float, phi_deg: float) -> None:
    r1_km = 6678.0
    r2_km = r1_km * r2_over_r1
    delta_theta_deg = 180.0 - phi_deg
    result = estimate_cost(r1_km, r2_km, delta_theta_deg, EARTH_MU)
    assert math.isfinite(result.delta_v_total_estimate_km_s)
    assert result.delta_v_total_estimate_km_s > 0.0
    assert math.isfinite(result.K_km_s) and result.K_km_s > 0.0


# ---------------------------------------------------------------------------
# curvature_coefficient_km_s must be well-conditioned under an extreme,
# but legitimate, choice of mu (e.g. a low-mu body such as a small moon)
# combined with a large radius -- neither axis alone is adversarial, the
# combination has not been tested elsewhere.
# ---------------------------------------------------------------------------
def test_curvature_coefficient_well_formed_for_low_mu_large_radius() -> None:
    r1_km, r2_km = 1.0e6, 2.0e6  # large radii (e.g. a distant/weak-field body)
    mu_small = 4.9e3  # order of magnitude of a small moon's GM, km^3/s^2
    K = curvature_coefficient_km_s(r1_km, r2_km, mu_small)
    assert math.isfinite(K) and K > 0.0


# ---------------------------------------------------------------------------
# Input validation at the numerical edge of "positive": a radius that is
# positive but smaller than any physically meaningful orbit (not zero,
# not negative -- those paths are already tested) must still be accepted
# as a pure algebraic positive number, since this module makes no
# assumption about physical units; this documents that behaviour
# explicitly rather than leaving it implicit.
# ---------------------------------------------------------------------------
def test_estimate_accepts_arbitrarily_small_positive_radii() -> None:
    result = estimate_cost(1e-9, 2e-9, 170.0, mu_km3_s2=1.0)
    assert math.isfinite(result.delta_v_total_estimate_km_s)
    assert result.delta_v_total_estimate_km_s >= 0.0


# ---------------------------------------------------------------------------
# curvature_coefficient_km_s's own radii-positivity guard (not exercised via
# estimate_cost or compare_to_exact elsewhere: those already validate radii
# before ever calling this function, so its own independent check line was
# otherwise reachable only through this function's direct, standalone use).
# ---------------------------------------------------------------------------
def test_curvature_coefficient_rejects_negative_or_zero_radii() -> None:
    with pytest.raises(DomainError):
        curvature_coefficient_km_s(-6678.0, 10017.0, EARTH_MU)
    with pytest.raises(DomainError):
        curvature_coefficient_km_s(6678.0, -10017.0, EARTH_MU)
    with pytest.raises(DomainError):
        curvature_coefficient_km_s(0.0, 10017.0, EARTH_MU)
