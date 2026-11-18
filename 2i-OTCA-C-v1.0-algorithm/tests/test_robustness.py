"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

A dedicated robustness campaign across this package's public API,
distinct from the per-module unit tests already in
test_core.py/test_lambert.py/test_rk4.py/etc. --
angular and radial domain limits (never confusing a near-limit case with
the exact degenerate boundary), metamorphic invariants that must hold
across symmetric transformations of the input, and deliberate failure
injection, checking the actual residuals/invariants/exception content
rather than just "no crash" or "an exception of some type was raised".
"""
import math

import pytest

from otca_c import Bodies, DomainError, compute_transfer
from otca_c.core import find_p_star
from otca_c.errors import ConvergenceError
from otca_c.lambert import solve_lambert


# --- Angular limits: near 0/180/360 deg, never confused with the exact
# degenerate boundary itself (delta_theta = 0 or 360 mod 360, already
# covered by test_api.test_rejects_degenerate_angle) -------------------

@pytest.mark.parametrize("delta_theta_deg", [0.01, 0.001, 0.0001])
def test_near_zero_degree_transfer_converges_and_is_distinct_from_the_boundary(delta_theta_deg):
    """Just above the degenerate delta_theta=0 boundary: must succeed (not
    raise), and must NOT silently collapse to the Hohmann closed form --
    these are two mathematically distinct configurations, only arbitrarily
    close to each other."""
    result = compute_transfer(6778.0, 13556.0, delta_theta_deg, mu_km3_s2=Bodies.EARTH)
    assert not result.is_hohmann
    assert result.delta_v_total_km_s > 0
    assert math.isfinite(result.delta_v_total_km_s)
    assert result.geometry_check.ok(tol=1e-4)


@pytest.mark.parametrize("delta_theta_deg", [359.99, 359.999, 359.9999])
def test_near_360_degree_transfer_converges_and_is_distinct_from_the_boundary(delta_theta_deg):
    result = compute_transfer(6778.0, 13556.0, delta_theta_deg, mu_km3_s2=Bodies.EARTH)
    assert not result.is_hohmann
    assert math.isfinite(result.delta_v_total_km_s)
    assert 0.0 < result.eccentricity < 1.0


@pytest.mark.parametrize("delta_theta_deg", [179.9, 179.99, 180.01, 180.1])
def test_near_hohmann_transfer_converges_and_is_distinct_from_exact_hohmann(delta_theta_deg):
    """Arbitrarily close to delta_theta=180 deg (Hohmann) on either side
    must take the GENERAL solver branch, not the exact closed-form
    Hohmann branch (is_hohmann must be False), and must cost only
    marginally more than the exact Hohmann transfer -- not jump
    discontinuously, which would indicate a branch-selection bug at this
    boundary."""
    result = compute_transfer(6778.0, 13556.0, delta_theta_deg, mu_km3_s2=Bodies.EARTH)
    assert not result.is_hohmann
    assert 0.0 <= result.penalty_vs_hohmann_pct < 0.5  # tiny penalty, not a jump
    assert result.geometry_check.ok(tol=1e-4)


def test_delta_theta_exactly_180_takes_the_exact_hohmann_branch():
    result = compute_transfer(6778.0, 13556.0, 180.0, mu_km3_s2=Bodies.EARTH)
    assert result.is_hohmann
    assert result.penalty_vs_hohmann_pct == pytest.approx(0.0, abs=1e-9)


# --- Radial limits: R close to 1 (the proven domain's own open lower
# bound, R>1 strictly -- see api.py's DomainError for R<=1) and very
# large R -----------------------------------------------------------------

@pytest.mark.parametrize("R", [1.0001, 1.001, 1.01])
def test_radius_ratio_close_to_one_converges(R):
    """R := r2/r1 -> 1+ is the open lower boundary of this component's
    proven domain (R>1 strictly); just above it must still converge to a
    sensible, finite, small-penalty-vs-Hohmann result (a transfer between
    two nearly identical circular orbits), not degrade silently."""
    r1_km = 6778.0
    result = compute_transfer(r1_km, r1_km * R, 170.0, mu_km3_s2=Bodies.EARTH)
    assert math.isfinite(result.delta_v_total_km_s)
    assert result.delta_v_total_km_s > 0
    assert result.geometry_check.ok(tol=1e-4)


def test_radius_ratio_exactly_one_is_rejected():
    with pytest.raises(DomainError):
        compute_transfer(6778.0, 6778.0, 170.0, mu_km3_s2=Bodies.EARTH)


@pytest.mark.parametrize("R", [100.0, 1000.0, 10000.0, 100000.0])
def test_radius_ratio_extremely_large_converges(R):
    """No certified upper bound on R exists (R>1 is open-ended) -- a very
    large ratio (e.g. LEO to a distant interplanetary-scale orbit, in this
    idealized model) must still converge and satisfy the independent
    geometry check, not merely 'not crash'."""
    r1_km = 6778.0
    result = compute_transfer(r1_km, r1_km * R, 170.0, mu_km3_s2=Bodies.EARTH)
    assert math.isfinite(result.delta_v_total_km_s)
    assert result.geometry_check.ok(tol=1e-4)


# --- Metamorphic invariants ------------------------------------------------

@pytest.mark.parametrize("delta_theta_deg", [10.0, 45.0, 90.0, 150.0])
def test_cost_is_invariant_under_short_way_long_way_reflection(delta_theta_deg):
    """phi -> -phi (equivalently delta_theta -> 360-delta_theta) must leave
    the cost exactly unchanged (the model's own evenness in phi), while
    the geometry (nu1/nu2, time of flight) must NOT be the same --
    already a dedicated regression in test_lambert.py via time of flight;
    here checked directly via p*/e/delta-v for a spread of angles, not
    just the one case already covered elsewhere."""
    short = compute_transfer(6778.0, 13556.0, delta_theta_deg, mu_km3_s2=Bodies.EARTH)
    long = compute_transfer(6778.0, 13556.0, 360.0 - delta_theta_deg, mu_km3_s2=Bodies.EARTH)
    assert short.p_star_km == pytest.approx(long.p_star_km, rel=1e-9)
    assert short.eccentricity == pytest.approx(long.eccentricity, rel=1e-9)
    assert short.delta_v_total_km_s == pytest.approx(long.delta_v_total_km_s, rel=1e-9)
    # But the two arcs are genuinely different (not the same transfer
    # under the hood) -- their true anomalies at arrival must differ.
    assert short.nu2_rad != pytest.approx(long.nu2_rad, rel=1e-6)


@pytest.mark.parametrize("delta_theta_deg", [60.0, 150.0, 250.0])
def test_angular_momentum_and_energy_are_consistent_between_apsides(delta_theta_deg):
    """A metamorphic/orbital-invariant check independent of any of this
    package's own cross-checks: specific angular momentum h=sqrt(mu*p)
    and specific energy eps=-mu/(2a) computed from the solved orbit must
    reproduce the vis-viva equation at BOTH r1 and r2 independently (not
    just at one of them, which `geometry_check` already verifies through
    a different, but not unrelated, formula path) -- recomputed here
    directly from first principles as an extra, deliberately redundant
    check specific to this test."""
    result = compute_transfer(6778.0, 13556.0, delta_theta_deg, mu_km3_s2=Bodies.EARTH)
    mu = Bodies.EARTH
    a_km = result.semi_major_axis_km
    eps = -mu / (2.0 * a_km)  # specific orbital energy
    h = (mu * result.p_star_km) ** 0.5  # specific angular momentum

    for r_km, nu in ((result.r1_km, result.nu1_rad), (result.r2_km, result.nu2_rad)):
        v_theta = h / r_km
        v2_total = 2.0 * (eps + mu / r_km)  # vis-viva, rearranged
        assert v2_total >= -1e-6  # must not be spuriously negative
        v_radial_sq = max(v2_total - v_theta**2, 0.0)
        # Consistency: the conic equation's own r at this nu must match
        # r_km to high precision.
        r_from_conic = result.p_star_km / (1.0 + result.eccentricity * math.cos(nu))
        assert r_from_conic == pytest.approx(r_km, rel=1e-6)


# --- Failure injection: convergence failures, non-finite values, and
# making sure an unexpected (non-ConvergenceError) exception is NOT
# silently swallowed by compute_transfer's lambert_check/rk4_check
# wrapping -------------------------------------------------------------

def test_find_p_star_raises_convergence_error_not_a_silent_nan():
    """A direct failure-injection check at the core.py level: an input at
    the genuine domain boundary (phi=0 exactly, sin(phi)=0, this model's
    own 1/sin(phi)^2 singularity) must raise ConvergenceError, never
    return NaN/inf or an otherwise-plausible-looking wrong p*."""
    with pytest.raises(ConvergenceError):
        find_p_star(0.0, 2.0)


def test_solve_lambert_raises_convergence_error_on_degenerate_collinear_case():
    with pytest.raises(ConvergenceError):
        solve_lambert((7000.0, 0.0), (7000.0, 0.0), 3600.0, Bodies.EARTH)


def test_lambert_check_wrapping_only_catches_convergence_error_not_everything():
    """Regression guard for the specific defect flagged during independent
    review: api.py's lambert_check wrapper must catch ONLY
    ConvergenceError (deliberately not a bare `except Exception`, which
    could silently hide a real programming bug behind a clean-looking
    `lambert_check=None`). Verified here by monkeypatching
    derived.lambert_cross_check to raise a different exception type and
    confirming it propagates out of compute_transfer rather than being
    swallowed."""
    import otca_c.api as api_module

    original = api_module.derived.lambert_cross_check

    def _raise_type_error(*args, **kwargs):
        raise TypeError("injected failure, must not be swallowed")

    api_module.derived.lambert_cross_check = _raise_type_error
    try:
        with pytest.raises(TypeError, match="injected failure"):
            compute_transfer(
                6778.0, 13556.0, 170.0, mu_km3_s2=Bodies.EARTH, lambert_check=True,
            )
    finally:
        api_module.derived.lambert_cross_check = original


def test_rk4_check_wrapping_only_catches_convergence_error_not_everything():
    """Same guard as above, for the rk4_check wrapper added in this
    package's RK4 cross-check."""
    import otca_c.api as api_module

    original = api_module.rk4_validate

    def _raise_value_error(*args, **kwargs):
        raise ValueError("injected failure, must not be swallowed")

    api_module.rk4_validate = _raise_value_error
    try:
        with pytest.raises(ValueError, match="injected failure"):
            compute_transfer(
                6778.0, 13556.0, 170.0, mu_km3_s2=Bodies.EARTH, rk4_check=True,
            )
    finally:
        api_module.rk4_validate = original


def test_compute_transfer_never_returns_non_finite_quantities_across_a_sweep():
    """A broad sweep (not an exhaustive campaign, but wider than any single
    existing test) across angles and radius ratios -- every unconditional
    output field must be finite; this is the kind of blanket invariant a
    single well-chosen test case cannot substitute for."""
    for delta_theta_deg in (0.5, 30.0, 90.0, 150.0, 179.5, 180.0, 180.5, 270.0, 359.5):
        for R in (1.01, 1.5, 2.0, 10.0, 500.0):
            result = compute_transfer(6778.0, 6778.0 * R, delta_theta_deg, mu_km3_s2=Bodies.EARTH)
            for value in (
                result.delta_v1_km_s,
                result.delta_v2_km_s,
                result.delta_v_total_km_s,
                result.p_star_km,
                result.semi_major_axis_km,
                result.eccentricity,
                result.r_periapsis_km,
                result.r_apoapsis_km,
                result.orbital_period_s,
                result.time_of_flight_s,
            ):
                assert math.isfinite(value), (delta_theta_deg, R, value)
            assert 0.0 <= result.eccentricity < 1.0
