"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Tests for the independent RK4 numerical-propagation cross-check
(otca_c.rk4, wired into otca_c.api.compute_transfer as `rk4_check`), a
third cross-check genuinely distinct from both the closed-form geometry
check and the Lambert boundary-value cross-check: an initial-value
numerical integrator that exercises time_of_flight_s against the actual
two-body dynamics, across a grid of transfer angles (including a
near-periapsis-passage high-eccentricity case and the near-360-degree
boundary) and several central bodies.
"""
import math

import pytest

from otca_c import Bodies, compute_transfer
from otca_c.errors import ConvergenceError
from otca_c.rk4 import _state_vector_from_elements, rk4_validate


CASES = [
    (6778.0, 13556.0, 170.0, Bodies.EARTH),
    (6778.0, 13556.0, 10.0, Bodies.EARTH),
    (6778.0, 13556.0, 190.0, Bodies.EARTH),  # long way, > 180 deg
    (6778.0, 13556.0, 180.0, Bodies.EARTH),  # exact Hohmann
    (6778.0, 42164.0, 160.0, Bodies.EARTH),
    (3789.9, 24190.9, 170.0, Bodies.MARS),
    (1000.0, 7000.0, 90.0, Bodies.MOON),
]


@pytest.mark.parametrize("r1_km,r2_km,delta_theta_deg,mu", CASES)
def test_rk4_check_matches_to_high_precision(r1_km, r2_km, delta_theta_deg, mu):
    result = compute_transfer(
        r1_km, r2_km, delta_theta_deg, mu_km3_s2=mu, rk4_check=True,
    )
    assert result.rk4_check is not None
    assert result.rk4_check.ok()
    # The integrator must not silently drift in energy/angular momentum
    # either -- this is a sanity check on the integrator itself,
    # independent of the arrival-state comparison above.
    assert result.rk4_check.max_energy_drift < 1e-6
    assert result.rk4_check.max_ang_mom_drift < 1e-6


def test_rk4_check_is_none_by_default():
    result = compute_transfer(6778.0, 13556.0, 170.0, mu_km3_s2=Bodies.EARTH)
    assert result.rk4_check is None


def test_rk4_check_adapts_step_count_for_a_near_periapsis_passage_case():
    """delta_theta=350 deg, R=2.0 around Earth has e=0.982 and its swept
    arc (170-350 deg relative to periapsis) passes directly through
    periapsis -- the sharp curvature there needs far more RK4 steps than
    this package's well-conditioned default (confirmed empirically: 8000
    steps alone gives a ~52 km position residual, badly outside
    tolerance; adaptive doubling must keep refining until it converges)."""
    result = compute_transfer(
        6778.0, 13556.0, 350.0, mu_km3_s2=Bodies.EARTH, rk4_check=True,
    )
    assert result.eccentricity > 0.9
    assert result.rk4_check is not None
    assert result.rk4_check.ok()
    assert result.rk4_check.n_steps > 8000


def test_rk4_check_is_none_near_the_359_999_degree_boundary():
    """This extreme near-360-degree, near-e=1 case (see
    test_core.test_find_p_star_near_360_degree_boundary_regression and
    test_derived.test_compute_transfer_warns_near_the_359_999_degree_boundary)
    is not expected to converge within this function's default
    max_doublings/n_steps budget -- rk4_check must degrade to None
    (an expected, reported non-availability), exactly like
    lambert_check=None for the exact Hohmann case, never a silently
    wrong or a hard-crashing result for an otherwise successful solve."""
    result = compute_transfer(
        6778.0, 13556.0, 359.999, mu_km3_s2=Bodies.EARTH, rk4_check=True,
    )
    assert result.rk4_check is None
    # The primary solve itself must still have succeeded.
    assert 0 < result.eccentricity < 1


def test_state_vector_round_trips_through_the_conic_equation():
    """A basic self-consistency check on the state-vector formulas
    themselves (independent of the integrator): at theta=0 the radius
    must be exactly r1 = p/(1+e*cos(nu1)), and the speed must match the
    vis-viva equation for that radius and the semi-major axis."""
    p_km, e, nu1_rad = 9033.391, 0.334498, -1.2
    mu = Bodies.EARTH
    a_km = p_km / (1.0 - e**2)
    pos, vel = _state_vector_from_elements(p_km, e, nu1_rad, 0.0, mu)
    r = (pos[0] ** 2 + pos[1] ** 2) ** 0.5
    expected_r = p_km / (1.0 + e * math.cos(nu1_rad))
    assert r == pytest.approx(expected_r, rel=1e-12)
    speed_sq = vel[0] ** 2 + vel[1] ** 2
    expected_speed_sq = mu * (2.0 / r - 1.0 / a_km)  # vis-viva
    assert speed_sq == pytest.approx(expected_speed_sq, rel=1e-9)


def test_rk4_validate_raises_on_non_elliptical_eccentricity():
    with pytest.raises(ConvergenceError):
        rk4_validate(
            p_star_km=9000.0, e=1.0, nu1_rad=0.0, delta_theta_rad=1.0,
            tof_s=3600.0, mu_km3_s2=Bodies.EARTH, r2_km=13556.0,
        )


def test_rk4_check_near_periapsis_passage_regression_pins_residual_and_order():
    """Regression test for the real case found during this module's own
    development: delta_theta=350 deg, R=2.0 around Earth (e=0.982) needed
    128000 steps (5 doublings from the n_steps=4000 default) to converge,
    where 8000 steps alone gave a ~52 km position residual -- this pins
    both the final converged residual AND the observed convergence order,
    not just that the call eventually returns a non-None result (which
    would pass even if a future change silently widened the tolerance or
    weakened the order check)."""
    result = compute_transfer(
        6778.0, 13556.0, 350.0, mu_km3_s2=Bodies.EARTH, rk4_check=True,
    )
    c = result.rk4_check
    assert c is not None
    assert c.n_steps == 128000
    assert c.position_residual_km < 1e-3
    assert c.velocity_residual_km_s < 1e-6
    # The order check fired at least once along the way (residuals were
    # above the noise floor before the final, fully-converged doubling);
    # when it did fire it must have been consistent with genuine 4th
    # order, not just "close enough to not raise".
    if c.observed_convergence_order is not None:
        assert c.observed_convergence_order > 2.5


def test_rk4_validate_raises_cleanly_when_the_step_count_cap_is_hit():
    """With max_doublings deliberately set too low for the 350 deg
    near-periapsis-passage case (which needs 5 doublings from n_steps=4000
    to converge -- see the regression test above), rk4_validate must raise
    ConvergenceError with a message that names the actual reason (tolerance
    not reached within the doubling budget) distinctly from the
    4th-order-violation message raised elsewhere -- a caller or a future
    maintainer reading the exception text must be able to tell "needs more
    steps" apart from "something is actually broken"."""
    result = compute_transfer(6778.0, 13556.0, 350.0, mu_km3_s2=Bodies.EARTH)
    with pytest.raises(ConvergenceError, match="did not reach tolerance"):
        rk4_validate(
            p_star_km=result.p_star_km,
            e=result.eccentricity,
            nu1_rad=result.nu1_rad,
            delta_theta_rad=result.nu2_rad - result.nu1_rad,
            tof_s=result.time_of_flight_s,
            mu_km3_s2=Bodies.EARTH,
            r2_km=result.r2_km,
            n_steps=4000,
            max_doublings=2,
        )


def test_rk4_validate_raises_distinctly_on_4th_order_violation_near_359_999():
    """The near-360-degree, near-e=1 case (1-e of order 1e-10) fails for a
    DIFFERENT reason than the step-count-cap case above: the departure
    state itself (p* ~ 2e-6 km, an almost perfectly radial/degenerate
    ellipse) is different enough that successive doublings do not even
    show the expected convergence order -- this must be reported with the
    "did not show the expected ~4th-order convergence" message, not
    silently folded into the generic "did not reach tolerance" message,
    so the two distinct failure modes stay distinguishable."""
    result = compute_transfer(6778.0, 13556.0, 359.999, mu_km3_s2=Bodies.EARTH)
    with pytest.raises(ConvergenceError, match="did not show the expected"):
        rk4_validate(
            p_star_km=result.p_star_km,
            e=result.eccentricity,
            nu1_rad=result.nu1_rad,
            delta_theta_rad=result.nu2_rad - result.nu1_rad,
            tof_s=result.time_of_flight_s,
            mu_km3_s2=Bodies.EARTH,
            r2_km=result.r2_km,
        )


def test_rk4_validate_raises_on_nonpositive_time_of_flight():
    with pytest.raises(ConvergenceError):
        rk4_validate(
            p_star_km=9000.0, e=0.3, nu1_rad=0.0, delta_theta_rad=1.0,
            tof_s=0.0, mu_km3_s2=Bodies.EARTH, r2_km=13556.0,
        )
