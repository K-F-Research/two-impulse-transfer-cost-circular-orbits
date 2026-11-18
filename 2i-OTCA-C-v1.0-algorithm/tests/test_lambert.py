"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Tests for the independent Lambert-problem cross-check (otca_c.lambert,
otca_c.derived.lambert_cross_check), across a grid of transfer angles
including both short-way and long-way (>180 deg) cases and several
central bodies, plus a dedicated regression test for the nu1/nu2 sign
convention for transfer angles exceeding 180 degrees (see api.py's
compute_transfer).
"""
import math

import pytest

from otca_c import Bodies, compute_transfer
from otca_c.lambert import solve_lambert
from otca_c.errors import ConvergenceError


CASES = [
    (6778.0, 13556.0, 170.0, Bodies.EARTH),
    (6778.0, 13556.0, 10.0, Bodies.EARTH),
    (6778.0, 13556.0, 190.0, Bodies.EARTH),  # long way, > 180 deg
    (6778.0, 13556.0, 350.0, Bodies.EARTH),  # near the 360 deg boundary
    (6778.0, 13556.0, 269.0, Bodies.EARTH),
    (6778.0, 13556.0, 5.0, Bodies.EARTH),
    (6778.0, 42164.0, 160.0, Bodies.EARTH),
    (3789.9, 24190.9, 170.0, Bodies.MARS),
    (1000.0, 7000.0, 90.0, Bodies.MOON),
]


@pytest.mark.parametrize("r1_km,r2_km,delta_theta_deg,mu", CASES)
def test_lambert_cross_check_matches_to_high_precision(r1_km, r2_km, delta_theta_deg, mu):
    result = compute_transfer(
        r1_km, r2_km, delta_theta_deg, mu_km3_s2=mu, lambert_check=True,
    )
    assert result.lambert_check is not None
    # A generous but still extremely tight tolerance: the near-360-degree
    # cases converge to a visibly larger (though still minuscule) residual
    # than the well-conditioned ones, due to the Lambert solver's own
    # numerical conditioning there (see api.py's compute_transfer
    # docstring) -- not a reduction in confidence in compute_transfer's
    # own result.
    assert result.lambert_check.ok(tol=1e-6)


def test_lambert_cross_check_hohmann_case():
    """Exactly dnu=180 deg is a genuine degeneracy of THIS solver's own
    A-formula (A -> 0), not of the underlying geometry -- compute_transfer
    returns lambert_check=None for the exact Hohmann case, the same way it
    returns geometry_check=None there (that case already has an exact
    closed form and never needed either cross-check)."""
    result = compute_transfer(
        6778.0, 13556.0, 180.0, mu_km3_s2=Bodies.EARTH, lambert_check=True,
    )
    assert result.is_hohmann
    assert result.lambert_check is None


def test_lambert_check_is_none_by_default():
    result = compute_transfer(6778.0, 13556.0, 170.0, mu_km3_s2=Bodies.EARTH)
    assert result.lambert_check is None


def test_long_way_transfer_has_longer_time_of_flight_than_short_way():
    """Regression test for the nu1/nu2 sign bug the Lambert cross-check
    uncovered: delta_theta=170 (short way) and delta_theta=190 (long way,
    the complementary arc of the identical ellipse) must NOT have the
    same time of flight, even though they have the identical cost
    (delta-v) by the model's own even-in-phi symmetry."""
    short_way = compute_transfer(6778.0, 13556.0, 170.0, mu_km3_s2=Bodies.EARTH)
    long_way = compute_transfer(6778.0, 13556.0, 190.0, mu_km3_s2=Bodies.EARTH)
    # Same ellipse (same p, e, a) -- the cost-minimization symmetry.
    assert short_way.p_star_km == pytest.approx(long_way.p_star_km, rel=1e-9)
    assert short_way.eccentricity == pytest.approx(long_way.eccentricity, rel=1e-9)
    # But the long-way arc (190 deg swept) must take noticeably longer
    # than the short-way arc (170 deg swept) of that same ellipse.
    assert long_way.time_of_flight_s > short_way.time_of_flight_s
    # And the two times of flight must sum to the full orbital period
    # (the two arcs partition the whole ellipse exactly once).
    assert short_way.time_of_flight_s + long_way.time_of_flight_s == pytest.approx(
        short_way.orbital_period_s, rel=1e-9
    )


def test_solve_lambert_matches_known_vis_viva_circular_limit():
    """Sanity check against a trivial closed-form case: a transfer with a
    tiny angle and p chosen so the orbit is nearly circular should give
    departure/arrival speeds close to the local circular speed."""
    mu = Bodies.EARTH
    r1_km = 7000.0
    r2_vec = (r1_km * math.cos(math.radians(1.0)), r1_km * math.sin(math.radians(1.0)))
    # time of flight for a near-circular orbit sweeping 1 degree:
    circular_speed = (mu / r1_km) ** 0.5
    period = 2 * math.pi * r1_km / circular_speed
    tof = period * (1.0 / 360.0)
    sol = solve_lambert((r1_km, 0.0), r2_vec, tof, mu, prograde=True)
    speed1 = (sol.v1_km_s[0] ** 2 + sol.v1_km_s[1] ** 2) ** 0.5
    assert speed1 == pytest.approx(circular_speed, rel=1e-3)


def test_solve_lambert_raises_on_degenerate_zero_angle():
    with pytest.raises(ConvergenceError):
        solve_lambert((7000.0, 0.0), (7000.0, 0.0), 3600.0, Bodies.EARTH)


# --- Regression guard for the clamp-biased-finite-difference defect -----
#
# An independent review of this module found that the former Newton +
# finite-difference root search could jump to the z=4*pi^2 clamp boundary
# on its very first step for transfer angles close to 360 deg, then crawl
# back down over dozens of iterations -- confirmed by tracing the actual
# z-sequence for delta_theta=350 deg, which needed 47-48 iterations (vs.
# 4-7 for a well-conditioned case) and spent ~9 of them pinned exactly at
# the clamp. The fix replaces Newton+finite-difference with an adaptive
# bracket + a bracketing root-finder (tof(z) is monotonic on the valid
# range -- verified numerically across every case below), the same
# pattern already used for p* in core.py. This test pins a performance
# ceiling that would fail under the old implementation and must keep
# passing under the new one, as a permanent guard against regressing back
# to clamp-biased Newton.
@pytest.mark.parametrize("delta_theta_deg", [350.0, 355.0, 358.0, 359.0, 359.9, 359.99, 269.0])
def test_lambert_near_360_boundary_converges_without_excessive_iterations(delta_theta_deg):
    result = compute_transfer(6778.0, 13556.0, delta_theta_deg, mu_km3_s2=Bodies.EARTH)
    dtheta_rad = math.radians(delta_theta_deg)
    r1_vec = (6778.0, 0.0)
    r2_vec = (13556.0 * math.cos(dtheta_rad), 13556.0 * math.sin(dtheta_rad))
    sol = solve_lambert(r1_vec, r2_vec, result.time_of_flight_s, Bodies.EARTH)
    # A bracketing root-finder on a monotonic function should need very
    # few iterations (tens, not dozens) regardless of how close the
    # geometry is to the 360-degree boundary -- the old Newton+clamp
    # implementation needed 47-48 here.
    assert sol.iterations < 15
    # A low iteration count alone does not prove the root is accurate
    # (flagged during independent review) -- also check the actual
    # residual of the time-of-flight equation, the finiteness of the
    # recovered velocities, and agreement with the primary solver's own
    # (p, e).
    from otca_c.lambert import _stumpff_c, _stumpff_s
    mu = Bodies.EARTH
    C = _stumpff_c(sol.z)
    S = _stumpff_s(sol.z)
    r1, r2 = math.hypot(*r1_vec), math.hypot(*r2_vec)
    cos_dnu = max(-1.0, min(1.0, (r1_vec[0] * r2_vec[0] + r1_vec[1] * r2_vec[1]) / (r1 * r2)))
    dnu = math.acos(cos_dnu)
    if r1_vec[0] * r2_vec[1] - r1_vec[1] * r2_vec[0] < 0:
        dnu = 2 * math.pi - dnu
    A = (2.0 * r1 * r2) ** 0.5 * math.cos(dnu / 2.0)
    y = r1 + r2 + A * (sol.z * S - 1.0) / (C**0.5)
    chi = (y / C) ** 0.5
    t_at_root = (chi**3 * S + A * y**0.5) / mu**0.5
    assert abs(t_at_root - result.time_of_flight_s) < 1e-6
    assert all(math.isfinite(c) for c in sol.v1_km_s + sol.v2_km_s)
    assert sol.p_km == pytest.approx(result.p_star_km, rel=1e-6)
    assert sol.e == pytest.approx(result.eccentricity, abs=1e-6)
