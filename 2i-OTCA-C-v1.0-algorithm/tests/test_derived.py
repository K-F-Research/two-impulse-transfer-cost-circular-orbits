"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Tests for the derived quantities in otca_c.derived, and their wiring into
compute_transfer (otca_c.api).
"""
import math

import pytest

from otca_c import Bodies, BodyRadii, compute_transfer
from otca_c.derived import (
    ECCENTRICITY_CONDITIONING_THRESHOLD,
    altitude_info,
    apoapsis_radius_km,
    eccentricity_conditioning_warning,
    orbital_period_s,
    periapsis_radius_km,
    propellant_info,
    time_of_flight_s,
)


def test_periapsis_apoapsis_radii_match_conic_definition():
    p, e = 10000.0, 0.3
    r_p = periapsis_radius_km(p, e)
    r_a = apoapsis_radius_km(p, e)
    a = p / (1 - e**2)
    assert r_p == pytest.approx(a * (1 - e), rel=1e-12)
    assert r_a == pytest.approx(a * (1 + e), rel=1e-12)


def test_orbital_period_matches_keplers_third_law():
    a_km = 10000.0
    mu = Bodies.EARTH
    T = orbital_period_s(a_km, mu)
    assert T == pytest.approx(2 * math.pi * (a_km**3 / mu) ** 0.5, rel=1e-12)


def test_time_of_flight_hohmann_is_half_period():
    """For an exact Hohmann transfer (periapsis to apoapsis), time of
    flight must be exactly half the orbital period -- a basic sanity
    check independent of the general formula's unwrapping logic."""
    a_km, e, mu = 10167.0, 1.0 / 3.0, Bodies.EARTH
    tof = time_of_flight_s(a_km, e, 0.0, math.pi, mu)
    period = orbital_period_s(a_km, mu)
    assert tof == pytest.approx(period / 2.0, rel=1e-9)


def test_time_of_flight_full_orbit_equals_period():
    a_km, e, mu = 10167.0, 1.0 / 3.0, Bodies.EARTH
    tof = time_of_flight_s(a_km, e, 0.0, 2 * math.pi, mu)
    period = orbital_period_s(a_km, mu)
    assert tof == pytest.approx(period, rel=1e-9)


def test_time_of_flight_monotonic_in_swept_angle():
    """Sweeping a larger angle must always take longer (no wraparound
    sign error) -- regression guard for the unwrapping logic."""
    a_km, e, mu = 10167.0, 1.0 / 3.0, Bodies.EARTH
    nu1 = 0.3
    tof_small = time_of_flight_s(a_km, e, nu1, nu1 + 1.0, mu)
    tof_large = time_of_flight_s(a_km, e, nu1, nu1 + 2.5, mu)
    assert 0 < tof_small < tof_large


def test_altitude_info_flags_intersecting_orbit():
    info = altitude_info(r_periapsis_km=6000.0, r_apoapsis_km=8000.0, body_radius_km=6378.137)
    assert info.warning is not None
    assert info.periapsis_altitude_km < 0


def test_altitude_info_normal_case_no_warning():
    info = altitude_info(r_periapsis_km=6778.0, r_apoapsis_km=13556.0, body_radius_km=6378.137)
    assert info.warning is None
    assert info.periapsis_altitude_km == pytest.approx(6778.0 - 6378.137)
    assert info.apoapsis_altitude_km == pytest.approx(13556.0 - 6378.137)


def test_propellant_info_tsiolkovsky():
    # Classic sanity check: delta_v = Isp*g0*ln(1/(1-mf)) inverted.
    isp = 300.0
    g0 = 9.80665
    delta_v_km_s = 2.0
    info = propellant_info(delta_v_km_s, isp)
    recovered_delta_v = -isp * g0 * math.log(1 - info.mass_fraction) / 1000.0
    assert recovered_delta_v == pytest.approx(delta_v_km_s, rel=1e-9)


def test_propellant_info_rejects_bad_isp():
    with pytest.raises(ValueError):
        propellant_info(1.0, 0.0)


def test_apoapsis_radius_rejects_hyperbolic():
    with pytest.raises(ValueError):
        apoapsis_radius_km(10000.0, 1.0)


# --- Wiring into compute_transfer --------------------------------------


def test_compute_transfer_unconditional_derived_fields_always_present():
    result = compute_transfer(6778.0, 13556.0, 170.0, mu_km3_s2=Bodies.EARTH)
    assert result.time_of_flight_s > 0
    assert result.orbital_period_s > result.time_of_flight_s
    assert result.r_periapsis_km < result.r_apoapsis_km
    assert result.altitude_info is None
    assert result.propellant_info_departure is None
    assert result.propellant_info_arrival is None


def test_compute_transfer_altitude_only_when_body_radius_supplied():
    result = compute_transfer(
        6778.0, 13556.0, 170.0, mu_km3_s2=Bodies.EARTH, body_radius_km=BodyRadii.EARTH,
    )
    assert result.altitude_info is not None
    assert result.altitude_info.warning is None


def test_compute_transfer_propellant_only_when_isp_supplied():
    result = compute_transfer(6778.0, 13556.0, 170.0, mu_km3_s2=Bodies.EARTH, isp_s=320.0)
    assert result.propellant_info_departure is not None
    assert result.propellant_info_arrival is not None
    assert 0 < result.propellant_info_departure.mass_fraction < 1
    assert 0 < result.propellant_info_arrival.mass_fraction < 1


def test_compute_transfer_hohmann_time_of_flight_is_half_period():
    result = compute_transfer(6778.0, 13556.0, 180.0, mu_km3_s2=Bodies.EARTH)
    assert result.is_hohmann
    assert result.time_of_flight_s == pytest.approx(result.orbital_period_s / 2.0, rel=1e-9)


def test_compute_transfer_altitude_warns_for_too_low_orbit():
    # r1 just above Earth's surface but the transfer's periapsis (==r1
    # here, since this is the departure radius itself) at a body radius
    # that is deliberately larger than r1 -- forces the warning branch.
    result = compute_transfer(
        6500.0, 13000.0, 170.0, mu_km3_s2=Bodies.EARTH, body_radius_km=6600.0,
    )
    assert result.altitude_info.warning is not None


# --- Eccentricity numerical-conditioning warnings (e -> 1) -----------------
#
# These are NOT physical-validity errors (e<1 genuinely holds -- the orbit
# is a real, physically valid ellipse); they flag that a FEW downstream
# quantities that divide by (1-e) or (1-e^2) -- r_apoapsis_km,
# semi_major_axis_km, orbital_period_s, time_of_flight_s -- can carry a
# much larger relative floating-point error than the rest of the output
# when 1-e is extremely small, per the repository's agreed three-level
# structure (physical validity / numerical validity / per-output
# reliability), structured by grandeur rather than a single global
# eccentricity threshold.

def test_eccentricity_conditioning_warning_does_not_fire_for_ordinary_cases():
    for e in (0.0, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999):
        assert eccentricity_conditioning_warning(e) is None


def test_eccentricity_conditioning_warning_fires_just_past_the_threshold():
    e = 1.0 - ECCENTRICITY_CONDITIONING_THRESHOLD / 10.0
    warning = eccentricity_conditioning_warning(e)
    assert warning is not None
    assert "r_apoapsis_km" in warning
    assert "semi_major_axis_km" in warning
    assert "orbital_period_s" in warning
    assert "time_of_flight_s" in warning
    # The quantities NOT affected by this specific conditioning issue must
    # be explicitly named as unaffected, not just silently absent --
    # a reader skimming this warning should not have to guess whether
    # delta-v and periapsis are also suspect.
    assert "delta_v1_km_s" in warning
    assert "r_periapsis_km" in warning


def test_eccentricity_conditioning_warning_is_none_outside_the_valid_e_domain():
    # Not this function's job to validate e's domain (0<=e<1) -- that is
    # already enforced elsewhere (e.g. time_of_flight_s raises ValueError
    # for e>=1); this function must simply decline to comment when e is
    # outside its own purview, rather than raise or fabricate a warning.
    assert eccentricity_conditioning_warning(-0.1) is None
    assert eccentricity_conditioning_warning(1.0) is None
    assert eccentricity_conditioning_warning(1.5) is None


def test_compute_transfer_has_no_numerical_warnings_for_an_ordinary_case():
    result = compute_transfer(6778.0, 13556.0, 170.0, mu_km3_s2=Bodies.EARTH)
    assert result.numerical_warnings == ()


def test_compute_transfer_warns_near_the_359_999_degree_boundary():
    # The known Delta_theta=359.999 deg, R=2.0 case (see
    # test_core.test_find_p_star_near_360_degree_boundary_regression) has
    # e extremely close to 1 (1-e of order 1e-10, well past this
    # threshold of 1e-6) -- compute_transfer must surface the warning
    # here, not just core.find_p_star's own successful convergence.
    result = compute_transfer(6778.0, 13556.0 * 1.0, 359.999, mu_km3_s2=Bodies.EARTH)
    assert 1.0 - result.eccentricity < ECCENTRICITY_CONDITIONING_THRESHOLD
    assert len(result.numerical_warnings) == 1
    assert "time_of_flight_s" in result.numerical_warnings[0]
