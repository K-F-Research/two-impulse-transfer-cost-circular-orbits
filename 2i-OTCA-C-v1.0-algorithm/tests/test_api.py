"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Tests for the high-level physical-units API (otca_c.api).
"""
import math

import pytest

from otca_c import Bodies, DomainError, compute_transfer


def test_leo_to_geo_like_case_runs():
    result = compute_transfer(
        r1_km=6778.0,  # ~400 km altitude LEO, Earth radius 6378 km
        r2_km=13556.0,  # R=2 case
        delta_theta_deg=170.0,
        mu_km3_s2=Bodies.EARTH,
    )
    assert result.delta_v_total_km_s > 0
    assert result.geometry_check is not None
    assert result.geometry_check.ok(tol=1e-5)
    assert result.penalty_vs_hohmann_pct > 0  # any non-Hohmann angle costs more


def test_hohmann_case_exact():
    result = compute_transfer(
        r1_km=6778.0, r2_km=13556.0, delta_theta_deg=180.0, mu_km3_s2=Bodies.EARTH,
    )
    assert result.is_hohmann
    assert result.penalty_vs_hohmann_pct == pytest.approx(0.0, abs=1e-9)


def test_symmetry_short_vs_long_way():
    """delta_theta and 360-delta_theta must give the identical cost (the
    model's own evenness in phi -- see api.py docstring)."""
    r1_result = compute_transfer(6778.0, 13556.0, 170.0)
    r2_result = compute_transfer(6778.0, 13556.0, 190.0)
    assert r1_result.delta_v_total_km_s == pytest.approx(
        r2_result.delta_v_total_km_s, rel=1e-9
    )


def test_rejects_inbound_transfer():
    with pytest.raises(DomainError):
        compute_transfer(13556.0, 6778.0, 170.0)


def test_rejects_degenerate_angle():
    with pytest.raises(DomainError):
        compute_transfer(6778.0, 13556.0, 0.0)
    with pytest.raises(DomainError):
        compute_transfer(6778.0, 13556.0, 360.0)


def test_other_central_bodies():
    mars = compute_transfer(
        r1_km=3789.5 + 400.0, r2_km=2 * (3789.5 + 400.0), delta_theta_deg=170.0,
        mu_km3_s2=Bodies.MARS,
    )
    assert mars.delta_v_total_km_s > 0
    # Mars has weaker gravity than Earth: same normalized problem (same R,
    # same phi) must give a strictly smaller real delta-v.
    earth = compute_transfer(
        r1_km=3789.5 + 400.0, r2_km=2 * (3789.5 + 400.0), delta_theta_deg=170.0,
        mu_km3_s2=Bodies.EARTH,
    )
    assert mars.delta_v_total_km_s < earth.delta_v_total_km_s
