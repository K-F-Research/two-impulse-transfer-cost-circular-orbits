"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Additional edge-case and integration tests for `otca_c.surrogate`, written
to widen coverage beyond the two existing test files
(`test_surrogate.py`: fixed reference values; `test_surrogate_properties.py`:
scaling/metamorphic invariants). These focus on boundaries, package-level
exposure, and error-path inputs not already exercised elsewhere, so as not
to duplicate what those two files already check.
"""
import math

import pytest

import otca_c
from otca_c.errors import DomainError
from otca_c.surrogate import (
    MAX_VALIDATED_PHI_DEG,
    SurrogateComparison,
    compare_to_exact,
    curvature_coefficient_km_s,
    estimate_cost,
)


# ---------------------------------------------------------------------------
# Package-level exposure: surrogate must be reachable from `otca_c` directly,
# not only from `otca_c.surrogate`, now that it has been merged into the
# public API (see __init__.py).
# ---------------------------------------------------------------------------
def test_surrogate_reexported_from_package_root() -> None:
    assert otca_c.estimate_cost is estimate_cost
    assert otca_c.compare_to_exact is compare_to_exact
    assert otca_c.curvature_coefficient_km_s is curvature_coefficient_km_s
    assert otca_c.MAX_VALIDATED_PHI_DEG == MAX_VALIDATED_PHI_DEG


# ---------------------------------------------------------------------------
# Exact Hohmann point (phi = 0): the surrogate's quadratic term vanishes
# identically there, so the estimate must equal the closed-form Hohmann
# total exactly (up to floating-point noise), with zero curvature
# contribution -- a degenerate but legitimate corner of the domain.
# ---------------------------------------------------------------------------
def test_estimate_at_exact_hohmann_equals_hohmann_reference() -> None:
    r1_km, r2_km = 6678.0, 10017.0  # R = 1.5
    result = estimate_cost(r1_km, r2_km, delta_theta_deg=180.0, mu_km3_s2=otca_c.Bodies.EARTH)
    assert result.phi_deg == pytest.approx(0.0, abs=1e-9)
    assert result.delta_v_total_estimate_km_s == pytest.approx(
        result.hohmann_delta_v_total_km_s, rel=1e-12
    )
    assert result.within_measured_domain is True


# ---------------------------------------------------------------------------
# Boundary just past MAX_VALIDATED_PHI_DEG: must flip to False, immediately
# next to the exact-boundary case already covered by test_surrogate.py's
# test_domain_flag_and_guard_behaviour (which only checks exactly at the
# threshold, not just beyond it).
# ---------------------------------------------------------------------------
def test_domain_flag_just_past_threshold() -> None:
    r1_km, r2_km = 6678.0, 10017.0
    just_inside_delta_theta = 180.0 - (MAX_VALIDATED_PHI_DEG - 1e-6)
    just_outside_delta_theta = 180.0 - (MAX_VALIDATED_PHI_DEG + 1e-6)
    inside = estimate_cost(r1_km, r2_km, just_inside_delta_theta, otca_c.Bodies.EARTH)
    outside = estimate_cost(r1_km, r2_km, just_outside_delta_theta, otca_c.Bodies.EARTH)
    assert inside.within_measured_domain is True
    assert outside.within_measured_domain is False


# ---------------------------------------------------------------------------
# Angle wrap-around: delta_theta_deg outside [0,360) must fold to the same
# physical angle as its canonical representative (e.g. 540 deg == 180 deg,
# -10 deg == 350 deg), since estimate_cost uses `% 360.0` internally. Not
# previously tested directly.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "raw_delta_theta, canonical_delta_theta",
    [(540.0, 180.0), (-10.0, 350.0), (720.0 + 170.0, 170.0), (-190.0, 170.0)],
)
def test_angle_wraparound_matches_canonical(raw_delta_theta: float, canonical_delta_theta: float) -> None:
    r1_km, r2_km = 6678.0, 10017.0
    a = estimate_cost(r1_km, r2_km, raw_delta_theta, otca_c.Bodies.EARTH)
    b = estimate_cost(r1_km, r2_km, canonical_delta_theta, otca_c.Bodies.EARTH)
    assert a.delta_v_total_estimate_km_s == pytest.approx(b.delta_v_total_estimate_km_s, rel=1e-12)
    assert a.phi_deg == pytest.approx(b.phi_deg, abs=1e-9)


# ---------------------------------------------------------------------------
# compare_to_exact's relative_error must be non-negative and finite across a
# broader grid of (R, phi) than the single-R check already in
# test_surrogate.py, and SurrogateComparison's fields must be internally
# consistent (estimate.delta_v_total_estimate_km_s matches what went into
# the stored relative_error).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("R", [1.05, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0])
@pytest.mark.parametrize("phi_deg", [1.0, 5.0, 15.0, 29.0])
def test_compare_to_exact_grid_well_formed(R: float, phi_deg: float) -> None:
    r1_km = 6678.0
    r2_km = r1_km * R
    delta_theta_deg = 180.0 - phi_deg
    comparison = compare_to_exact(r1_km, r2_km, delta_theta_deg, otca_c.Bodies.EARTH)
    assert isinstance(comparison, SurrogateComparison)
    assert math.isfinite(comparison.relative_error)
    assert comparison.relative_error >= 0.0
    recomputed_err = abs(
        comparison.estimate.delta_v_total_estimate_km_s - comparison.exact_delta_v_total_km_s
    ) / comparison.exact_delta_v_total_km_s
    assert comparison.relative_error == pytest.approx(recomputed_err, rel=1e-12)


# ---------------------------------------------------------------------------
# Error paths not already covered: NaN/inf inputs, and r1_km == r2_km
# (R == 1 exactly, a degenerate circular case distinct from the already
# tested r2 < r1 case).
# ---------------------------------------------------------------------------
def test_estimate_cost_rejects_equal_radii() -> None:
    with pytest.raises(DomainError):
        estimate_cost(6678.0, 6678.0, 170.0)


def test_estimate_cost_rejects_non_finite_inputs() -> None:
    with pytest.raises(DomainError):
        estimate_cost(float("nan"), 10017.0, 170.0)
    with pytest.raises(DomainError):
        estimate_cost(6678.0, float("inf"), 170.0)


def test_curvature_coefficient_rejects_equal_radii() -> None:
    with pytest.raises(DomainError):
        curvature_coefficient_km_s(6678.0, 6678.0, otca_c.Bodies.EARTH)


# ---------------------------------------------------------------------------
# Consistency between estimate_cost and curvature_coefficient_km_s: the K
# embedded inside a SurrogateEstimate must be exactly what the standalone
# function returns for the same inputs (no duplicated/drifted computation
# path between the two public entry points).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("R", [1.1, 1.5, 2.0, 4.0, 8.0])
def test_estimate_K_matches_standalone_curvature_function(R: float) -> None:
    r1_km = 6678.0
    r2_km = r1_km * R
    est = estimate_cost(r1_km, r2_km, delta_theta_deg=175.0, mu_km3_s2=otca_c.Bodies.EARTH)
    K_standalone = curvature_coefficient_km_s(r1_km, r2_km, otca_c.Bodies.EARTH)
    assert est.K_km_s == pytest.approx(K_standalone, rel=1e-12)
