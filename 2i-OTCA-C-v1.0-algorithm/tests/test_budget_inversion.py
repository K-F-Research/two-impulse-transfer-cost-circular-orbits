"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Tests for `budget_inversion.py`, including an independent cross-check
against a plain bisection on phi that does not share any code path with
`invert_phi_for_budget`'s own safeguarded-Newton solver (other than the
underlying `core.find_p_star`/`core.cost`, which is exercised throughout
the rest of this package's own test suite).
"""
import math
import random

import pytest

from otca_c import core, hohmann
from otca_c.budget_inversion import invert_phi_for_budget
from otca_c.errors import ConvergenceError, DomainError
from otca_c.units import EARTH_MU_KM3_S2, circular_speed_km_s


def _independent_bisection_phi(R: float, budget_norm: float, tol: float = 1e-13) -> float:
    """A plain bisection on phi, sharing no logic with
    `budget_inversion._invert_phi_newton` other than `core.find_p_star`/
    `core.cost` themselves."""

    def g(phi: float) -> float:
        p = core.find_p_star(phi, R)
        return core.cost(p, phi, R) - budget_norm

    lo, hi = 1e-9, math.pi - 1e-9
    g_lo, g_hi = g(lo), g(hi)
    assert g_lo < 0 < g_hi
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        g_mid = g(mid)
        if abs(g_mid) < tol or (hi - lo) < 1e-14:
            return mid
        if g_mid < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@pytest.mark.parametrize("R", [1.05, 1.5, 2.0, 4.887, 11.324, 50.0])
@pytest.mark.parametrize("budget_factor", [1.01, 1.3, 2.0, 3.0])
def test_matches_independent_bisection(R, budget_factor):
    r1_km = 6678.0
    r2_km = r1_km * R
    vc1 = circular_speed_km_s(r1_km, EARTH_MU_KM3_S2)
    hohmann_km_s = hohmann.solve_hohmann(R).delta_v_total * vc1
    budget_km_s = hohmann_km_s * budget_factor

    result = invert_phi_for_budget(r1_km, r2_km, budget_km_s)

    budget_norm = budget_km_s / vc1
    phi_bisect = _independent_bisection_phi(R, budget_norm)

    assert math.isclose(result.phi_max_rad, phi_bisect, rel_tol=1e-8, abs_tol=1e-9)
    # Fewer iterations than plain bisection would need to reach the same
    # tolerance (bisection needs ~45 halvings of an O(1)-radian bracket to
    # reach 1e-13); demonstrates the Theorem-3-based solver's actual
    # benefit, not just its correctness.
    assert result.iterations < 20


def test_exact_hohmann_boundary_not_solved_numerically():
    r1_km, r2_km = 6678.0, 13356.0
    R = r2_km / r1_km
    vc1 = circular_speed_km_s(r1_km, EARTH_MU_KM3_S2)
    hohmann_km_s = hohmann.solve_hohmann(R).delta_v_total * vc1

    result = invert_phi_for_budget(r1_km, r2_km, hohmann_km_s)
    assert result.method == "exact_hohmann_boundary"
    assert result.phi_max_rad == 0.0
    assert result.delta_theta_min_deg == 180.0
    assert result.iterations == 0


def test_budget_below_hohmann_is_infeasible():
    r1_km, r2_km = 6678.0, 13356.0
    R = r2_km / r1_km
    vc1 = circular_speed_km_s(r1_km, EARTH_MU_KM3_S2)
    hohmann_km_s = hohmann.solve_hohmann(R).delta_v_total * vc1

    with pytest.raises(DomainError):
        invert_phi_for_budget(r1_km, r2_km, hohmann_km_s * 0.999)


def test_invalid_radii_raise_domain_error():
    with pytest.raises(DomainError):
        invert_phi_for_budget(6678.0, 6678.0, 10.0)  # r2 == r1
    with pytest.raises(DomainError):
        invert_phi_for_budget(6678.0, 1000.0, 10.0)  # r2 < r1
    with pytest.raises(DomainError):
        invert_phi_for_budget(-1.0, 6678.0, 10.0)  # negative radius


def test_regression_against_published_grid_round_trip():
    """For a handful of (R, delta_theta) pairs from the independently
    archived comparison grid, invert the already-known cost and recover
    the original delta_theta -- a round-trip check against externally
    archived data, not only against this module's own forward solver."""
    cases = [
        (1.5, 170.0), (1.5, 150.0),
        (2.0, 170.0), (2.0, 160.0),
        (5.0, 175.0), (10.0, 170.0),
    ]
    r1_km = 6678.0
    for R, delta_theta_deg in cases:
        r2_km = r1_km * R
        from otca_c.api import compute_transfer
        forward = compute_transfer(r1_km, r2_km, delta_theta_deg, validate=False)
        result = invert_phi_for_budget(r1_km, r2_km, forward.delta_v_total_km_s)
        assert math.isclose(
            result.delta_theta_min_deg, delta_theta_deg, rel_tol=1e-6, abs_tol=1e-6
        )


def test_extreme_near_degenerate_budget_does_not_hang():
    """A very large budget pushes phi_max close to the pi boundary (the
    degenerate antipodal-ish transfer); the solver must still terminate
    (possibly via the bisection fallback near the steep boundary) rather
    than raise ConvergenceError or loop indefinitely."""
    r1_km, r2_km = 6678.0, 13356.0
    result = invert_phi_for_budget(r1_km, r2_km, delta_v_budget_km_s=1000.0)
    assert result.delta_theta_min_deg < 1.0
    assert abs(result.residual_km_s) < 1e-2
