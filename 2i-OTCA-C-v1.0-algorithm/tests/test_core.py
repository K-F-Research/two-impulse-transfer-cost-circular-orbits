"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Unit tests for the normalized-unit physics core (otca_c.core).
"""
import math

import pytest

from otca_c import core
from otca_c.errors import ConvergenceError


def test_p_star_matches_stationary_condition():
    """p* must make C_p vanish, for a spread of (phi, R)."""
    for phi_deg in (10, 45, 90, 135, 170):
        for R in (1.2, 2.0, 6.314, 50.0):
            phi = math.radians(phi_deg)
            p_star = core.find_p_star(phi, R)
            assert abs(core.cost_derivative(p_star, phi, R)) < 1e-8


def test_p_star_is_a_local_minimum():
    """Perturbing p away from p* in either direction must increase the
    cost (second-order check, central-difference)."""
    phi, R = math.radians(60.0), 2.0
    p_star = core.find_p_star(phi, R)
    c0 = core.cost(p_star, phi, R)
    h = 1e-5
    c_minus = core.cost(p_star - h, phi, R)
    c_plus = core.cost(p_star + h, phi, R)
    assert c_minus > c0
    assert c_plus > c0


def test_r2_curvature_sign_matches_theorem4():
    """Independent regression against the certified result of Theorem 4
    (circular/inflection_uniqueness_R2/): for R=2, the curvature of V_2
    changes sign at phi_c approx 164.3710517 deg -- i.e. just below that
    angle the stationary point is in one curvature regime, just above it,
    the other. We check this indirectly: the cost function's second
    derivative with respect to phi (finite difference through p*(phi))
    changes sign across phi_c, consistent with the certified value.
    This ties this operational implementation back to the repository's
    own certified mathematics, rather than leaving it a free-standing,
    unverified reimplementation.
    """
    R = 2.0
    phi_c_deg = 164.3710517

    def v2(phi_deg):
        phi = math.radians(phi_deg)
        p_star = core.find_p_star(phi, R)
        return core.cost(p_star, phi, R)

    h = 0.01
    below = phi_c_deg - 2.0
    above = phi_c_deg + 2.0
    d2_below = (v2(below - h) - 2 * v2(below) + v2(below + h)) / h**2
    d2_above = (v2(above - h) - 2 * v2(above) + v2(above + h)) / h**2
    assert d2_below * d2_above < 0, (
        "expected a curvature sign change near the Theorem 4 certified "
        f"phi_c={phi_c_deg} deg, got d2_below={d2_below!r}, d2_above={d2_above!r}"
    )


def test_find_p_star_self_corrects_a_starting_bracket_that_misses_the_root():
    """find_p_star's starting p_lo/p_hi are only a SEED, not a fixed
    bracket: a starting bracket that does not itself contain the root
    must now be adaptively contracted/expanded until it does, rather than
    raising just because the initial guess happened to miss -- this is a
    deliberate behavior change from an earlier version of this function
    that treated any non-bracketing p_lo/p_hi as a hard failure."""
    phi, R = math.radians(90.0), 2.0
    # [1e-3, 1e-2] does not contain the true root at phi=90deg, R=2 -- the
    # adaptive search must still find it rather than raise.
    p_star = core.find_p_star(phi, R, p_lo=1e-3, p_hi=1e-2)
    # Cross-check against the default-bracket call, which does happen to
    # contain the root directly.
    assert p_star == pytest.approx(core.find_p_star(phi, R), rel=1e-9)


def test_find_p_star_raises_on_genuine_domain_boundary():
    """phi=0 exactly is the Hohmann/domain boundary (sin(phi)=0, the
    model's own 1/sin(phi)**2 singularity) -- not reachable through
    api.compute_transfer (routed to the exact Hohmann closed form
    instead), but find_p_star itself must still fail loudly rather than
    silently expand its bracket forever or return a nonsensical value."""
    with pytest.raises(ConvergenceError):
        core.find_p_star(0.0, 2.0)


def test_find_p_star_near_360_degree_boundary_regression():
    """Regression test for a real failure found during independent
    review: delta_theta=359.999 deg (R=2.0) corresponds to
    phi=3.1415752... (extremely close to pi), whose true p* (~3e-10) was
    below the previous fixed p_lo=1e-9 -- C_p(1e-9) was already positive
    there (the wrong sign per Lemma 1's p->0+ prediction), causing a
    spurious ConvergenceError even though the root genuinely exists and
    is unique. The adaptive bracket search must find it."""
    phi = math.pi - math.radians(359.999)  # signed phi, matches api.py's own convention
    p_star = core.find_p_star(phi, 2.0)
    assert 0 < p_star < 1e-8
    # C_p is extremely steep in this near-boundary region (it swings from
    # ~-1e7 to ~+1e8 across a change in p of a few times 1e-10), so an
    # absolute tolerance on C_p(p_star) itself is not a meaningful check
    # here -- instead verify directly that p_star is actually a local
    # minimum of the cost (lower than its immediate neighbors), which is
    # what "being a root of C_p" is supposed to mean.
    c_star = core.cost(p_star, phi, 2.0)
    assert c_star <= core.cost(p_star * 1.001, phi, 2.0)
    assert c_star <= core.cost(p_star * 0.999, phi, 2.0)


@pytest.mark.parametrize("R", [1.01, 1.5, 2.0, 6.314, 10.0, 100.0])
def test_solve_runs_across_radius_ratios(R):
    sol = core.solve(math.radians(170.0), R)
    assert sol.p_star > 0
    assert 0 < sol.e < 1
    assert sol.a > 0
    assert sol.C > 0
