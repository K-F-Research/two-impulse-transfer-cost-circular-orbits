"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Checks that the general solver (otca_c.core) smoothly approaches the exact
closed-form Hohmann transfer (otca_c.hohmann) as phi -> 0, and that
Theorem 2 (circular/hohmann_local_minimum/: Hohmann is a strict local
minimum of the cost) is reflected numerically.
"""
import math

import pytest

from otca_c import core, hohmann


@pytest.mark.parametrize("R", [1.1, 1.5, 2.0, 6.314, 50.0])
def test_general_solver_converges_to_hohmann_as_phi_to_0(R):
    h = hohmann.solve_hohmann(R)
    for phi_deg in (1.0, 0.1, 0.01):
        sol = core.solve(math.radians(phi_deg), R)
        assert sol.C == pytest.approx(h.delta_v_total, abs=5e-3)


@pytest.mark.parametrize("R", [1.5, 2.0, 6.314])
def test_hohmann_is_strict_local_minimum(R):
    """Theorem 1: moving slightly away from phi=0 (either direction, via
    the model's symmetry in phi) must strictly increase the optimized
    cost."""
    h = hohmann.solve_hohmann(R)
    for phi_deg in (0.5, 2.0, 5.0):
        sol = core.solve(math.radians(phi_deg), R)
        assert sol.C > h.delta_v_total


def test_hohmann_rejects_r_not_greater_than_1():
    with pytest.raises(ValueError):
        hohmann.solve_hohmann(1.0)
    with pytest.raises(ValueError):
        hohmann.solve_hohmann(0.5)
