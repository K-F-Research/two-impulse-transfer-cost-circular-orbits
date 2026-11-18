"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Tests for the independent classical-mechanics cross-check (otca_c.geometry)
-- the project's stand-in for comparison against an external solver (no
maintained, installable Lambert-solver package was available in this
environment; see the folder README).
"""
import math

import pytest

from otca_c import core
from otca_c.geometry import independent_check


@pytest.mark.parametrize("phi_deg", [5, 30, 60, 90, 120, 150, 175])
@pytest.mark.parametrize("R", [1.1, 1.5, 2.0, 6.314, 20.0])
def test_geometry_self_consistent(phi_deg, R):
    phi = math.radians(phi_deg)
    sol = core.solve(phi, R)
    check = independent_check(sol.p_star, phi, R, sol.e, sol.S1, sol.S2)
    assert check.ok(tol=1e-6), (
        f"geometry residuals too large at phi={phi_deg} deg, R={R}: {check}"
    )
