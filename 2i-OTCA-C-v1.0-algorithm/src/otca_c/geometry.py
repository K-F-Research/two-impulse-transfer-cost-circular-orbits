"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: none (independent cross-check via classical
two-body relations, deliberately not a re-derivation of this
repository's own results -- see the module description below).

Independent geometric/physical cross-check of a solved transfer, using only
textbook two-body relations (the conic section equation, conservation of
angular momentum, and the vis-viva equation) -- not the project's own
S1/S2/C definitions. This is the validation layer requested for this
component: before trusting a solved transfer operationally, re-derive its
basic physical consistency from first principles that do not depend on
`core.py`'s own algebra, rather than from a third-party solver (no
maintained, installable Lambert-solver package was available in this
environment -- see the folder README for the note on that choice).

Convention (consistent with `core.py` and
`docs/mathematical_conventions.md`): r1 = 1, mu = 1, R := r2/r1, phi := pi -
Delta_theta, Delta_theta being the true-anomaly angle actually traversed
from departure to arrival.
"""
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class GeometryCheck:
    """Residuals of the independent re-derivation; all should be ~0 (at
    floating-point precision) for a correctly solved transfer."""

    r1_residual: float
    r2_residual: float
    delta_v1_residual: float
    delta_v2_residual: float
    nu1_rad: float
    nu2_rad: float

    @property
    def max_abs_residual(self) -> float:
        return max(
            abs(self.r1_residual),
            abs(self.r2_residual),
            abs(self.delta_v1_residual),
            abs(self.delta_v2_residual),
        )

    def ok(self, tol: float = 1e-8) -> bool:
        return self.max_abs_residual < tol


def independent_check(p: float, phi: float, R: float, e: float, S1: float, S2: float) -> GeometryCheck:
    """Re-derive, from the conic equation r=p/(1+e*cos(nu)), conservation
    of angular momentum h=sqrt(mu*p) (mu=1 here), and the vis-viva equation
    v^2 = mu*(2/r - 1/a), the two radii and the two impulsive delta-v
    magnitudes -- entirely independently of core.py's S1/S2 formulas -- and
    compare them against the values core.py reports.

    This does NOT re-derive p* or e themselves (those still come from
    core.solve); it checks that, given p* and e, the resulting orbit is
    geometrically and physically consistent with r1=1, r2=R and the
    requested transfer angle, and that the delta-v magnitudes implied by
    elementary mechanics match S1, S2.
    """
    c, s = math.cos(phi), math.sin(phi)
    delta_theta = math.pi - phi  # true-anomaly angle actually traversed

    x = p - 1.0
    y = p / R - 1.0
    # e*sin(nu1), from y = e*cos(nu1+Delta_theta) expanded (see core.py's
    # module docstring derivation); x = e*cos(nu1) by definition of x.
    e_sin_nu1 = -(x * c + y) / s
    nu1 = math.atan2(e_sin_nu1, x)
    nu2 = nu1 + delta_theta

    r1_from_conic = p / (1 + e * math.cos(nu1))
    r2_from_conic = p / (1 + e * math.cos(nu2))

    a = p / (1 - e**2)
    h = p**0.5  # sqrt(mu*p), mu=1

    def speed_components(r: float, nu: float):
        v_theta = h / r
        v2_total = 2.0 / r - 1.0 / a  # vis-viva, mu=1
        v_radial_sq = v2_total - v_theta**2
        # numerically v_radial_sq can be a tiny negative at nu~0/pi due to
        # floating point; clamp rather than raise for a quantity that is
        # analytically exactly zero there.
        v_radial = math.copysign((max(v_radial_sq, 0.0)) ** 0.5, math.sin(nu))
        return v_theta, v_radial

    v_theta1, v_radial1 = speed_components(r1_from_conic, nu1)
    v_theta2, v_radial2 = speed_components(r2_from_conic, nu2)

    vc1 = 1.0
    vc2 = 1.0 / (R ** 0.5)

    delta_v1_sq_from_mechanics = (v_theta1 - vc1) ** 2 + v_radial1**2
    delta_v2_sq_from_mechanics = (v_theta2 - vc2) ** 2 + v_radial2**2

    return GeometryCheck(
        r1_residual=r1_from_conic - 1.0,
        r2_residual=r2_from_conic - R,
        delta_v1_residual=delta_v1_sq_from_mechanics - S1,
        delta_v2_residual=delta_v2_sq_from_mechanics - S2,
        nu1_rad=nu1,
        nu2_rad=nu2,
    )
