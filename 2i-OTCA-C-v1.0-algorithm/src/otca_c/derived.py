"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: none (classical two-body relations -- Kepler's
equation, the Tsiolkovsky rocket equation -- independent of this
repository's own proven results; `lambert_cross_check` below delegates
to `lambert.py`, itself independent of them too).

Derived quantities beyond the core cost/geometry solve, split strictly into
two groups:

1. Quantities computable with full rigor from data the solver already has
   (the transfer orbit's own a, e, p, nu1, nu2, and the central body's mu,
   which the caller always supplies): time of flight (Kepler's equation),
   the transfer orbit's own period, and its periapsis/apoapsis radii. These
   are always computed, unconditionally.

2. Quantities that need a genuinely separate physical input the orbit
   solve itself never implies: an altitude needs the central body's
   physical radius (not derivable from mu); a propellant mass fraction
   needs the engine's specific impulse (Isp). Both are optional
   parameters here, returned as None when not supplied, rather than
   guessed or defaulted -- see `units.BodyRadii` and `compute_transfer`'s
   docstring in `api.py`.
"""
import math
from dataclasses import dataclass
from typing import Optional

from .lambert import solve_lambert


# Numerical (not physical) safety threshold: below 1-e this small, a few
# downstream quantities built from 1/(1-e) -- r_apoapsis_km = p/(1-e),
# semi_major_axis_km = p/(1-e^2), and anything derived from a (the
# orbital period, the time of flight) -- can amplify an ordinary
# floating-point-level error in e into a visibly larger relative error in
# the output. This is a conditioning threshold of THIS API's floating-
# point implementation, not a claim that the physical model breaks down
# there: the orbit remains a perfectly valid ellipse (e<1 genuinely
# holds) above and below this value, and the cost/delta-v values
# (computed without dividing by 1-e anywhere) are not affected by this at
# all -- see `api.py`'s own per-quantity documentation of which fields
# this warning actually covers. The exact value below is a deliberately
# conservative engineering choice, not derived from a rigorous backward-
# error analysis of the whole solve chain (phi/R -> p* -> e); independent
# review flagged such an analysis as a possible future improvement.
ECCENTRICITY_CONDITIONING_THRESHOLD = 1e-6


def eccentricity_conditioning_warning(e: float) -> Optional[str]:
    """None if e is comfortably below 1 (1-e well above the threshold);
    otherwise a message naming the threshold and the affected quantities,
    for a caller to attach to its result rather than silently trust a
    numerically fragile apoapsis/period/time-of-flight figure."""
    if not (0.0 <= e < 1.0):
        return None  # not this function's job to validate e's domain
    if (1.0 - e) < ECCENTRICITY_CONDITIONING_THRESHOLD:
        return (
            f"eccentricity e={e!r} is extremely close to 1 (1-e="
            f"{1.0 - e:.3e}, below this API's numerical conditioning "
            f"threshold of {ECCENTRICITY_CONDITIONING_THRESHOLD!r}): "
            "r_apoapsis_km, semi_major_axis_km, orbital_period_s and "
            "time_of_flight_s divide by (1-e) or (1-e^2) and can carry a "
            "much larger relative error here than elsewhere in the "
            "domain, even though the underlying orbit is a genuine, "
            "physically valid ellipse. delta_v1_km_s, delta_v2_km_s, "
            "delta_v_total_km_s and r_periapsis_km do not divide by "
            "(1-e) and are not affected by this warning."
        )
    return None


# Standard gravity, used only for the Tsiolkovsky rocket equation (Isp is
# conventionally quoted in seconds, which already bakes in g0 = 9.80665
# m/s^2 by definition of "specific impulse in seconds").
STANDARD_GRAVITY_M_S2 = 9.80665


def periapsis_radius_km(p_km: float, e: float) -> float:
    """r_p = p / (1 + e), exact for any conic; no extra input needed."""
    return p_km / (1.0 + e)


def apoapsis_radius_km(p_km: float, e: float) -> float:
    """r_a = p / (1 - e); only defined (finite) for a bound ellipse, e<1,
    which is always the case here since the transfer orbit's two apsides
    coincide with two distinct finite circular-orbit radii."""
    if e >= 1.0:
        raise ValueError(f"apoapsis radius undefined for e={e!r} >= 1 (not an ellipse)")
    return p_km / (1.0 - e)


def orbital_period_s(a_km: float, mu_km3_s2: float) -> float:
    """Kepler's third law: T = 2*pi*sqrt(a^3/mu). This is the period of the
    full transfer ellipse if the spacecraft were left on it indefinitely --
    not the duration of the transfer itself (see `time_of_flight_s`)."""
    if a_km <= 0:
        raise ValueError(f"semi-major axis must be positive, got {a_km!r} km")
    return 2.0 * math.pi * (a_km**3 / mu_km3_s2) ** 0.5


def _true_to_eccentric(nu_rad: float, e: float) -> float:
    """Eccentric anomaly E from true anomaly nu, for 0<=e<1, via the
    standard half-angle form E = 2*atan2(sqrt(1-e)*sin(nu/2),
    sqrt(1+e)*cos(nu/2)). atan2 keeps this well-defined for any nu, and
    returns a value in (-pi, pi]; the caller unwraps across branch cuts
    (see `time_of_flight_s`) rather than relying on this function's raw
    branch."""
    return 2.0 * math.atan2(
        (1.0 - e) ** 0.5 * math.sin(nu_rad / 2.0),
        (1.0 + e) ** 0.5 * math.cos(nu_rad / 2.0),
    )


def time_of_flight_s(a_km: float, e: float, nu1_rad: float, nu2_rad: float, mu_km3_s2: float) -> float:
    """Time to travel from true anomaly nu1 to nu2 (nu2 > nu1, both in
    radians, nu2 - nu1 the transfer angle actually swept -- these are NOT
    folded into [0, 2*pi) by the caller, so this function unwraps the
    corresponding eccentric/mean anomalies onto the same continuous branch
    before differencing, rather than risking a silent 2*pi wraparound
    error for a transfer angle that exceeds pi.

    Uses Kepler's equation M = E - e*sin(E) and mean motion n =
    sqrt(mu/a^3): t = (M2 - M1) / n.
    """
    if not (0.0 <= e < 1.0):
        raise ValueError(f"time_of_flight_s requires an ellipse, 0<=e<1; got e={e!r}")
    n = (mu_km3_s2 / a_km**3) ** 0.5  # mean motion, rad/s

    def mean_anomaly_unwrapped(nu: float) -> float:
        # Number of full revolutions nu has swept past, tracked explicitly
        # so E (and hence M) is placed on the *same* continuous branch as
        # nu itself, instead of being silently folded back into (-pi, pi].
        k = math.floor((nu + math.pi) / (2.0 * math.pi))
        nu_folded = nu - k * 2.0 * math.pi
        E = _true_to_eccentric(nu_folded, e)
        M = E - e * math.sin(E)
        return M + k * 2.0 * math.pi

    M1 = mean_anomaly_unwrapped(nu1_rad)
    M2 = mean_anomaly_unwrapped(nu2_rad)
    return (M2 - M1) / n


@dataclass(frozen=True)
class AltitudeInfo:
    """Periapsis/apoapsis altitude above the central body's *surface*,
    requiring the body's physical radius -- a second, independent input
    beyond mu (see `units.BodyRadii`). `warning` is set (non-None) if the
    transfer orbit's periapsis radius is at or below the body radius,
    meaning the idealized point-mass trajectory would intersect the body
    -- a physically impossible transfer that this component's point-mass
    model cannot itself detect without this extra input."""

    body_radius_km: float
    periapsis_altitude_km: float
    apoapsis_altitude_km: float
    warning: Optional[str]


def altitude_info(r_periapsis_km: float, r_apoapsis_km: float, body_radius_km: float) -> AltitudeInfo:
    if body_radius_km <= 0:
        raise ValueError(f"body_radius_km must be positive, got {body_radius_km!r}")
    peri_alt = r_periapsis_km - body_radius_km
    apo_alt = r_apoapsis_km - body_radius_km
    warning = None
    if peri_alt <= 0:
        warning = (
            f"periapsis radius ({r_periapsis_km:.3f} km) is at or below the "
            f"supplied body radius ({body_radius_km:.3f} km): this idealized "
            "point-mass transfer orbit would intersect the central body and "
            "is not physically flyable as computed."
        )
    return AltitudeInfo(
        body_radius_km=body_radius_km,
        periapsis_altitude_km=peri_alt,
        apoapsis_altitude_km=apo_alt,
        warning=warning,
    )


@dataclass(frozen=True)
class PropellantInfo:
    """Propellant mass fraction from the Tsiolkovsky rocket equation,
    requiring the engine's specific impulse (Isp, in seconds) -- an input
    this orbital-mechanics component never has on its own. `delta_v_km_s`
    is echoed back for traceability; `mass_fraction` = 1 - exp(-delta_v /
    (Isp*g0)) is the fraction of the *initial* (wet) mass that must be
    propellant to achieve delta_v_km_s with a single engine of this Isp
    (two separate burns, at departure and arrival, each with their own
    delta-v, in general need a fresh mass_fraction applied to the mass
    remaining *after* the first burn -- this function reports the
    fraction for one delta-v in isolation; combining two burns' fractions
    is the caller's responsibility and is NOT done here, to avoid
    silently assuming a specific burn order or staging strategy)."""

    delta_v_km_s: float
    isp_s: float
    mass_fraction: float


def propellant_info(delta_v_km_s: float, isp_s: float) -> PropellantInfo:
    if isp_s <= 0:
        raise ValueError(f"isp_s must be positive, got {isp_s!r}")
    delta_v_m_s = delta_v_km_s * 1000.0
    mass_fraction = 1.0 - math.exp(-delta_v_m_s / (isp_s * STANDARD_GRAVITY_M_S2))
    return PropellantInfo(delta_v_km_s=delta_v_km_s, isp_s=isp_s, mass_fraction=mass_fraction)


@dataclass(frozen=True)
class LambertCrossCheck:
    """Result of cross-checking a solved transfer against an independent
    Lambert-problem solve (`lambert.solve_lambert`) fed this transfer's own
    time of flight -- a genuinely separate classical method (see
    `lambert.py`'s module docstring), not a re-reading of this package's
    own S1/S2 algebra (that is `geometry.GeometryCheck`'s job)."""

    p_km_residual: float
    eccentricity_residual: float
    delta_v1_km_s_residual: float
    delta_v2_km_s_residual: float
    lambert_iterations: int

    @property
    def max_abs_residual(self) -> float:
        return max(
            abs(self.p_km_residual),
            abs(self.eccentricity_residual),
            abs(self.delta_v1_km_s_residual),
            abs(self.delta_v2_km_s_residual),
        )

    def ok(self, tol: float = 1e-6) -> bool:
        return self.max_abs_residual < tol


def lambert_cross_check(
    r1_km: float,
    r2_km: float,
    delta_theta_rad: float,
    tof_s: float,
    mu_km3_s2: float,
    p_star_km: float,
    e: float,
    delta_v1_km_s: float,
    delta_v2_km_s: float,
) -> LambertCrossCheck:
    """Re-solve the transfer from scratch via an independent Lambert
    solver (departure/arrival position + the already-computed time of
    flight as the three inputs Lambert's problem actually takes), and
    compare its orbit (p, e) and its implied delta-v magnitudes against
    the ones this package's own cost-minimization solve produced.

    `delta_theta_rad` is the actual swept angle (0, 2*pi), consistent with
    `api.compute_transfer`'s own folded `delta_theta_deg_folded`
    convention -- NOT the signed `phi` model parameter.
    """
    r1_vec = (r1_km, 0.0)
    r2_vec = (r2_km * math.cos(delta_theta_rad), r2_km * math.sin(delta_theta_rad))
    sol = solve_lambert(r1_vec, r2_vec, tof_s, mu_km3_s2, prograde=True)

    vc1 = (mu_km3_s2 / r1_km) ** 0.5
    vc2 = (mu_km3_s2 / r2_km) ** 0.5
    vc1_vec = (0.0, vc1)
    vc2_vec = (-math.sin(delta_theta_rad) * vc2, math.cos(delta_theta_rad) * vc2)

    lambert_dv1 = (
        (sol.v1_km_s[0] - vc1_vec[0]) ** 2 + (sol.v1_km_s[1] - vc1_vec[1]) ** 2
    ) ** 0.5
    lambert_dv2 = (
        (vc2_vec[0] - sol.v2_km_s[0]) ** 2 + (vc2_vec[1] - sol.v2_km_s[1]) ** 2
    ) ** 0.5

    return LambertCrossCheck(
        p_km_residual=sol.p_km - p_star_km,
        eccentricity_residual=sol.e - e,
        delta_v1_km_s_residual=lambert_dv1 - delta_v1_km_s,
        delta_v2_km_s_residual=lambert_dv2 - delta_v2_km_s,
        lambert_iterations=sol.iterations,
    )
