"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: none (independent numerical propagation of the
two-body equations of motion, used only as a cross-check method -- see
the module description below).

Independent validation of a solved transfer by direct numerical
propagation of the two-body equations of motion (classical fixed-step
RK4), as a THIRD, genuinely distinct cross-check alongside
`geometry.independent_check` and `derived.lambert_cross_check`:

- `geometry.independent_check` re-reads the same solved (p*, e) through
  different textbook formulas (conic equation, vis-viva), but never
  integrates anything over time.
- `derived.lambert_cross_check` (`lambert.py`) is a from-scratch
  *boundary-value* solve: given r1, r2 and the time of flight, it
  iteratively searches (root-finding on a Stumpff-function series) for
  the orbit connecting them. It shares no algebra with this package's own
  solve, but it is still an algebraic/semi-analytic method, and its own
  input is this package's own `time_of_flight_s`.
- This module is an *initial-value* numerical integrator: it takes the
  transfer orbit's departure state vector (position and velocity, in
  Cartesian coordinates, obtained directly from p*, e, nu1 via the
  standard closed-form conic-to-Cartesian formulas below -- no
  root-finding at all), and simply integrates Newton's law of
  gravitation (`d^2 r/dt^2 = -mu * r / |r|^3`) forward for the already-
  computed `time_of_flight_s`, then checks whether the numerically
  propagated state lands on the already-known arrival position and
  velocity. It shares no formulas with the Lambert solver's z/Stumpff
  machinery, and unlike `geometry.independent_check`, it genuinely
  exercises `time_of_flight_s` (Kepler's equation) against the actual
  dynamics, not just the conic shape of the orbit.

A convergence study (Richardson-style step doubling) accompanies every
call: a single RK4 run at one step count proves nothing about accuracy by
itself (RK4 converges at 4th order, but only in the step-size regime
where its own truncation error dominates numerical noise) -- so this
module always runs at both `n_steps` and `2*n_steps` and reports the
observed order of convergence, raising `ConvergenceError` if it is not at
least roughly consistent with genuine 4th-order behavior (allowing for
the point, near full convergence, where floating-point noise rather than
truncation error dominates both runs).

Scope, stated explicitly (an independent review of this component's
earlier cross-checks asked that this be made unambiguous): this module
does NOT independently re-derive p* or e -- those are taken as already
solved (exactly like `geometry.independent_check` takes them as given).
What it DOES independently exercise is `time_of_flight_s` (Kepler's
equation) and the orbit's actual dynamical consistency, by checking that
literally integrating Newton's law from the known departure state for
the known duration lands on the known arrival state -- something neither
`geometry_check` (no time integration at all) nor `lambert_check` (an
algebraic boundary-value method, with its own, different, root-finding
machinery) actually tests. A low energy/angular-momentum drift along the
trajectory is reported as a sanity check on the integrator's own
numerical behavior, but is NOT by itself sufficient evidence that the
propagated trajectory is correct (a systematically biased integrator
could in principle conserve both while still drifting from the true
path) -- `ok()` / the tolerance checks inside `rk4_validate` compare the
propagated FINAL STATE against the independently-computed analytic
arrival state, which is the actual pass/fail criterion; the drift
figures are diagnostic, not sufficient.
"""
import math
from dataclasses import dataclass
from typing import Optional, Tuple

from .errors import ConvergenceError

Vec2 = Tuple[float, float]

# Below this absolute step-to-step residual (km or km/s), floating-point
# noise -- not RK4 truncation error -- dominates, and the 4th-order
# convergence ratio becomes meaningless (division of two numbers both at
# the noise floor). Below this floor, both runs already agree with the
# analytic arrival state to a precision no verification method run in
# ordinary double precision could usefully improve on.
_NOISE_FLOOR_KM = 1e-6
_NOISE_FLOOR_KM_S = 1e-9


def _state_vector_from_elements(
    p_km: float, e: float, nu1_rad: float, theta_rad: float, mu_km3_s2: float
) -> Tuple[Vec2, Vec2]:
    """Position and velocity, in the package's own (r1_km, 0) / (r2_km *
    cos(delta_theta), r2_km * sin(delta_theta)) planar frame (see
    `derived.lambert_cross_check`), at the point of the transfer orbit
    whose true anomaly is `nu1_rad + theta_rad` -- `theta_rad` being the
    angle swept *from departure* (0 at departure, `delta_theta_rad` --
    the actual folded sweep angle -- at arrival), NOT the true anomaly
    itself (nu1 is in general nonzero; it is an output of the solve, not
    a frame choice).

    Obtained directly from the conic equation r = p/(1+e*cos(nu)) and the
    angular-momentum/vis-viva relations h = sqrt(mu*p),
    v_radial = sqrt(mu/p)*e*sin(nu) -- no root-finding, no reference to
    `core.py`'s S1/S2 algebra or to `lambert.py`'s Stumpff functions.
    """
    nu = nu1_rad + theta_rad
    r_km = p_km / (1.0 + e * math.cos(nu))
    h_km2_s = (mu_km3_s2 * p_km) ** 0.5
    v_tangential_km_s = h_km2_s / r_km
    v_radial_km_s = (mu_km3_s2 / p_km) ** 0.5 * e * math.sin(nu)

    c, s = math.cos(theta_rad), math.sin(theta_rad)
    x_km = r_km * c
    y_km = r_km * s
    vx_km_s = v_radial_km_s * c - v_tangential_km_s * s
    vy_km_s = v_radial_km_s * s + v_tangential_km_s * c
    return (x_km, y_km), (vx_km_s, vy_km_s)


def _accel(pos: Vec2, mu_km3_s2: float) -> Vec2:
    x, y = pos
    r3 = (x * x + y * y) ** 1.5
    return (-mu_km3_s2 * x / r3, -mu_km3_s2 * y / r3)


def _rk4_step(
    pos: Vec2, vel: Vec2, mu_km3_s2: float, h_s: float
) -> Tuple[Vec2, Vec2]:
    """One classical 4th-order Runge-Kutta step on the first-order system
    (dpos/dt = vel, dvel/dt = accel(pos)); `h_s` is the step size in
    seconds."""

    def deriv(p: Vec2, v: Vec2) -> Tuple[Vec2, Vec2]:
        return v, _accel(p, mu_km3_s2)

    k1p, k1v = deriv(pos, vel)
    k2p, k2v = deriv(
        (pos[0] + 0.5 * h_s * k1p[0], pos[1] + 0.5 * h_s * k1p[1]),
        (vel[0] + 0.5 * h_s * k1v[0], vel[1] + 0.5 * h_s * k1v[1]),
    )
    k3p, k3v = deriv(
        (pos[0] + 0.5 * h_s * k2p[0], pos[1] + 0.5 * h_s * k2p[1]),
        (vel[0] + 0.5 * h_s * k2v[0], vel[1] + 0.5 * h_s * k2v[1]),
    )
    k4p, k4v = deriv(
        (pos[0] + h_s * k3p[0], pos[1] + h_s * k3p[1]),
        (vel[0] + h_s * k3v[0], vel[1] + h_s * k3v[1]),
    )
    new_pos = (
        pos[0] + (h_s / 6.0) * (k1p[0] + 2 * k2p[0] + 2 * k3p[0] + k4p[0]),
        pos[1] + (h_s / 6.0) * (k1p[1] + 2 * k2p[1] + 2 * k3p[1] + k4p[1]),
    )
    new_vel = (
        vel[0] + (h_s / 6.0) * (k1v[0] + 2 * k2v[0] + 2 * k3v[0] + k4v[0]),
        vel[1] + (h_s / 6.0) * (k1v[1] + 2 * k2v[1] + 2 * k3v[1] + k4v[1]),
    )
    return new_pos, new_vel


def _propagate(
    pos0: Vec2, vel0: Vec2, mu_km3_s2: float, tof_s: float, n_steps: int
) -> Tuple[Vec2, Vec2, float, float]:
    """Propagate for `tof_s` seconds in `n_steps` equal RK4 steps; also
    returns the max relative drift of specific energy and of the (scalar,
    z-component of the) specific angular momentum observed along the
    trajectory, as an independent sanity check ON THE INTEGRATOR ITSELF
    (these are exactly conserved by the true continuous dynamics; a
    fixed-step RK4 integrator conserves them only approximately, with the
    drift shrinking as the step size shrinks)."""
    h_s = tof_s / n_steps
    pos, vel = pos0, vel0
    r0 = (pos0[0] ** 2 + pos0[1] ** 2) ** 0.5
    energy0 = 0.5 * (vel0[0] ** 2 + vel0[1] ** 2) - mu_km3_s2 / r0
    ang_mom0 = pos0[0] * vel0[1] - pos0[1] * vel0[0]
    max_energy_drift = 0.0
    max_ang_mom_drift = 0.0
    for _ in range(n_steps):
        pos, vel = _rk4_step(pos, vel, mu_km3_s2, h_s)
        r = (pos[0] ** 2 + pos[1] ** 2) ** 0.5
        if r <= 0 or not math.isfinite(r):
            raise ConvergenceError(
                "RK4 propagation diverged (non-finite or non-positive "
                f"radius) after a step with h_s={h_s!r} s, n_steps={n_steps!r}"
            )
        energy = 0.5 * (vel[0] ** 2 + vel[1] ** 2) - mu_km3_s2 / r
        ang_mom = pos[0] * vel[1] - pos[1] * vel[0]
        if energy0 != 0:
            max_energy_drift = max(max_energy_drift, abs((energy - energy0) / energy0))
        if ang_mom0 != 0:
            max_ang_mom_drift = max(max_ang_mom_drift, abs((ang_mom - ang_mom0) / ang_mom0))
    return pos, vel, max_energy_drift, max_ang_mom_drift


@dataclass(frozen=True)
class RK4Validation:
    """Result of propagating the transfer orbit's own departure state
    vector forward by numerical integration (classical two-body dynamics
    only -- no reference to this package's own cost/geometry/Lambert
    algebra) and comparing the result to the already-known arrival state.

    `position_residual_km`/`velocity_residual_km_s` are the final-run
    (finest step size) discrepancies against the analytic arrival state.
    `observed_convergence_order` is the empirical order estimated from
    halving the step size (should be close to 4 for classical RK4, in the
    regime where truncation error dominates floating-point noise --
    `None` when both runs were already at the noise floor, which is
    itself a form of convergence, just not one this ratio can measure).
    `max_energy_drift`/`max_ang_mom_drift` are conservation-law checks on
    the integrator itself, independent of the arrival-state comparison.
    """

    position_residual_km: float
    velocity_residual_km_s: float
    n_steps: int
    observed_convergence_order: Optional[float]
    max_energy_drift: float
    max_ang_mom_drift: float

    def ok(self, position_tol_km: float = 1e-3, velocity_tol_km_s: float = 1e-6) -> bool:
        return (
            self.position_residual_km < position_tol_km
            and self.velocity_residual_km_s < velocity_tol_km_s
        )


def rk4_validate(
    p_star_km: float,
    e: float,
    nu1_rad: float,
    delta_theta_rad: float,
    tof_s: float,
    mu_km3_s2: float,
    r2_km: float,
    n_steps: int = 4000,
    position_tol_km: float = 1e-3,
    velocity_tol_km_s: float = 1e-6,
    max_doublings: int = 8,
) -> RK4Validation:
    """Independently validate a solved transfer by numerically
    propagating its own departure state for `tof_s` seconds and comparing
    to its own (analytically computed) arrival state.

    `n_steps` is only a STARTING step count, not a fixed one: a transfer
    whose swept arc passes close to periapsis of a highly eccentric orbit
    (confirmed empirically -- e.g. delta_theta=350 deg, R=2.0 around
    Earth has e=0.982 and a periapsis passage squarely inside the swept
    170-350 deg arc) needs many more, much finer steps than an ordinary
    well-conditioned case to resolve the sharp curvature there; a single
    fixed step count that works for one geometry can be grossly
    insufficient for another. This function therefore DOUBLES `n_steps`
    (up to `max_doublings` times) until either (a) the finer run's
    residual against the analytic arrival state is within tolerance, or
    (b) the residual is already at the floating-point noise floor -- at
    each doubling checking that the observed convergence order is
    consistent with RK4's genuine 4th-order behavior (not just that the
    residual happens to be small), exactly the same "verify the expected
    sign/rate, do not just hope for the best" discipline already used by
    `core.find_p_star`'s and `lambert.solve_lambert`'s own adaptive
    searches. Raises `ConvergenceError` (never a silently inaccurate
    result) if tolerance is not reached within `max_doublings`, or if any
    doubling fails to show the expected convergence order.

    `delta_theta_rad` is the actual folded sweep angle in (0, 2*pi),
    exactly as `derived.lambert_cross_check` expects it -- NOT the signed
    `phi` model parameter.
    """
    if not (0.0 <= e < 1.0):
        raise ConvergenceError(
            f"rk4_validate requires an ellipse, 0<=e<1 (needed for a "
            f"finite, bounded departure state to propagate); got e={e!r}"
        )
    if tof_s <= 0:
        raise ConvergenceError(f"time of flight must be positive, got tof_s={tof_s!r}")

    pos0, vel0 = _state_vector_from_elements(p_star_km, e, nu1_rad, 0.0, mu_km3_s2)
    pos_arrival, vel_arrival = _state_vector_from_elements(
        p_star_km, e, nu1_rad, delta_theta_rad, mu_km3_s2
    )

    def dist(u: Vec2, v: Vec2) -> float:
        return ((u[0] - v[0]) ** 2 + (u[1] - v[1]) ** 2) ** 0.5

    pos_a, vel_a, energy_drift_a, ang_mom_drift_a = _propagate(
        pos0, vel0, mu_km3_s2, tof_s, n_steps
    )
    pos_residual_a = dist(pos_a, pos_arrival)

    for attempt in range(max_doublings):
        finer_steps = 2 * n_steps
        pos_b, vel_b, energy_drift_b, ang_mom_drift_b = _propagate(
            pos0, vel0, mu_km3_s2, tof_s, finer_steps
        )
        pos_residual_b = dist(pos_b, pos_arrival)
        vel_residual_b = dist(vel_b, vel_arrival)

        observed_order = None
        if pos_residual_a > _NOISE_FLOOR_KM and pos_residual_b > _NOISE_FLOOR_KM:
            # RK4's global error scales as h^4; halving h should shrink
            # the error by a factor of about 2^4 = 16 -- log2 of the
            # observed ratio estimates the empirical order directly.
            ratio = pos_residual_a / pos_residual_b
            if ratio <= 1.0:
                raise ConvergenceError(
                    "doubling the RK4 step count did not reduce the "
                    f"residual against the analytic arrival state "
                    f"({pos_residual_a!r} km at n_steps={n_steps!r} vs "
                    f"{pos_residual_b!r} km at n_steps={finer_steps!r} -- "
                    f"ratio={ratio!r}), inconsistent with a correctly "
                    "implemented convergent integrator."
                )
            observed_order = math.log2(ratio)
            if observed_order < 2.5:
                raise ConvergenceError(
                    "RK4 step-doubling did not show the expected ~4th-order "
                    f"convergence (observed order={observed_order:.2f} from "
                    f"residuals {pos_residual_a!r} km at n_steps={n_steps!r} "
                    f"vs {pos_residual_b!r} km at n_steps={finer_steps!r}); "
                    "the integrator or the analytic arrival state used for "
                    "comparison may be wrong -- this is reported rather "
                    "than trusted silently."
                )

        if (
            pos_residual_b < position_tol_km
            and vel_residual_b < velocity_tol_km_s
        ) or (
            pos_residual_b <= _NOISE_FLOOR_KM
            and vel_residual_b <= _NOISE_FLOOR_KM_S
        ):
            return RK4Validation(
                position_residual_km=pos_residual_b,
                velocity_residual_km_s=vel_residual_b,
                n_steps=finer_steps,
                observed_convergence_order=observed_order,
                max_energy_drift=max(energy_drift_a, energy_drift_b),
                max_ang_mom_drift=max(ang_mom_drift_a, ang_mom_drift_b),
            )

        # Not yet within tolerance (and convergence order, where
        # measurable, was consistent with genuine 4th-order behavior) --
        # this geometry needs finer resolution still; double again.
        n_steps = finer_steps
        pos_a, vel_a, energy_drift_a, ang_mom_drift_a = pos_b, vel_b, energy_drift_b, ang_mom_drift_b
        pos_residual_a = pos_residual_b

    raise ConvergenceError(
        f"RK4 validation did not reach tolerance (position_tol_km="
        f"{position_tol_km!r}, velocity_tol_km_s={velocity_tol_km_s!r}) "
        f"within {max_doublings} step-count doublings (reached n_steps="
        f"{n_steps!r}, position_residual_km={pos_residual_a!r}); this "
        "geometry (e.g. a very high eccentricity with a periapsis "
        "passage inside the swept arc) may need a larger starting "
        "n_steps or a higher max_doublings than this call used."
    )
