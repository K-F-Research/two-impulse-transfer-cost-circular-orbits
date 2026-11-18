"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: Theorem 1 (global monotonicity of V_R(phi) in
phi -- establishes that Hohmann, phi=0, is the unique global minimum
of the cost, so a budget below it is infeasible, and that the
inversion solved below has at most one root) and Theorem 3 (strict
convexity of V_R in z := cot^2(phi/2), for every R>1 -- see module
docstring below for exactly how this is used).

WHAT THIS MODULE IS, IN ONE SENTENCE
-------------------------------------
`api.compute_transfer` and `core.solve` answer "given a transfer angle,
what is its cost?". This module answers the inverse question: "given a
delta-v budget, what is the LARGEST angular deviation from Hohmann
(equivalently, the smallest `delta_theta`) that budget can afford?".

WHY THIS IS SAFE: THEOREM 1 (EXISTENCE AND UNIQUENESS OF THE INVERSE)
-----------------------------------------------------------------------
Theorem 1 proves V_R(phi) is strictly increasing on the open interval
phi in (0, pi), from V_R(0) (the Hohmann minimum) up to +infinity as
phi -> pi. So for any budget strictly greater than the Hohmann cost,
there is exactly one phi_max in (0, pi) with V_R(phi_max) = budget, and
inverting V_R is a well-posed, safe root-find -- never a search for a
root that might not exist or might not be unique.

Two degenerate cases are handled directly, without any root-finding:
budget < Hohmann cost is infeasible (DomainError: Theorem 1 makes
Hohmann the global minimum, so no phi can cost less); budget very close
to the Hohmann cost (within `hohmann_boundary_tol_km_s`) is resolved
EXACTLY as phi_max=0, rather than asked of a Newton iteration whose
target would be z=cot^2(phi/2) -> +infinity as phi -> 0 -- a numerically
hostile target for no reason, when the exact closed-form answer is
already known from `hohmann.py`.

WHY THEOREM 3 MAKES THE SOLVER BETTER, NOT JUST "ALSO POSSIBLE"
-------------------------------------------------------------------
Theorem 1 alone already justifies a safe bisection on phi -- this would
work, but only converges linearly. Theorem 3 proves V_R is additionally
*strictly convex* in z := cot^2(phi/2), for every R>1: so g(z) :=
V_R(z) - budget has a strictly monotonic derivative dV/dz, which is
exactly the structural property a safeguarded Newton iteration needs to
be well-behaved. Convexity alone does not, by itself and without further
conditions on the starting point and step strategy, guarantee that a
bare Newton iteration converges globally and monotonically from any
starting point -- so this module does not claim that. What it implements
is:

    safeguarded Newton with monotone bracketing, using the strict
    convexity established by Theorem 3.

Concretely: `_invert_phi_newton` maintains a bracket [phi_lo, phi_hi]
in phi (not z) with g_phi(phi_lo) < 0 < g_phi(phi_hi), guaranteed to
exist and to be unique by Theorem 1. Each iteration proposes a Newton
step computed in z-space (`dV/dz` obtained exactly via the envelope
theorem and the chain rule -- see `_dV_dz`, no finite difference
anywhere), converts the proposed z back to a candidate phi, and ACCEPTS
it only if it strictly tightens the bracket; otherwise it falls back to
a bisection midpoint in phi. This is what makes the method safe from an
arbitrary starting point, not the convexity by itself. When the Newton
step is accepted (the common case once the bracket is reasonably
tight), convergence is quadratic; the bisection fallback guarantees the
bracket still shrinks geometrically even on the iterations where it
is not.

Defensive, run-time checks (never a silent wrong answer): at every
iteration this function checks that `dV/dz` is finite and nonzero
before using it for a Newton step (if not, it falls back to bisection
for that step rather than dividing by a degenerate derivative), and
checks after each new evaluation that the bracket invariant
g_phi(phi_lo) < 0 < g_phi(phi_hi) still holds (raising ConvergenceError,
not silently continuing, if it does not -- this should never happen for
a correctly posed problem and would indicate an implementation defect).
"""
import math
from dataclasses import dataclass

from . import core, hohmann
from .errors import ConvergenceError, DomainError
from .units import EARTH_MU_KM3_S2, circular_speed_km_s


@dataclass(frozen=True)
class BudgetInversionResult:
    """The largest angular deviation from Hohmann that a given delta-v
    budget affords, in real physical units, together with diagnostics."""

    # Inputs, echoed back for traceability.
    r1_km: float
    r2_km: float
    mu_km3_s2: float
    delta_v_budget_km_s: float

    # Derived model parameter (dimensionless).
    R: float

    # The answer.
    phi_max_rad: float
    delta_theta_min_deg: float
    p_star_km: float

    # Reference and diagnostics.
    hohmann_delta_v_total_km_s: float
    method: str  # "exact_hohmann_boundary" or "safeguarded_newton_in_z"
    iterations: int
    newton_steps_accepted: int
    bisection_steps_taken: int
    residual_km_s: float  # V_R(phi_max) - budget, should be ~0


def _dV_dz(phi: float, R: float, p_star: float) -> float:
    """dV_R/dz at the given phi (z := cot^2(phi/2)), computed exactly: the
    envelope theorem gives dV/dphi = C_phi(p*,phi,R) (core.cost_partial_phi,
    since C_p(p*,phi,R)=0 by construction of p*), and the chain rule gives
    dz/dphi = -cos(phi/2) / sin(phi/2)**3 in closed form. No finite
    difference is used anywhere in this module."""
    dV_dphi = core.cost_partial_phi(p_star, phi, R)
    half = phi / 2.0
    s, c = math.sin(half), math.cos(half)
    dz_dphi = -c / (s ** 3)
    return dV_dphi / dz_dphi


def _z_of_phi(phi: float) -> float:
    return (math.cos(phi / 2.0) / math.sin(phi / 2.0)) ** 2


def _phi_of_z(z: float) -> float:
    """Inverse of z = cot^2(phi/2) for phi in (0, pi): cot(phi/2) = sqrt(z)
    (positive root, since phi/2 in (0, pi/2) makes cot(phi/2) > 0 there)."""
    return 2.0 * math.atan2(1.0, z ** 0.5)


def _invert_phi_newton(
    R: float,
    budget_norm: float,
    phi_lo: float,
    phi_hi: float,
    xtol_rad: float = 1e-13,
    ftol_norm: float = 1e-13,
    max_iterations: int = 200,
) -> "tuple[float, float, int, int, int]":
    """Safeguarded Newton (in z) with monotone bracketing (in phi), as
    described in this module's docstring. Returns (phi_root, p_star,
    iterations, newton_steps_accepted, bisection_steps_taken)."""

    def g_phi(phi: float) -> "tuple[float, float]":
        p_star = core.find_p_star(phi, R)
        V = core.cost(p_star, phi, R)
        return V - budget_norm, p_star

    g_lo, p_lo = g_phi(phi_lo)
    g_hi, p_hi = g_phi(phi_hi)
    if not (g_lo < 0.0 < g_hi):
        raise ConvergenceError(
            f"initial bracket invariant failed: g(phi_lo={phi_lo!r})="
            f"{g_lo!r}, g(phi_hi={phi_hi!r})={g_hi!r}; expected g_lo<0<g_hi "
            "(Theorem 1's monotonicity). This indicates the caller's "
            "bracket was not actually a valid bracket for this R/budget."
        )

    newton_accepted = 0
    bisections = 0
    phi = 0.5 * (phi_lo + phi_hi)
    p_star = None
    for iteration in range(1, max_iterations + 1):
        g_mid, p_star = g_phi(phi)
        if abs(g_mid) < ftol_norm or (phi_hi - phi_lo) < xtol_rad:
            return phi, p_star, iteration, newton_accepted, bisections

        if g_mid < 0.0:
            phi_lo, g_lo = phi, g_mid
        else:
            phi_hi, g_hi = phi, g_mid
        if not (g_lo < 0.0 < g_hi):
            raise ConvergenceError(
                f"bracket invariant violated after update at iteration "
                f"{iteration} (phi={phi!r}, g={g_mid!r}); this should "
                "never happen for a correctly posed problem and indicates "
                "an implementation defect, not a genuine domain edge case."
            )

        # Propose a Newton step in z-space, using the exact envelope-
        # theorem derivative; fall back to bisection if the derivative is
        # degenerate or the proposed point does not strictly tighten the
        # bracket (see module docstring: this is what makes the method
        # safe, not the convexity alone).
        proposed_phi = None
        try:
            slope = _dV_dz(phi, R, p_star)
            if math.isfinite(slope) and slope != 0.0:
                z_mid = _z_of_phi(phi)
                z_newton = z_mid - g_mid / slope
                if math.isfinite(z_newton) and z_newton > 0.0:
                    candidate = _phi_of_z(z_newton)
                    if math.isfinite(candidate) and phi_lo < candidate < phi_hi:
                        proposed_phi = candidate
                        newton_accepted += 1
        except (ZeroDivisionError, ValueError, OverflowError):
            proposed_phi = None

        if proposed_phi is None:
            proposed_phi = 0.5 * (phi_lo + phi_hi)
            bisections += 1
        phi = proposed_phi

    raise ConvergenceError(
        f"safeguarded Newton did not converge within {max_iterations} "
        f"iterations for R={R!r}, budget_norm={budget_norm!r}; final "
        f"bracket width was {phi_hi - phi_lo!r} rad."
    )


def invert_phi_for_budget(
    r1_km: float,
    r2_km: float,
    delta_v_budget_km_s: float,
    mu_km3_s2: float = EARTH_MU_KM3_S2,
    hohmann_boundary_tol_km_s: float = 1e-9,
) -> BudgetInversionResult:
    """Find the largest angular deviation from Hohmann that a given
    delta-v budget affords: the unique phi_max in (0, pi) with
    V_R(phi_max) = delta_v_budget_km_s (converted to normalized units),
    or the exact phi_max=0 boundary case when the budget is (numerically)
    exactly the Hohmann cost.

    Parameters
    ----------
    r1_km, r2_km : departure and arrival circular-orbit radii, in km.
        Same convention as `api.compute_transfer` (requires r2_km >
        r1_km).
    delta_v_budget_km_s : the total delta-v budget, in km/s.
    mu_km3_s2 : gravitational parameter of the central body, km^3/s^2.
        Defaults to Earth.
    hohmann_boundary_tol_km_s : if the budget is within this tolerance of
        the exact Hohmann cost, phi_max is reported as exactly 0 (the
        closed-form `hohmann.py` boundary case) rather than solved
        numerically -- see module docstring for why this avoids asking a
        Newton iteration in z to converge toward z -> +infinity.

    Raises
    ------
    DomainError
        If the radii are invalid, or the budget is strictly below the
        Hohmann cost (infeasible: Theorem 1 makes Hohmann the global
        minimum of the cost over the whole domain).
    ConvergenceError
        If the safeguarded Newton/bisection solve unexpectedly fails to
        converge despite the input being inside the domain Theorem 1
        guarantees a unique root on (see `_invert_phi_newton`).
    """
    if not (r1_km > 0 and r2_km > 0):
        raise DomainError(f"radii must be positive: r1_km={r1_km!r}, r2_km={r2_km!r}")
    if not r2_km > r1_km:
        raise DomainError(
            f"this function requires r2_km > r1_km (R:=r2/r1>1); got "
            f"r1_km={r1_km!r}, r2_km={r2_km!r}."
        )
    if not delta_v_budget_km_s > 0:
        raise DomainError(f"delta_v_budget_km_s must be positive, got {delta_v_budget_km_s!r}")

    R = r2_km / r1_km
    vc1 = circular_speed_km_s(r1_km, mu_km3_s2)
    hohmann_ref = hohmann.solve_hohmann(R)
    hohmann_total_km_s = hohmann_ref.delta_v_total * vc1

    if delta_v_budget_km_s < hohmann_total_km_s - hohmann_boundary_tol_km_s:
        raise DomainError(
            f"delta_v_budget_km_s={delta_v_budget_km_s!r} is below the "
            f"Hohmann cost ({hohmann_total_km_s!r} km/s) for R={R!r}: "
            "infeasible, since Theorem 1 establishes the Hohmann transfer "
            "(phi=0) as the unique GLOBAL minimum of the cost over the "
            "entire domain -- no angular deviation can cost less."
        )

    if abs(delta_v_budget_km_s - hohmann_total_km_s) <= hohmann_boundary_tol_km_s:
        return BudgetInversionResult(
            r1_km=r1_km, r2_km=r2_km, mu_km3_s2=mu_km3_s2,
            delta_v_budget_km_s=delta_v_budget_km_s, R=R,
            phi_max_rad=0.0, delta_theta_min_deg=180.0,
            p_star_km=hohmann_ref.p * r1_km,
            hohmann_delta_v_total_km_s=hohmann_total_km_s,
            method="exact_hohmann_boundary",
            iterations=0, newton_steps_accepted=0, bisection_steps_taken=0,
            residual_km_s=delta_v_budget_km_s - hohmann_total_km_s,
        )

    budget_norm = delta_v_budget_km_s / vc1

    # Bracket phi away from both domain boundaries: phi_lo close to 0 (but
    # not equal to it -- the Hohmann boundary is already handled above),
    # phi_hi expanded toward pi (geometrically, mirroring core.py's own
    # bracket-expansion pattern for find_p_star) until V_R(phi_hi) exceeds
    # the budget; Theorem 1 guarantees this expansion terminates.
    phi_lo = 1e-6
    phi_hi = math.pi - 1e-6
    max_expansions = 300
    expansions = 0
    while True:
        p_hi = core.find_p_star(phi_hi, R)
        if core.cost(p_hi, phi_hi, R) > budget_norm:
            break
        expansions += 1
        if expansions >= max_expansions or (math.pi - phi_hi) < 1e-15:
            raise ConvergenceError(
                f"could not bracket phi_max for R={R!r}, "
                f"budget_km_s={delta_v_budget_km_s!r} after {expansions} "
                "expansions toward phi=pi; this input may be on or past a "
                "genuine domain boundary."
            )
        phi_hi = math.pi - (math.pi - phi_hi) / 2.0

    phi_root, p_star, iterations, newton_accepted, bisections = _invert_phi_newton(
        R, budget_norm, phi_lo, phi_hi
    )
    V_root = core.cost(p_star, phi_root, R)

    return BudgetInversionResult(
        r1_km=r1_km, r2_km=r2_km, mu_km3_s2=mu_km3_s2,
        delta_v_budget_km_s=delta_v_budget_km_s, R=R,
        phi_max_rad=phi_root,
        delta_theta_min_deg=180.0 - math.degrees(phi_root),
        p_star_km=p_star * r1_km,
        hohmann_delta_v_total_km_s=hohmann_total_km_s,
        method="safeguarded_newton_in_z",
        iterations=iterations,
        newton_steps_accepted=newton_accepted,
        bisection_steps_taken=bisections,
        residual_km_s=(V_root - budget_norm) * vc1,
    )
