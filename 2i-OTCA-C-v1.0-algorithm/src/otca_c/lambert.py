"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: none (independent classical Lambert-problem
solver, used only as a cross-check method -- see the module
description below).

An independent, textbook Lambert-problem solver (universal-variable
formulation), used as a *genuinely separate* cross-check of a solved
transfer -- not a re-reading of this package's own S1/S2 algebra (that is
`geometry.py`'s job), but a different classical method answering a
differently-posed question: "given a fixed departure point, arrival point,
and time of flight, what orbit connects them?" (Lambert's problem), versus
this package's own question, "given a fixed transfer angle, what is the
minimum-cost orbit?" (an optimization over p).

When fed the time of flight this package's own solution implies
(`derived.time_of_flight_s`), a correctly solved transfer must make this
independently-implemented Lambert solver recover the identical transfer
orbit (same p, e) and the identical two impulsive delta-v magnitudes. This
is the comparison against "classical/reference methods" requested for this
component, in addition to (not instead of) the internal consistency check
already in `geometry.py`.

Reference: R. H. Battin's and (more directly, the algorithm followed here)
D. Vallado's "Fundamentals of Astrodynamics and Applications", 4th ed.,
Algorithm 58 ("Lambertuniv"), using the Stumpff functions C(z)/S(z) -- the
standard, numerically robust universal-variable formulation that covers
ellipse (z>0), parabola (z=0) and hyperbola (z<0) with one unified
iteration, avoiding the coordinate/quadrant special-casing of older
eccentric/hyperbolic-anomaly-based Lambert methods.

Root-finding strategy: the time-of-flight equation t(z) is solved with an
adaptive bracket search followed by `scipy.optimize.brentq` (the same
proven pattern as `core.find_p_star`'s p* search), not a hand-rolled
Newton iteration. An earlier version of this module used Newton with a
finite-difference derivative and clamped z to a fixed domain; an
independent review found that this clamp could bias the finite-difference
derivative and stall convergence for transfer angles near the 0/360 deg
boundary (confirmed: z jumping straight to the clamp and needing ~9
further iterations to recover, for a total of 47-48 iterations instead of
the few an encadrement method needs -- see
`tests/test_lambert.py::test_lambert_near_360_boundary_converges_without_excessive_iterations`,
kept as a permanent regression guard). The bracketing approach never
clamps or silently projects an invalid evaluation onto the domain
boundary: an evaluation outside the valid domain is simply treated as
undefined, and failure to find a valid bracket raises `ConvergenceError`
explicitly rather than degrading silently.

This implementation only targets this package's own domain: a single
(zero-revolution), prograde, non-degenerate transfer between two distinct
circular orbits (0 < transfer angle < 360 deg, excluding exactly 180 deg,
which is a degeneracy of this solver's own A-formula -- not of the
underlying geometry -- and is routed around entirely via the exact
Hohmann closed form in `hohmann.py`/`api.py` instead). It is not a
general-purpose Lambert solver for retrograde transfers or multiple
revolutions.

Domain of z actually enforced here: `-4*pi**2 < z < 4*pi**2`. The upper
bound is a true mathematical boundary (C(z) has its first zero at
z=4*pi**2, the single-rev/multi-rev transition). The LOWER bound is not:
the universal-variable formulation's hyperbolic branch (z<0) has no
comparable natural stopping point, and z=-4*pi**2 is simply a practical
cutoff adequate for this package's own elliptical transfers, not a claim
that every hyperbolic Lambert solution lies above it. This module is a
cross-check for this package's own (always-elliptical, see `derived.py`'s
own e<1 assumption) transfers, not a general-purpose Lambert library, and
this asymmetry in how the two bounds are justified should not be read as
an oversight.
"""
import math
from dataclasses import dataclass

from scipy.optimize import brentq

from .errors import ConvergenceError


def _stumpff_c(z: float) -> float:
    """C(z): 1/2 - z/24 + z^2/720 - ... for small |z| (series, numerically
    safer near z=0 than the closed forms' 0/0 cancellation); closed forms
    otherwise."""
    if abs(z) < 1e-6:
        return 1.0 / 2.0 - z / 24.0 + z**2 / 720.0
    if z > 0:
        sq = z**0.5
        return (1.0 - math.cos(sq)) / z
    sq = (-z) ** 0.5
    return (1.0 - math.cosh(sq)) / z


def _stumpff_s(z: float) -> float:
    """S(z): 1/6 - z/120 + z^2/5040 - ... for small |z|; closed forms
    otherwise."""
    if abs(z) < 1e-6:
        return 1.0 / 6.0 - z / 120.0 + z**2 / 5040.0
    if z > 0:
        sq = z**0.5
        return (sq - math.sin(sq)) / sq**3
    sq = (-z) ** 0.5
    return (math.sinh(sq) - sq) / sq**3


@dataclass(frozen=True)
class LambertSolution:
    """Result of solving Lambert's problem for one departure/arrival
    position pair and a fixed time of flight."""

    p_km: float
    a_km: float
    e: float
    v1_km_s: tuple
    v2_km_s: tuple
    z: float
    iterations: int


def solve_lambert(
    r1_vec_km: tuple,
    r2_vec_km: tuple,
    tof_s: float,
    mu_km3_s2: float,
    prograde: bool = True,
    max_iter: int = 100,
    tol_s: float = 1e-6,
) -> LambertSolution:
    """Solve Lambert's problem: given two position vectors (km, in the
    orbital plane, as (x, y) tuples) and a time of flight (s), find the
    transfer orbit and the two velocity vectors at r1 and r2.

    Single-revolution, short-or-long way selected by `prograde` (True:
    counter-clockwise, the convention used consistently by this package's
    own transfer-angle definition in `api.py`).

    Raises ConvergenceError if the universal-variable Newton iteration
    does not converge within `max_iter` steps.
    """
    r1 = (r1_vec_km[0] ** 2 + r1_vec_km[1] ** 2) ** 0.5
    r2 = (r2_vec_km[0] ** 2 + r2_vec_km[1] ** 2) ** 0.5

    cross_z = r1_vec_km[0] * r2_vec_km[1] - r1_vec_km[1] * r2_vec_km[0]
    cos_dnu = (r1_vec_km[0] * r2_vec_km[0] + r1_vec_km[1] * r2_vec_km[1]) / (r1 * r2)
    cos_dnu = max(-1.0, min(1.0, cos_dnu))
    dnu = math.acos(cos_dnu)
    # Resolve the true (0, 2*pi) swept angle from the sign of the cross
    # product, consistent with a prograde (counter-clockwise) transfer.
    if prograde:
        if cross_z < 0:
            dnu = 2.0 * math.pi - dnu
    else:
        if cross_z > 0:
            dnu = 2.0 * math.pi - dnu

    # A true geometric degeneracy, independent of A's own magnitude: r1_vec
    # and r2_vec exactly colinear and same-pointing (dnu=0 or dnu=2*pi).
    # The orbital-plane normal (r1 x r2) is undefined there, so no single
    # transfer plane exists -- unlike the dnu=pi (Hohmann) case below,
    # which is a degeneracy of THIS solver's A-formula, not of the
    # underlying geometry (a transfer plane is well-defined at dnu=pi;
    # this package simply never routes that angle through this solver,
    # using the exact Hohmann closed form instead -- see api.py).
    if abs(cross_z) < 1e-12 * r1 * r2 and cos_dnu > 0:
        raise ConvergenceError(
            f"Lambert solver: r1_vec and r2_vec are colinear and "
            f"same-pointing (dnu={math.degrees(dnu):.6f} deg ~ 0 or 360 deg) "
            "-- the transfer plane is geometrically undefined; this "
            "degenerate zero-angle case is already excluded upstream by "
            "compute_transfer's own delta_theta normalization."
        )

    # A, in the numerically stable form A = sqrt(2*r1*r2)*cos(dnu/2).
    #
    # This is NOT the textbook A = sin(dnu)*sqrt(r1*r2/(1-cos(dnu))) typed
    # in literally -- that form subtracts two nearly-equal values
    # (1-cos(dnu)) when dnu is close to 0 or 2*pi, which is exactly this
    # package's problematic region (transfer angles near the 0/360 deg
    # boundary). The two are algebraically identical, by the half-angle
    # identities sin(dnu) = 2*sin(dnu/2)*cos(dnu/2) and 1-cos(dnu) =
    # 2*sin(dnu/2)^2, PROVIDED sin(dnu/2) is taken with its true sign
    # rather than |sin(dnu/2)|:
    #   sin(dnu)/sqrt(1-cos(dnu)) = 2 sin(dnu/2)cos(dnu/2) / (sqrt(2)*sin(dnu/2))
    #                             = sqrt(2)*cos(dnu/2)
    # (the sin(dnu/2) factor cancels cleanly with NO extra sign() needed,
    # because dnu/2 in (0,pi) for dnu in (0,2*pi) means sin(dnu/2) >= 0
    # throughout -- an earlier draft of this derivation incorrectly
    # inserted an extra sign(sin(dnu)) factor here, which double-flips
    # the sign for dnu>pi; this was caught and corrected during
    # independent review). This form instead has its own cancellation
    # near dnu=pi (cos(dnu/2)->0), which is this package's Hohmann case
    # -- already solved in closed form and routed around this solver
    # entirely at exactly dnu=pi (see module docstring); it is not a
    # universal fix, just the better trade-off for the region this
    # package actually exercises Lambert in (near 0/360 deg).
    A = (2.0 * r1 * r2) ** 0.5 * math.cos(dnu / 2.0)
    if abs(A) < 1e-9 * (r1 * r2) ** 0.5:
        raise ConvergenceError(
            f"Lambert solver: A~0 (dnu={math.degrees(dnu):.6f} deg is 0, 180, "
            "or 360 deg -- a degenerate or exactly-Hohmann transfer plane); "
            "this case is excluded upstream by compute_transfer's own "
            "delta_theta normalization and its exact Hohmann closed form."
        )

    sqrt_mu = mu_km3_s2**0.5
    # Hard domain bounds: C(z) has its first zero at z=4*pi^2 (the
    # single-rev/multi-rev boundary); z<0 is the hyperbolic branch, with
    # no comparable hard bound for this package's elliptical use case.
    # IMPORTANT: unlike the previous version of this function, these
    # bounds are never used to clamp an evaluation point -- an evaluation
    # at or beyond z_hi is simply treated as invalid (returns None,
    # exactly like y<0), so the bracketing search and root-finder below
    # only ever see genuinely evaluated points, never a value silently
    # projected onto the boundary.
    z_hi = 4.0 * math.pi**2

    def tof_of(z: float):
        if not (-4.0 * math.pi**2 < z < z_hi):
            return None
        C = _stumpff_c(z)
        if C <= 0.0:
            return None
        S = _stumpff_s(z)
        y = r1 + r2 + A * (z * S - 1.0) / (C**0.5)
        if y < 0:
            return None
        chi = (y / C) ** 0.5
        return (chi**3 * S + A * y**0.5) / sqrt_mu

    # tof(z) is monotonically increasing in z throughout the valid
    # (y>=0, C>0) region for this package's single-revolution elliptical
    # transfers (verified numerically across the full angular domain
    # during development -- not re-derived analytically here, but the
    # bracketing method below does not actually require global
    # monotonicity: it only needs ONE sign change in the bracket it
    # hands to brentq, which the expansion loop below searches for
    # explicitly and fails loudly if it cannot find). This replaces the
    # previous Newton-iteration-with-clamped-finite-difference approach,
    # which could have its derivative estimate biased by the clamp near
    # z_hi (confirmed during review: for transfer angles approaching 360
    # deg, Newton's very first step jumped onto the z_hi clamp and
    # needed ~9 further iterations pinned there before recovering) --
    # using the same proven bracket-expansion + bracketing-root-finder
    # pattern as `core.find_p_star` removes that failure mode entirely
    # rather than patching around it.
    f0 = tof_of(0.0)
    if f0 is None:
        raise ConvergenceError(
            f"Lambert solver: z=0 (parabolic) is itself invalid for this "
            f"geometry (r1={r1:.3f} km, r2={r2:.3f} km, "
            f"dnu={math.degrees(dnu):.3f} deg) -- cannot even start the "
            "bracket search; please report this as a bug."
        )
    target = tof_s

    def signed(z):
        t = tof_of(z)
        return None if t is None else (t - target)

    f_at_0 = f0 - target
    if abs(f_at_0) <= tol_s:
        z_lo, z_hi_bracket = 0.0, 0.0
    else:
        # Expand geometrically from z=0 in the direction implied by the
        # sign at z=0 (tof increasing in z), doubling the step each time,
        # until a sign change is found or the domain bound is reached.
        # IMPORTANT: the candidate step is explicitly clipped to stay
        # strictly inside the open domain (-4*pi**2, z_hi) on every
        # iteration, rather than letting a doubling step jump clean over
        # the last untested region next to the boundary -- an earlier
        # version of this loop could skip from z=32 straight to z=64,
        # missing a valid sign change that might exist between 32 and
        # the true boundary at 4*pi**2~39.48, found during independent
        # review. Clipping means the last few steps shrink instead of
        # doubling forever, which is exactly the point: it guarantees
        # every point up to (but not touching) the boundary gets tested
        # before giving up.
        direction_sign = 1.0 if f_at_0 < 0 else -1.0
        z_bound = z_hi if direction_sign > 0 else -4.0 * math.pi**2
        margin = 1e-9 * abs(z_bound) if z_bound != 0 else 1e-9
        step = 1.0
        z_prev, f_prev = 0.0, f_at_0
        z_lo = z_hi_bracket = None
        for _ in range(200):
            z_candidate = direction_sign * step
            z_next = max(-4.0 * math.pi**2 + margin, min(z_bound - margin, z_candidate))
            if z_next == z_prev:
                # The clip has pinned us to the same point already
                # tested (we've reached the boundary with no sign
                # change found anywhere in between) -- stop rather than
                # loop forever re-testing the same point.
                break
            f_next = signed(z_next)
            if f_next is not None and (f_prev is None or f_next * f_prev <= 0):
                lo, hi = sorted((z_prev, z_next))
                z_lo, z_hi_bracket = lo, hi
                break
            if f_next is not None:
                z_prev, f_prev = z_next, f_next
            step *= 2.0
        if z_lo is None:
            raise ConvergenceError(
                f"Lambert solver: could not bracket a root for tof_s={tof_s!r} s "
                f"(r1={r1:.3f} km, r2={r2:.3f} km, dnu={math.degrees(dnu):.3f} "
                "deg) within the valid domain of z -- this time of flight may "
                "be outside what a single-revolution transfer can achieve for "
                "this geometry."
            )

    def f(z):
        t = tof_of(z)
        if t is None:
            # z_lo and z_hi_bracket were both confirmed valid (non-None)
            # above, by construction of the expansion loop -- so brentq
            # must only ever probe this function at points strictly
            # between two valid, sign-changing evaluations. An invalid
            # interior point here would mean tof(z) is NOT the
            # well-behaved, single-valid-region function this bracketing
            # strategy assumes, which is a real anomaly worth surfacing
            # loudly rather than papering over with a signed-infinity
            # sentinel (as an earlier version of this function did).
            raise ConvergenceError(
                f"Lambert solver: internal bracket [{z_lo!r}, "
                f"{z_hi_bracket!r}] contains an invalid interior point at "
                f"z={z!r} even though both of its own endpoints were "
                "valid -- tof(z) is not well-behaved on this bracket; "
                "please report this as a bug rather than trust the result."
            )
        return t - target

    if z_lo == z_hi_bracket:
        z_root = z_lo
        iterations = 0
    else:
        try:
            z_root, root_info = brentq(
                f, z_lo, z_hi_bracket, xtol=1e-12, rtol=1e-14, maxiter=max_iter,
                full_output=True,
            )
        except ConvergenceError:
            raise  # already specific and informative (e.g. the internal
            # bracket-invalidity check inside f() above) -- do not
            # re-wrap it with a more generic message.
        except (ValueError, RuntimeError) as exc:
            raise ConvergenceError(
                f"Lambert solver's bracketing root-finder failed "
                f"(r1={r1:.3f} km, r2={r2:.3f} km, dnu={math.degrees(dnu):.3f} "
                f"deg, tof_s={tof_s!r}): {exc}"
            ) from exc
        if not root_info.converged:
            raise ConvergenceError(
                f"Lambert solver's bracketing root-finder did not converge "
                f"within {max_iter} iterations (r1={r1:.3f} km, r2={r2:.3f} "
                f"km, dnu={math.degrees(dnu):.3f} deg, tof_s={tof_s!r})."
            )
        iterations = root_info.iterations

    C = _stumpff_c(z_root)
    S = _stumpff_s(z_root)
    y = r1 + r2 + A * (z_root * S - 1.0) / (C**0.5)

    f_coef = 1.0 - y / r1
    g = A * (y / mu_km3_s2) ** 0.5
    gdot = 1.0 - y / r2

    v1x = (r2_vec_km[0] - f_coef * r1_vec_km[0]) / g
    v1y = (r2_vec_km[1] - f_coef * r1_vec_km[1]) / g
    v2x = (gdot * r2_vec_km[0] - r1_vec_km[0]) / g
    v2y = (gdot * r2_vec_km[1] - r1_vec_km[1]) / g

    # Semi-latus rectum from the specific angular momentum h = r x v
    # (more directly robust than re-deriving it from chi/C(z)):
    h = r1_vec_km[0] * v1y - r1_vec_km[1] * v1x
    p_km = h**2 / mu_km3_s2
    v1_sq = v1x**2 + v1y**2
    denom = 2.0 / r1 - v1_sq / mu_km3_s2
    if not (denom > 0):
        raise ConvergenceError(
            f"Lambert solver: vis-viva denominator (2/r1 - v1^2/mu) = "
            f"{denom!r} is not positive -- the recovered velocity implies a "
            "non-elliptical (parabolic or hyperbolic) orbit, outside this "
            "package's always-elliptical domain (see module docstring); "
            "refusing to silently report a>0 or e<1 in that case."
        )
    a_km = 1.0 / denom
    ratio = p_km / a_km
    # A tiny (1e-10 relative) overshoot of ratio above 1.0 is tolerated as
    # ordinary floating-point roundoff and clamped exactly to 1.0 (e=0,
    # a circular orbit) -- anything beyond that tolerance is a genuine
    # numerical anomaly and must be surfaced, not silently reinterpreted
    # as a plausible-looking circular orbit (an earlier version of this
    # line did exactly that, unconditionally: `e=0.0 if ratio>1 else
    # sqrt(1-ratio)`, found during independent review to be a silent
    # failure mode that could mask a real bug behind a falsely-clean e=0).
    if ratio > 1.0:
        if ratio - 1.0 > 1e-10:
            raise ConvergenceError(
                f"Lambert solver: p/a = {ratio!r} exceeds 1.0 by more than "
                "the tolerated floating-point roundoff (1e-10) -- this is "
                "not a circular orbit artifact, it indicates an actual "
                "inconsistency in the recovered solution; refusing to "
                "silently report e=0."
            )
        ratio = 1.0
    e = (1.0 - ratio) ** 0.5

    return LambertSolution(
        p_km=p_km,
        a_km=a_km,
        e=e,
        v1_km_s=(v1x, v1y),
        v2_km_s=(v2x, v2y),
        z=z_root,
        iterations=iterations,
    )
