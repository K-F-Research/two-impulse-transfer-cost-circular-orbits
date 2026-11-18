"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: none of the manuscript's own numbered theorems
directly. The correctness of the root-finder below rests on the
repository's supplementary existence/uniqueness result
(`circular/existence_uniqueness_minimizer/` -- not a numbered lemma of
the manuscript; see `STATUS.md`). Theorem 1 (global monotonicity in
phi) and Theorem 2 (local curvature at Hohmann) are not needed here:
this module solves for p* at a single, fixed phi and does not depend
on how the resulting cost varies as phi changes.

Normalized-unit physics core of 2i-OTCA-C.

All quantities here are in the repository's own normalized units
(`docs/mathematical_conventions.md`): `r1 = 1`, `mu = 1`, `R := r2/r1 > 1`,
`phi := pi - Delta_theta in (0, pi)` the transfer-angle variable (Hohmann is
`phi = 0`). `p` is the semi-latus rectum of the transfer orbit, `p* (phi,R)`
its minimizing value (global existence and uniqueness established by a supplementary
repository result -- not a numbered lemma of the manuscript -- in
`circular/existence_uniqueness_minimizer/`; see `STATUS.md`).

This module reproduces, without modification, the exact primitive cost
definitions verified independently in
`circular/global_monotonicity/verify_monotonicity_independent.py`
(`transfer_cost`) -- the same closed-form `S1(p,phi,R)`, `S2(p,phi,R)`,
`C := sqrt(S1) + sqrt(S2)` used throughout this repository's proofs. The
symbolic derivative `C_p` needed for the root-finder below is derived here
directly from that same definition with SymPy (not hand-transcribed), then
compiled to an ordinary fast numeric function -- the same pattern used by
`circular/inflection_uniqueness_R2/mesh_certificate_R2.py` and
`krawczyk_certificate_R2.py`, so a reader who has already checked those
proof scripts recognizes the construction immediately.

Why the root-finder is safe (not just "it usually works"): the
existence/uniqueness result (`circular/existence_uniqueness_minimizer/
lemma1_closure_stationary_point.py` -- a supplementary repository result,
not a numbered lemma of the manuscript; see `STATUS.md`)
proves `C_p(p,phi,R) < 0` as `p -> 0+` and `C_p(p,phi,R) -> +infinity` as
`p -> +infinity`, for every `phi in (0,pi)` and `R>1` -- so a sign change,
and by the same lemma's interior argument exactly one, is guaranteed to
exist in `(0, +infinity)`. `find_p_star` below verifies this sign change
numerically on a generous bracket before calling the bracketed root-finder,
and raises `ConvergenceError` (never a silent wrong answer) if the expected
sign pattern is not observed for the given input.
"""
import math
from dataclasses import dataclass

import numpy as np
import sympy as sp
from scipy.optimize import brentq

from .errors import ConvergenceError

_p, _phi, _R = sp.symbols("p phi R", positive=True)
_c = sp.cos(_phi)
_s = sp.sin(_phi)
_x = _p - 1
_y = _p / _R - 1
_e2 = _x**2 + (_x * _c + _y) ** 2 / _s**2
_v1_sq = 2 + (_e2 - 1) / _p
_v2_sq = 2 / _R + (_e2 - 1) / _p
_v_theta1 = sp.sqrt(_p)
_v_theta2 = sp.sqrt(_p) / _R
_v_c1 = sp.Integer(1)
_v_c2 = 1 / sp.sqrt(_R)
_S1 = _v1_sq + _v_c1**2 - 2 * _v_c1 * _v_theta1
_S2 = _v2_sq + _v_c2**2 - 2 * _v_c2 * _v_theta2
_C = sp.sqrt(_S1) + sp.sqrt(_S2)
_Cp = sp.diff(_C, _p)
_Cphi = sp.diff(_C, _phi)
_e = sp.sqrt(_e2)

# Compiled to ordinary (non-interval) floating point for speed: this module
# is the operational/engineering component, not a certification script --
# the certificates for the specific case R=2 live in
# `circular/inflection_uniqueness_R2/` and use interval arithmetic instead.
_f_S1 = sp.lambdify((_p, _phi, _R), _S1, "numpy")
_f_S2 = sp.lambdify((_p, _phi, _R), _S2, "numpy")
_f_C = sp.lambdify((_p, _phi, _R), _C, "numpy")
_f_Cp = sp.lambdify((_p, _phi, _R), _Cp, "numpy")
_f_Cphi = sp.lambdify((_p, _phi, _R), _Cphi, "numpy")
_f_e = sp.lambdify((_p, _phi, _R), _e, "numpy")
_f_x = sp.lambdify((_p, _phi, _R), _x, "numpy")
_f_y = sp.lambdify((_p, _phi, _R), _y, "numpy")


@dataclass(frozen=True)
class NormalizedSolution:
    """The transfer, fully described in normalized units (r1=1, mu=1)."""

    phi: float
    R: float
    p_star: float
    e: float
    a: float
    S1: float
    S2: float
    C: float  # total cost = sqrt(S1) + sqrt(S2) = delta_v1 + delta_v2, in units of vc1


def cost(p: float, phi: float, R: float) -> float:
    """C(p,phi,R) = sqrt(S1) + sqrt(S2), the (not-yet-minimized) total
    two-impulse cost at a given semi-latus rectum p."""
    return float(_f_C(p, phi, R))


def cost_derivative(p: float, phi: float, R: float) -> float:
    """C_p(p,phi,R) = dC/dp, used by the root-finder; exposed for testing
    and for callers who want to inspect the stationary-point condition
    directly."""
    return float(_f_Cp(p, phi, R))


def cost_partial_phi(p: float, phi: float, R: float) -> float:
    """C_phi(p,phi,R) = partial C/partial phi at a FIXED p (not the total
    derivative dV/dphi of the optimized cost V_R(phi) = min_p C(p,phi,R)).

    Exposed for `budget_inversion.py`: by the envelope theorem, since
    C_p(p*,phi,R) = 0 at the optimizer p* := find_p_star(phi,R), the total
    derivative of V_R at phi equals this *partial* derivative evaluated at
    p=p*:

        dV_R/dphi (phi) = C_phi(p*(phi,R), phi, R).

    This lets a caller get dV_R/dphi exactly (from the same symbolic
    derivation as C_p above, not a finite difference) without having to
    differentiate through the implicit p*(phi,R) itself.
    """
    return float(_f_Cphi(p, phi, R))


def find_p_star(
    phi: float,
    R: float,
    p_lo: float = 1e-9,
    p_hi: float = 1e9,
    xtol: float = 1e-13,
    max_expansions: int = 300,
) -> float:
    """Find the unique minimizer p*(phi,R) of C(.,phi,R) by bracketed
    root-finding on C_p = 0 (Brent's method), after numerically confirming
    the sign pattern the existence/uniqueness result guarantees (C_p<0 near p_lo, C_p>0 near p_hi).

    `p_lo`/`p_hi` are only a STARTING bracket, not a fixed one: a
    reproducible failure was found during independent review for phi very
    close to pi (equivalently delta_theta very close to 0 or 360 deg) --
    e.g. phi=3.1415752..., R=2.0 (delta_theta=359.999 deg) has its true
    p* near 4e-10, well below the previous fixed p_lo=1e-9, so
    C_p(1e-9)=+6.3e8 was already positive (the WRONG sign) rather than
    negative as that existence/uniqueness result predicts near p->0+. This is
    not a failure of that result (which guarantees the sign pattern holds in the limit
    p->0+/p->+infinity for every phi in the open interval, not at any
    particular fixed finite p), it was a failure of assuming one fixed
    finite bracket is close enough to those limits for every phi and R
    this package can be asked to solve.

    This function now geometrically CONTRACTS p_lo toward 0 and EXPANDS
    p_hi toward +infinity (by repeated factors of 10 in each direction,
    independently) until the sign pattern is observed, each step checking
    that the evaluation is finite (not NaN/inf) before accepting it as a
    valid bracket endpoint -- an overflow or NaN is treated the same as a
    wrong sign (keep searching, do not accept it as "the" bracket). The
    search is bounded (`max_expansions` per direction, hard floor/ceiling
    at 1e-300/1e300) and raises ConvergenceError explicitly, distinguishing
    "no sign change found" from "an expansion limit was reached", rather
    than silently widening forever or accepting a non-finite endpoint.
    """
    hard_floor, hard_ceiling = 1e-300, 1e300

    def finite_sign(p: float):
        """Returns C_p(p) if finite, else None (never a fabricated sign).
        Division/invalid-value warnings are expected and already handled
        here (via the isfinite check below) when this probes an actual
        domain boundary (e.g. phi=0 exactly) -- silenced locally rather
        than left to print spurious noise for entirely legitimate,
        correctly-handled inputs."""
        with np.errstate(divide="ignore", invalid="ignore"):
            v = _f_Cp(p, phi, R)
        return float(v) if math.isfinite(v) else None

    g_lo = finite_sign(p_lo)
    expansions_lo = 0
    while (g_lo is None or g_lo >= 0) and p_lo > hard_floor and expansions_lo < max_expansions:
        p_lo /= 10.0
        g_lo = finite_sign(p_lo)
        expansions_lo += 1
    if g_lo is None or g_lo >= 0:
        raise ConvergenceError(
            "could not find p_lo with a finite, negative C_p (the "
            "existence/uniqueness result's "
            f"p->0+ limit) for phi={phi!r}, R={R!r} after {expansions_lo} "
            f"contractions down to p_lo={p_lo!r} (hard floor {hard_floor!r}); "
            "this input may be on or past a genuine domain boundary."
        )

    g_hi = finite_sign(p_hi)
    expansions_hi = 0
    while (g_hi is None or g_hi <= 0) and p_hi < hard_ceiling and expansions_hi < max_expansions:
        p_hi *= 10.0
        g_hi = finite_sign(p_hi)
        expansions_hi += 1
    if g_hi is None or g_hi <= 0:
        raise ConvergenceError(
            "could not find p_hi with a finite, positive C_p (the "
            "existence/uniqueness result's "
            f"p->+infinity limit) for phi={phi!r}, R={R!r} after "
            f"{expansions_hi} expansions up to p_hi={p_hi!r} (hard ceiling "
            f"{hard_ceiling!r}); this input may be on or past a genuine "
            "domain boundary."
        )

    return brentq(lambda p: _f_Cp(p, phi, R), p_lo, p_hi, xtol=xtol, rtol=1e-14)


def solve(phi: float, R: float) -> NormalizedSolution:
    """Solve the full normalized-unit transfer problem at a given (phi, R):
    find p*, and return every quantity needed downstream (eccentricity,
    semi-major axis, individual and total costs)."""
    p_star = find_p_star(phi, R)
    e = float(_f_e(p_star, phi, R))
    a = p_star / (1 - e**2)
    S1 = float(_f_S1(p_star, phi, R))
    S2 = float(_f_S2(p_star, phi, R))
    C = (S1 ** 0.5) + (S2 ** 0.5)
    return NormalizedSolution(phi=phi, R=R, p_star=p_star, e=e, a=a, S1=S1, S2=S2, C=C)
