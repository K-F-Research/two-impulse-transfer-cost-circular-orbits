"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: Theorem 2 of the manuscript (the closed-form
local curvature coefficient K(r1,r2) at the Hohmann configuration; see
below for the exact derivation this module reproduces).

Fast analytic surrogate for the two-impulse transfer cost, near the Hohmann
configuration (Theorem 2 of the accompanying manuscript).

WHAT THIS MODULE IS, IN ONE SENTENCE
-------------------------------------
`core.solve` (used by `api.compute_transfer`) finds the exact minimum-cost
transfer by numerically solving for the optimal semi-latus rectum p* at a
given transfer angle. That is correct for every angle in the model's whole
domain, but it is a numerical minimization, not a closed-form evaluation.
This module instead evaluates a *closed-form quadratic approximation*,

    V_R(phi) ~= V_R(0) + K(r1,r2) * phi**2,

which is cheap (one algebraic expression, no root-finding) but, unlike
`core.solve`, is only accurate near phi=0 (the Hohmann configuration) --
see `MAX_VALIDATED_PHI_DEG` below for the precise, honest boundary of
"near".

WHERE K(r1,r2) COMES FROM
--------------------------
`K(r1,r2)` is Theorem 2's closed-form local curvature coefficient, derived
in the manuscript's Section 4.6-4.7 and tracked in
`docs/mathematical_conventions.md` at the repository root (the Taylor
coefficient, `K(r1,r2) = (1/2) V_R''(0)` -- NOT the bare second derivative
`K_circ(R) = V_R''(0)`; the two differ by exactly a factor of 2, as that
glossary states). The closed form below (`_K_normalized`) reproduces,
term for term, the derivation of Section 4.7: with `t := sqrt(R)`,
`u := sqrt(1+t**2)`,

    N := 2*sqrt(2)*t + (1+t)*u
    M := sqrt(2)*t*(t+1) + u**3
    a_K := u**7 * N / (4 * t**6 * (t**2-1))
    c_K := (t**2-1) * N / (16 * t**2 * u)
    b_K := (t-1) * u * M / (4 * t**4)
    K(1,R) := c_K - b_K**2 / (4*a_K)

This is NOT a numerical fit: it is the exact algebraic expression the
manuscript proves equals the true Taylor coefficient, for every R>1. It
was checked against the manuscript's own worked table (Section 4.8:
K(1,R) at R=1.2, 1.5, 2, 3, 5, 10) in `tests/test_surrogate.py`.

K(r1,r2) FOR GENERAL (NOT NORMALIZED) INPUTS
----------------------------------------------
The manuscript's Lemmas 4.4-4.5 establish that K(r1,r2) depends on the
pair (r1,r2) only through R:=r2/r1 and an overall speed scale: at fixed
mu, K(c*r1, c*r2) = K(r1,r2) / sqrt(c). Combined with the model's own
normalized units (`core.py`, `units.py`: speeds are multiples of
vc1 = sqrt(mu/r1)), this means the real-unit curvature coefficient is

    K_km_s(r1_km, r2_km, mu_km3_s2) = K(1, R) * sqrt(mu_km3_s2 / r1_km)

-- exactly the same "multiply the normalized quantity by vc1" pattern
`api.compute_transfer` already uses for delta_v1_km_s/delta_v2_km_s (see
`units.circular_speed_km_s`), applied here to a curvature coefficient
instead of a cost.

THE HONEST LIMIT OF THIS SURROGATE (READ BEFORE USING)
---------------------------------------------------------
The quadratic model is a *truncated* Taylor expansion; the manuscript's
own Section 7.2 benchmark (`circular/hohmann_local_minimum/
surrogate_benchmark/benchmark_quadratic_surrogate.py`) measured its
relative error ONLY at R=1.5, for angular deviations from 5 deg to 30 deg,
finding 0.0003% at 5 deg growing to 0.24% at 30 deg (consistent with the
theorem's own O(phi**4) remainder). That benchmark is the ONLY
quantitative evidence this repository has for "how good is the
approximation, at what angle" -- it was not re-measured here for other
values of R. Accordingly:

- `MAX_VALIDATED_PHI_DEG = 30.0` is used as a conservative default warning
  threshold, for ANY R, not because the error is proven to stay small up
  to 30 deg for every R (it is not), but because 30 deg is the largest
  angle this repository has actually measured, for one representative R.
- `estimate_cost` always reports whether the requested angle is within
  this threshold, and never raises -- it is the caller's decision how to
  treat a flagged, out-of-validated-range estimate, exactly like
  `TransferResult.numerical_warnings` in `api.py` is a warning, not an
  exception.
- `compare_to_exact` is provided specifically so a caller can verify the
  surrogate's actual accuracy at their own (r1,r2,phi), against this
  package's own certified-domain exact solver (`core.solve`), rather than
  trusting the 30-degree default blindly for a case this repository never
  measured.
"""
import math
from dataclasses import dataclass
from typing import Optional

from . import core, hohmann
from .errors import DomainError
from .units import EARTH_MU_KM3_S2, circular_speed_km_s

#: Largest angular deviation from Hohmann, in degrees, for which this
#: repository has direct, measured evidence of the surrogate's accuracy
#: (at R=1.5; see this module's docstring). Not a proven universal bound
#: for every R -- a conservative, explicitly-sourced default threshold.
MAX_VALIDATED_PHI_DEG = 30.0


def _K_normalized(R: float) -> float:
    """K(1,R), the Taylor coefficient of Theorem 2, via the closed form of
    the manuscript's Section 4.6-4.7 (see module docstring). Normalized
    units: r1=1, mu=1."""
    if not R > 1.0:
        raise DomainError(
            f"K(1,R) is established by Theorem 2 for R>1 only (r1=1 is the "
            f"smaller radius by this component's convention); got R={R!r}. "
            "Use Lemma 4.4 (K(r1,r2)=K(r2,r1)) yourself if you need R<1."
        )
    t = R ** 0.5
    u = (1.0 + t * t) ** 0.5
    N = 2.0 * (2.0 ** 0.5) * t + (1.0 + t) * u
    M = (2.0 ** 0.5) * t * (1.0 + t) + u ** 3
    a_K = (u ** 7) * N / (4.0 * t ** 6 * (t * t - 1.0))
    c_K = (t * t - 1.0) * N / (16.0 * t * t * u)
    b_K = (t - 1.0) * u * M / (4.0 * t ** 4)
    K = c_K - (b_K * b_K) / (4.0 * a_K)

    # Defensive check, not a physical possibility if this function is
    # correct: Theorem 2 proves K(1,R) > 0 for every R>1 (manuscript
    # Section 4.7, "manifest positivity" of 4*a_K*c_K - b_K**2). A
    # non-finite or non-positive result here means this implementation
    # has a bug (e.g. a sign or algebra error introduced by a future
    # edit), not a legitimate edge case of the model -- raised rather
    # than silently returned, following this package's own
    # never-return-a-silently-wrong-value discipline (see errors.py).
    if not (math.isfinite(K) and K > 0.0):
        raise DomainError(
            f"internal consistency check failed: K(1,R={R!r}) evaluated to "
            f"{K!r}, but Theorem 2 proves this must be a finite, strictly "
            "positive number for every R>1 (manuscript Section 4.7). This "
            "indicates a defect in this function's implementation, not a "
            "valid result -- please report it rather than using this value."
        )
    return K


def curvature_coefficient_km_s(
    r1_km: float, r2_km: float, mu_km3_s2: float = EARTH_MU_KM3_S2
) -> float:
    """K(r1,r2) in real units -- the quadratic Taylor coefficient of
    Theorem 2.

    Parameters
    ----------
    r1_km, r2_km : float
        Departure and arrival circular-orbit radii, in km, measured from
        the central body's centre. Requires ``r2_km > r1_km`` (R>1 is the
        domain Theorem 2 is proven on; see `api.compute_transfer` for the
        same convention and how to handle an inbound transfer).
    mu_km3_s2 : float, optional
        Gravitational parameter of the central body, km^3/s^2. Defaults
        to Earth (`otca_c.units.Bodies.EARTH`).

    Returns
    -------
    float
        K(r1,r2), in km/s per rad**2: the coefficient such that
        ``V_R(phi) ~= V_R(0) + K(r1,r2) * phi**2`` for phi near 0 (see
        module docstring for the exact closed form and its domain of
        validity).

    Raises
    ------
    DomainError
        If the radii are non-positive, or r2_km does not exceed r1_km.
    """
    if not (r1_km > 0 and r2_km > 0):
        raise DomainError(f"radii must be positive: r1_km={r1_km!r}, r2_km={r2_km!r}")
    if not r2_km > r1_km:
        raise DomainError(
            f"this function requires r2_km > r1_km (R:=r2/r1>1, the domain "
            f"Theorem 2 is proven on); got r1_km={r1_km!r}, r2_km={r2_km!r}."
        )
    R = r2_km / r1_km
    return _K_normalized(R) * circular_speed_km_s(r1_km, mu_km3_s2)


@dataclass(frozen=True)
class SurrogateEstimate:
    """A fast, closed-form estimate of the two-impulse transfer cost near
    the Hohmann configuration (Theorem 2's quadratic surrogate), together
    with everything a caller needs to judge whether to trust it."""

    r1_km: float
    r2_km: float
    delta_theta_deg: float
    phi_deg: float
    mu_km3_s2: float

    K_km_s: float  # curvature_coefficient_km_s(r1_km, r2_km, mu_km3_s2)
    hohmann_delta_v_total_km_s: float
    delta_v_total_estimate_km_s: float

    #: True iff |phi_deg| <= MAX_VALIDATED_PHI_DEG (see module docstring:
    #: this is a conservative default threshold measured at one R, not a
    #: universal proven bound for every R).
    within_measured_domain: bool
    domain_note: str


def estimate_cost(
    r1_km: float,
    r2_km: float,
    delta_theta_deg: float,
    mu_km3_s2: float = EARTH_MU_KM3_S2,
) -> SurrogateEstimate:
    """Fast closed-form estimate of the two-impulse transfer cost, via
    Theorem 2's quadratic surrogate V_R(phi) ~= V_R(0) + K(r1,r2)*phi**2.

    This function does NOT run `core.solve`'s root-finder: it is cheap (a
    handful of algebraic operations), but only accurate near the Hohmann
    angle (delta_theta_deg near 180 deg) -- see the module docstring for
    exactly how near this repository has actually measured. It never
    raises merely for being far from Hohmann: it always returns an
    estimate, flagged via ``within_measured_domain`` and ``domain_note``;
    use `compare_to_exact` to check the real error at your own inputs
    rather than trusting that flag alone for anything beyond a quick
    sanity check.

    Parameters
    ----------
    r1_km, r2_km : float
        Departure and arrival circular-orbit radii, in km. Same
        convention as `api.compute_transfer` (requires ``r2_km >
        r1_km``); the two functions are safe to call side by side on the
        same inputs.
    delta_theta_deg : float
        Transfer angle actually traversed, in degrees, strictly between
        0 and 360 (mod 360). 180 deg is exactly Hohmann.
    mu_km3_s2 : float, optional
        Gravitational parameter of the central body, km^3/s^2. Defaults
        to Earth.

    Returns
    -------
    SurrogateEstimate
        The estimated total delta-v, the K used to get it, and an
        explicit flag plus explanatory note on whether this angle is
        within this repository's measured-accuracy range.

    Raises
    ------
    DomainError
        If the radii are invalid, or delta_theta_deg is the degenerate
        zero-angle transfer.
    """
    if not (r1_km > 0 and r2_km > 0):
        raise DomainError(f"radii must be positive: r1_km={r1_km!r}, r2_km={r2_km!r}")
    if not r2_km > r1_km:
        raise DomainError(
            f"this function requires r2_km > r1_km (R:=r2/r1>1); got "
            f"r1_km={r1_km!r}, r2_km={r2_km!r}."
        )

    folded = delta_theta_deg % 360.0
    if folded == 0.0:
        raise DomainError(
            f"delta_theta_deg={delta_theta_deg!r} is the degenerate zero-angle "
            "transfer; the physical domain is 0 < delta_theta < 360 deg."
        )
    phi_rad = math.pi - math.radians(folded)
    phi_deg = math.degrees(phi_rad)

    R = r2_km / r1_km
    K_norm = _K_normalized(R)
    vc1 = circular_speed_km_s(r1_km, mu_km3_s2)
    K_km_s = K_norm * vc1

    hohmann_ref = hohmann.solve_hohmann(R)
    hohmann_total_km_s = hohmann_ref.delta_v_total * vc1

    estimate_km_s = hohmann_total_km_s + K_km_s * (phi_rad ** 2)

    within = abs(phi_deg) <= MAX_VALIDATED_PHI_DEG
    note = (
        f"|phi|={abs(phi_deg):.3f} deg is within this repository's measured "
        f"accuracy range (<= {MAX_VALIDATED_PHI_DEG:.1f} deg, measured at "
        "R=1.5 -- see module docstring); the surrogate's actual error at "
        "your own R has not been separately measured."
        if within
        else
        f"|phi|={abs(phi_deg):.3f} deg EXCEEDS this repository's measured "
        f"accuracy range (<= {MAX_VALIDATED_PHI_DEG:.1f} deg, measured at "
        "R=1.5 only): this estimate is a mathematically valid evaluation "
        "of the Theorem 2 surrogate, but its error at this angle has not "
        "been measured and may be large (the surrogate's own remainder "
        "term grows as phi**4). Use compare_to_exact(), or api.compute_"
        "transfer() directly, before relying on this number."
    )

    return SurrogateEstimate(
        r1_km=r1_km,
        r2_km=r2_km,
        delta_theta_deg=delta_theta_deg,
        phi_deg=phi_deg,
        mu_km3_s2=mu_km3_s2,
        K_km_s=K_km_s,
        hohmann_delta_v_total_km_s=hohmann_total_km_s,
        delta_v_total_estimate_km_s=estimate_km_s,
        within_measured_domain=within,
        domain_note=note,
    )


@dataclass(frozen=True)
class SurrogateComparison:
    """Side-by-side comparison of the fast surrogate (`estimate_cost`)
    against this package's own certified-domain exact solver
    (`api.compute_transfer`), at the same inputs."""

    estimate: SurrogateEstimate
    exact_delta_v_total_km_s: float
    relative_error: float  # |estimate - exact| / exact, dimensionless


def compare_to_exact(
    r1_km: float,
    r2_km: float,
    delta_theta_deg: float,
    mu_km3_s2: float = EARTH_MU_KM3_S2,
) -> SurrogateComparison:
    """Run BOTH `estimate_cost` (this module) and `api.compute_transfer`
    (the exact solver) at the same inputs, and report the surrogate's
    actual relative error here, now, for this specific case -- rather
    than relying on the generic 30-degree default threshold, which this
    repository only measured at one representative R (see module
    docstring). Recommended before trusting the fast surrogate for any
    R, or any angle, this repository's own benchmark did not already
    cover.

    Parameters
    ----------
    r1_km, r2_km, delta_theta_deg, mu_km3_s2
        Same convention as `estimate_cost` and `api.compute_transfer`.

    Returns
    -------
    SurrogateComparison
        The surrogate's own result, the exact solver's total delta-v,
        and the surrogate's measured relative error at this exact input
        (not an assumption carried over from a different case).

    Raises
    ------
    DomainError
        If either underlying call raises it (invalid radii/angle).
    ConvergenceError
        If the exact solver (`core.find_p_star`) unexpectedly fails to
        converge for this input, despite it being inside the proven
        domain -- see `core.py`.

    Notes
    -----
    `validate=False` is passed to `api.compute_transfer` here: this
    function's own purpose already is an independent cross-check (of the
    surrogate against the exact solver), so `api.compute_transfer`'s own
    internal geometry self-check is a redundant, not a required,
    assurance for this specific call. Request the full `compute_transfer`
    result directly (with the default `validate=True`, or
    `lambert_check=True`/`rk4_check=True`) if the exact value's own
    cross-checks are also needed, rather than only the comparison.
    """
    # Imported here, not at module level, to avoid a hard dependency from
    # this lightweight module on api.py's larger import surface
    # (geometry/rk4/derived) for callers who only ever need the surrogate.
    from .api import compute_transfer

    est = estimate_cost(r1_km, r2_km, delta_theta_deg, mu_km3_s2)
    exact = compute_transfer(
        r1_km, r2_km, delta_theta_deg, mu_km3_s2, validate=False
    )
    exact_total = exact.delta_v_total_km_s
    rel_err = abs(est.delta_v_total_estimate_km_s - exact_total) / exact_total
    return SurrogateComparison(
        estimate=est, exact_delta_v_total_km_s=exact_total, relative_error=rel_err
    )
