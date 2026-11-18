"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: Theorem 4 (certified curvature sign change of the
curvature determinant H(p,phi) along the physical stationary branch,
at R=2 exactly). This module is a narrow diagnostic, not a general
solver: its scope is strictly R=2. The extension of Theorem 4's
sign-change picture to general R>1 is Conjecture 2
(`circular/curvature_sign_change_general_R/`), which remains an OPEN
CONJECTURE -- see `STATUS.md`. This module never extrapolates to R!=2.

WHAT THIS MODULE IS, IN ONE SENTENCE
-------------------------------------
For a transfer at R=2 exactly, classify whether the current angular
deviation phi sits before or after the certified curvature sign-change
point w_c (equivalently phi_c), i.e. whether the curvature determinant
H(p,phi) := C_phiphi*C_pp - C_pphi^2 is negative or positive there --
useful for a second-order sensitivity assessment of the cost to angular
deviation errors, at R=2.

WHERE w_c COMES FROM (NOT RECOMPUTED HERE)
---------------------------------------------
`w_c` below is not an independently rederived constant: it is copied,
literally, from the Krawczyk-certified isolated zero K_1 established in
`circular/theorem_4/03_central_mesh/` and assembled, with the other
four regions, in
`circular/theorem_4/06_global_audit/global_coverage_summary.md`. This
module's own test suite pins `W_C_CERTIFIED` against that certificate's
own recorded value, so that a future re-certification with a different
value is caught rather than silently diverging from this module's copy
(see `tests/test_curvature_regime.py`).

THIS MODULE DOES NOT RE-ESTABLISH THEOREM 4
-----------------------------------------------
It does not re-run any interval-arithmetic certification; it only (1)
solves for p* at the caller's (phi, R=2) via the package's own existing
solver (`core.solve`, already relied upon by `api.compute_transfer`),
(2) computes w := sqrt(p*), and (3) compares w to the already-certified
w_c. The classification is therefore only as trustworthy as `core.solve`
itself (already exercised across this package's whole test suite) plus
the one certified constant copied in below.
"""
import math
from dataclasses import dataclass
from typing import Literal

from . import core
from .errors import DomainError
from .units import EARTH_MU_KM3_S2

#: The Krawczyk-certified curvature sign-change point w_c := sqrt(p_c),
#: copied literally from `circular/theorem_4/03_central_mesh/
#: mesh_certificate_R2.py` (consistent with the value recorded in
#: `circular/theorem_4/06_global_audit/global_coverage_summary.md`).
#: NOT recomputed here -- see module docstring.
W_C_CERTIFIED = 0.30698081390187630871

#: phi_c in degrees, the angular deviation corresponding to W_C_CERTIFIED
#: on the physical stationary branch at R=2. This is NOT derivable from
#: W_C_CERTIFIED by a closed-form inversion (unlike budget_inversion.py's
#: z := cot^2(phi/2), the relation between w := sqrt(p) and phi here runs
#: through the actual stationary-branch solve p*(phi, R=2), which has no
#: simple closed form) -- so this is its own literal certified constant,
#: copied from the same source (`mesh_certificate_R2.py`'s own
#: `phi0_root = mp.mpf('164.3710517') * mp.pi / 180`), consistent with
#: `circular/theorem_4/06_global_audit/global_coverage_summary.md`.
PHI_C_CERTIFIED_DEG = 164.3710517

#: Relative tolerance for checking R == 2 exactly (floating-point
#: equality is never required of a caller-supplied r1_km/r2_km pair).
_R_EQUALS_2_RTOL = 1e-9


@dataclass(frozen=True)
class CurvatureRegime:
    """Classification of a R=2 transfer's position relative to the
    certified curvature sign-change point w_c."""

    r1_km: float
    r2_km: float
    delta_theta_deg: float
    mu_km3_s2: float

    phi_rad: float
    w: float
    p_star_km: float

    #: One of "pre_inflection" (H<0, w<w_c), "at_inflection" (w==w_c
    #: within `w_tol`), "post_inflection" (H>0, w>w_c).
    regime: Literal["pre_inflection", "at_inflection", "post_inflection"]
    distance_to_wc: float  # w - W_C_CERTIFIED, signed
    w_c_reference: float = W_C_CERTIFIED
    phi_c_reference_deg: float = PHI_C_CERTIFIED_DEG


def classify_curvature_regime(
    r1_km: float,
    r2_km: float,
    delta_theta_deg: float,
    mu_km3_s2: float = EARTH_MU_KM3_S2,
    w_tol: float = 1e-9,
) -> CurvatureRegime:
    """Classify a transfer at R=2 as pre-, at-, or post-inflection
    relative to Theorem 4's certified curvature sign-change point.

    Parameters
    ----------
    r1_km, r2_km : departure and arrival circular-orbit radii, in km.
        Must satisfy r2_km / r1_km == 2 (within `_R_EQUALS_2_RTOL`) --
        this is a hard requirement, not a default: see `DomainError`
        below.
    delta_theta_deg : transfer angle actually traversed, in degrees,
        strictly between 0 and 360 (mod 360). Same convention as
        `api.compute_transfer`.
    mu_km3_s2 : gravitational parameter of the central body, km^3/s^2.
    w_tol : tolerance, in units of w := sqrt(p) (dimensionless, p in
        units of r1), within which w is reported as "at_inflection"
        rather than strictly pre- or post-.

    Raises
    ------
    DomainError
        If R := r2_km/r1_km is not 2 (within tolerance): Theorem 4 is
        proven for R=2 only; its extension to general R>1 is the OPEN
        Conjecture 2 (`circular/curvature_sign_change_general_R/`).
        This function never extrapolates to another R.
        Also raised for invalid radii or the degenerate zero-angle
        transfer, same convention as `api.compute_transfer`.
    """
    if not (r1_km > 0 and r2_km > 0):
        raise DomainError(f"radii must be positive: r1_km={r1_km!r}, r2_km={r2_km!r}")
    R = r2_km / r1_km
    if abs(R - 2.0) > _R_EQUALS_2_RTOL * 2.0:
        raise DomainError(
            f"this diagnostic is certified for R=2 only (got R={R!r} from "
            f"r1_km={r1_km!r}, r2_km={r2_km!r}): Theorem 4's curvature "
            "sign change is proven at R=2 exactly; its extension to "
            "general R>1 is Conjecture 2 "
            "(`circular/curvature_sign_change_general_R/`), which remains "
            "an OPEN CONJECTURE (see STATUS.md). This function never "
            "extrapolates beyond the certified R=2 case."
        )

    folded = delta_theta_deg % 360.0
    if folded == 0.0:
        raise DomainError(
            f"delta_theta_deg={delta_theta_deg!r} is the degenerate zero-angle "
            "transfer; the physical domain is 0 < delta_theta < 360 deg."
        )
    phi = math.pi - math.radians(folded)
    phi_abs = abs(phi)

    if math.isclose(phi_abs, 0.0, abs_tol=1e-12):
        # Exact Hohmann: p* = hohmann.solve_hohmann(2.0).p in closed form.
        from . import hohmann
        p_star = hohmann.solve_hohmann(2.0).p
    else:
        p_star = core.find_p_star(phi_abs, 2.0)

    w = p_star ** 0.5
    distance = w - W_C_CERTIFIED

    if abs(distance) <= w_tol:
        regime: Literal["pre_inflection", "at_inflection", "post_inflection"] = "at_inflection"
    elif distance < 0.0:
        regime = "pre_inflection"
    else:
        regime = "post_inflection"

    return CurvatureRegime(
        r1_km=r1_km, r2_km=r2_km, delta_theta_deg=delta_theta_deg,
        mu_km3_s2=mu_km3_s2,
        phi_rad=phi, w=w, p_star_km=p_star * r1_km,
        regime=regime, distance_to_wc=distance,
    )
