"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: the repository's supplementary existence/
uniqueness result (via `core.py`, not a numbered lemma of the
manuscript -- see `STATUS.md`) and Theorem 2 (via `hohmann.py`'s exact
phi=0 closed form). Theorem 1 (global monotonicity in phi) is not
needed for the correctness of a single solve; see this module's own
`compute_transfer` docstring for exactly what domain claim rests on
what result.

High-level, physical-units public API for 2i-OTCA-C. See the package
docstring (`__init__.py`) and this folder's `README.md` for the full
mathematical contract: exact domain of validity, angle convention, and the
theorems this component relies on.
"""
import math
from dataclasses import dataclass
from typing import Optional, Tuple

from . import core, derived, hohmann
from .derived import AltitudeInfo, LambertCrossCheck, PropellantInfo
from .errors import ConvergenceError, DomainError
from .geometry import GeometryCheck, independent_check
from .rk4 import RK4Validation, rk4_validate
from .units import EARTH_MU_KM3_S2, circular_speed_km_s


@dataclass(frozen=True)
class TransferResult:
    """A fully solved two-impulse circular-to-circular transfer, in real
    physical units."""

    # Inputs, echoed back for traceability.
    r1_km: float
    r2_km: float
    delta_theta_deg: float
    mu_km3_s2: float

    # Derived model parameters (dimensionless).
    R: float
    phi_rad: float
    is_hohmann: bool

    # Transfer-orbit geometry.
    p_star_km: float
    semi_major_axis_km: float
    eccentricity: float

    # Costs.
    delta_v1_km_s: float
    delta_v2_km_s: float
    delta_v_total_km_s: float

    # Reference and derived comparison quantities.
    hohmann_delta_v_total_km_s: float
    penalty_vs_hohmann_pct: float

    # Independent cross-check (None only if validate=False was passed).
    geometry_check: Optional[GeometryCheck]

    # Unconditional derived quantities (always computed -- need nothing
    # beyond the orbit already solved and the mu already supplied).
    nu1_rad: float
    nu2_rad: float
    r_periapsis_km: float
    r_apoapsis_km: float
    orbital_period_s: float
    time_of_flight_s: float

    # Optional derived quantities -- None unless the caller explicitly
    # supplied the extra physical input they require (see
    # `compute_transfer`'s docstring: never guessed or defaulted).
    altitude_info: Optional[AltitudeInfo] = None
    propellant_info_departure: Optional[PropellantInfo] = None
    propellant_info_arrival: Optional[PropellantInfo] = None

    # Independent cross-check against a from-scratch Lambert-problem
    # solve (a genuinely separate classical method -- see
    # `lambert.py`/`derived.lambert_cross_check`), as opposed to
    # `geometry_check`'s internal textbook-consistency re-derivation.
    # None unless `lambert_check=True` was passed (it is not run by
    # default: it is a from-scratch root-find, meaningfully more
    # expensive than `geometry_check`'s closed-form re-derivation, and
    # this package does not spend that cost on every call unasked).
    lambert_check: Optional[LambertCrossCheck] = None

    # Independent cross-check against direct RK4 numerical propagation of
    # the two-body equations of motion from this transfer's own departure
    # state vector (see `rk4.py`'s module docstring for exactly how this
    # differs from both `geometry_check` and `lambert_check`). None
    # unless `rk4_check=True` was passed (off by default: it is the most
    # expensive of this package's three cross-checks, an initial-value
    # numerical integration run twice -- at two step counts -- for a
    # convergence study on every call).
    rk4_check: Optional[RK4Validation] = None

    # Per-output numerical-conditioning warnings (NOT a physical-validity
    # error -- those are already raised as exceptions elsewhere; this is
    # the repository's agreed three-level structure: (1) physical domain
    # validity, enforced by raising; (2)/(3) numerical conditioning and
    # per-output reliability, surfaced here instead of raising, since the
    # underlying orbit can be a genuine, physically valid ellipse even
    # when a specific downstream quantity (one that divides by (1-e) or
    # (1-e^2)) is numerically ill-conditioned at this eccentricity. Empty
    # unless at least one such warning applies; see
    # `derived.eccentricity_conditioning_warning`.
    numerical_warnings: Tuple[str, ...] = ()

    def summary(self) -> str:
        lines = [
            f"Transfer r1={self.r1_km:.3f} km -> r2={self.r2_km:.3f} km "
            f"(R={self.R:.6f}), delta_theta={self.delta_theta_deg:.3f} deg"
            + (" [Hohmann]" if self.is_hohmann else ""),
            f"  p* = {self.p_star_km:.3f} km, a = {self.semi_major_axis_km:.3f} km, "
            f"e = {self.eccentricity:.6f}",
            f"  delta_v1 = {self.delta_v1_km_s:.6f} km/s, "
            f"delta_v2 = {self.delta_v2_km_s:.6f} km/s, "
            f"total = {self.delta_v_total_km_s:.6f} km/s",
            f"  Hohmann reference total = {self.hohmann_delta_v_total_km_s:.6f} km/s "
            f"(penalty: {self.penalty_vs_hohmann_pct:+.3f} %)",
        ]
        if self.geometry_check is not None:
            lines.append(
                f"  independent geometry check: max |residual| = "
                f"{self.geometry_check.max_abs_residual:.3e} "
                f"({'OK' if self.geometry_check.ok() else 'FAILED'})"
            )
        lines.append(
            f"  time of flight = {self.time_of_flight_s / 3600.0:.3f} h, "
            f"transfer-orbit period = {self.orbital_period_s / 3600.0:.3f} h"
        )
        lines.append(
            f"  transfer orbit r_periapsis = {self.r_periapsis_km:.3f} km, "
            f"r_apoapsis = {self.r_apoapsis_km:.3f} km"
        )
        if self.altitude_info is not None:
            lines.append(
                f"  altitude: periapsis = {self.altitude_info.periapsis_altitude_km:.3f} km, "
                f"apoapsis = {self.altitude_info.apoapsis_altitude_km:.3f} km"
                + (f"  [WARNING: {self.altitude_info.warning}]" if self.altitude_info.warning else "")
            )
        if self.propellant_info_departure is not None:
            lines.append(
                f"  propellant mass fraction (departure burn, Isp="
                f"{self.propellant_info_departure.isp_s:.0f} s) = "
                f"{self.propellant_info_departure.mass_fraction * 100.0:.3f} %"
            )
        if self.propellant_info_arrival is not None:
            lines.append(
                f"  propellant mass fraction (arrival burn, Isp="
                f"{self.propellant_info_arrival.isp_s:.0f} s) = "
                f"{self.propellant_info_arrival.mass_fraction * 100.0:.3f} %"
            )
        if self.lambert_check is not None:
            lines.append(
                f"  independent Lambert-solver cross-check: max |residual| = "
                f"{self.lambert_check.max_abs_residual:.3e} "
                f"({'OK' if self.lambert_check.ok() else 'FAILED'}, "
                f"{self.lambert_check.lambert_iterations} iterations)"
            )
        if self.rk4_check is not None:
            order = self.rk4_check.observed_convergence_order
            order_str = f"{order:.2f}" if order is not None else "n/a (already at noise floor)"
            lines.append(
                f"  independent RK4 propagation cross-check: position "
                f"residual = {self.rk4_check.position_residual_km:.3e} km, "
                f"velocity residual = {self.rk4_check.velocity_residual_km_s:.3e} "
                f"km/s ({'OK' if self.rk4_check.ok() else 'FAILED'}, "
                f"{self.rk4_check.n_steps} steps, observed convergence order="
                f"{order_str})"
            )
        for warning in self.numerical_warnings:
            lines.append(f"  [NUMERICAL CONDITIONING WARNING: {warning}]")
        return "\n".join(lines)


def _normalize_delta_theta(delta_theta_deg: float) -> float:
    """Fold an input angle into (0, 360) and raise DomainError at the
    degenerate zero-transfer-angle boundary (delta_theta = 0 or 360,
    equivalently phi = +-pi)."""
    folded = delta_theta_deg % 360.0
    if folded == 0.0:
        raise DomainError(
            f"delta_theta_deg={delta_theta_deg!r} is a degenerate zero-angle "
            "transfer (delta_theta = 0 mod 360 deg); the physical domain is "
            "the open interval 0 < delta_theta < 360 deg."
        )
    return folded


def compute_transfer(
    r1_km: float,
    r2_km: float,
    delta_theta_deg: float,
    mu_km3_s2: float = EARTH_MU_KM3_S2,
    validate: bool = True,
    geometry_tolerance: float = 1e-6,
    body_radius_km: Optional[float] = None,
    isp_s: Optional[float] = None,
    lambert_check: bool = False,
    rk4_check: bool = False,
    rk4_n_steps: int = 4000,
) -> TransferResult:
    """Solve the optimal two-impulse circular-to-circular transfer.

    Parameters
    ----------
    r1_km, r2_km : departure and arrival circular-orbit radii, in km,
        measured from the central body's centre (not altitude). Requires
        r2_km > r1_km -- this component is proven only for R := r2/r1 > 1
        (an outbound transfer); for an inbound transfer (r2 < r1), swap the
        two radii to solve the equivalent outbound problem and interpret
        delta_v1/delta_v2 accordingly (this repository does not certify
        that the resulting delta-v magnitudes are identical for the
        reversed-direction transfer; no assumption of that symmetry is
        made by this function).
    delta_theta_deg : the transfer angle actually traversed (true anomaly
        swept from departure to arrival), in degrees, strictly between 0
        and 360 (mod 360). 180 deg is exactly the Hohmann transfer.
        Convention: `phi := pi - Delta_theta` (radians) -- see
        `docs/mathematical_conventions.md` at the repository root. Per the
        model's own symmetry (the cost depends only on |phi|, since the
        underlying formulas are even in phi), `delta_theta` and `360 -
        delta_theta` give the identical cost -- the "short way" and "long
        way" transfer at the same angular distance from the Hohmann
        configuration cost the same in this idealized two-body model.
    mu_km3_s2 : gravitational parameter of the central body, km^3/s^2.
        Defaults to Earth (`otca_c.units.Bodies.EARTH`); pass another value
        from `Bodies`, or your own, for a different central body.
    validate : if True (default), also run the independent geometric/
        physical cross-check (`geometry.independent_check`) and attach its
        result; set False to skip it (e.g. in a tight loop where the
        caller already trusts the solver and wants to avoid the extra
        evaluation).
    geometry_tolerance : the residual tolerance used to decide
        `geometry_check.ok()` is NOT used to raise an exception by this
        function -- the check result is always just attached to the
        returned object for the caller to inspect; this project does not
        silently hide a failed self-check, but also does not block a
        caller who wants the raw result regardless.
    body_radius_km : the central body's own physical (mean equatorial)
        radius, in km -- e.g. `otca_c.units.BodyRadii.EARTH`. Optional:
        periapsis/apoapsis *altitude* needs this second, independent
        physical constant (an orbit computation never implies it from mu
        alone), so `altitude_info` on the result stays None unless this
        is explicitly supplied here -- never guessed or defaulted.
    isp_s : the engine's specific impulse, in seconds. Optional:
        propellant mass fraction (Tsiolkovsky rocket equation) needs this
        input, which this orbital-mechanics component never has on its
        own, so `propellant_info_departure`/`_arrival` stay None unless
        this is explicitly supplied here.
    lambert_check : if True, additionally re-solve the transfer from
        scratch with an independent Lambert-problem solver
        (`lambert.solve_lambert`, the "universal variable" formulation --
        see that module's docstring) fed this transfer's own time of
        flight, and attach the comparison as `lambert_check`. This is a
        genuinely separate classical method (not a re-reading of this
        package's own S1/S2 algebra, unlike `geometry_check`), closer to
        what the project discussed as "comparison against classical/
        reference methods" -- but it is a from-scratch root-find, and
        noticeably more expensive than `geometry_check`'s closed-form
        re-derivation, so it is off by default. Near the extreme ends of
        the angular domain (transfer angles within a few degrees of 0 or
        360 deg), this specific numerical method needs more iterations
        and its residual, while always converging in this package's
        tests, is correspondingly larger (still many orders of magnitude
        below any physically meaningful threshold) -- a property of this
        verification method's own conditioning near that geometry, not an
        indication of reduced confidence in `compute_transfer`'s own
        result, whose validity is established independently by this
        repository's existence/uniqueness result and Theorems 1-2 (see
        `STATUS.md`) over the full open domain 0 < delta_theta < 360 deg,
        R > 1 (no narrower certified subrange exists for the core
        cost-minimization solve itself).
    rk4_check : if True, additionally validate this transfer by directly,
        numerically integrating the two-body equations of motion
        (classical fixed-step RK4, `rk4.rk4_validate`) from this
        transfer's own departure state vector for its own time of flight,
        and compare the result to the already-known arrival state. This
        is a THIRD, independent method, distinct from both
        `geometry_check` (re-reads the same solved orbit through
        different closed-form textbook relations, no time integration at
        all) and `lambert_check` (a from-scratch boundary-value
        root-find): an initial-value numerical integrator that exercises
        `time_of_flight_s` against the actual dynamics. It is the most
        expensive of the three cross-checks (two full RK4 propagations,
        for a built-in step-doubling convergence study, on every call),
        so it is off by default. `rk4_check` is None if `rk4_check=True`
        raised `ConvergenceError` internally (e.g. the convergence study
        did not show the expected ~4th order for this geometry within
        `rk4_n_steps`) -- never a silently wrong answer.
    rk4_n_steps : the (coarser of the two) step counts used by the RK4
        convergence study; the finer run uses `2*rk4_n_steps`. The
        default is adequate for the transfer geometries exercised by this
        package's own test suite; a very long time of flight or a very
        high eccentricity may need more steps to reach the convergence
        study's expected ~4th-order regime (see `rk4.py`).

    Raises
    ------
    DomainError : if the inputs fall outside the domain this component's
        underlying theorems (`circular/existence_uniqueness_minimizer/`,
        `circular/global_monotonicity/`) are proven on.
    ConvergenceError : if the root-finder unexpectedly fails to locate p*
        despite the input being inside the proven domain (see `core.py`).
    """
    if not (r1_km > 0 and r2_km > 0):
        raise DomainError(f"radii must be positive: r1_km={r1_km!r}, r2_km={r2_km!r}")
    if not r2_km > r1_km:
        raise DomainError(
            f"this component requires r2_km > r1_km (R:=r2/r1>1 is the proven "
            f"domain of this package's own correctness results, Theorems 1-2 "
            f"and the existence/uniqueness result; see STATUS.md); got "
            f"r1_km={r1_km!r}, "
            f"r2_km={r2_km!r}. Swap r1_km and r2_km to solve the equivalent "
            "outbound transfer if you intended an inbound leg."
        )
    R = r2_km / r1_km

    delta_theta_deg_folded = _normalize_delta_theta(delta_theta_deg)
    delta_theta_rad = math.radians(delta_theta_deg_folded)
    phi = math.pi - delta_theta_rad
    phi_abs = abs(phi)  # model is even in phi; see docstring above.

    is_hohmann = math.isclose(phi_abs, 0.0, abs_tol=1e-12)

    if is_hohmann:
        h = hohmann.solve_hohmann(R)
        p_star, e, a = h.p, h.e, h.a
        S1, S2 = h.delta_v1**2, h.delta_v2**2
        geom_check = None
        if validate:
            # The Hohmann closed form and the general conic/vis-viva check
            # coincide at phi=0 (radial velocity is exactly zero at both
            # apsides); run it anyway at a tiny but nonzero phi surrogate
            # is unnecessary complexity for a case already exact in closed
            # form -- skip, and say so explicitly rather than fake a check.
            geom_check = None
        # Exact Hohmann geometry: departure at periapsis (nu=0), arrival
        # at apoapsis (nu=pi) -- not an approximation, the exact closed
        # form for this special case.
        nu1, nu2 = 0.0, math.pi
    else:
        if not (0.0 < phi_abs < math.pi):
            raise DomainError(
                f"phi={phi!r} rad (from delta_theta_deg={delta_theta_deg!r}) is "
                "outside the open physical domain (0, pi); this occurs only at "
                "the degenerate boundary already excluded by delta_theta "
                "normalization above -- please report this as a bug if reached."
            )
        sol = core.solve(phi_abs, R)
        p_star, e, a = sol.p_star, sol.e, sol.a
        S1, S2 = sol.S1, sol.S2
        # IMPORTANT: the cost/geometry solve above (p*, e, a, S1, S2) is
        # correctly computed from phi_abs, since the underlying cost
        # formula is even in phi (depends only on cos(phi) and sin(phi)^2
        # -- an elementary algebraic property of the S1/S2 primitives
        # themselves, core.py's _e2/_S1/_S2, not a conclusion of Theorem 2
        # or any other numbered theorem). But the *arc actually swept* --
        # and hence the true anomalies nu1/nu2 and the time of flight --
        # is NOT symmetric: delta_theta and 360-delta_theta are the two
        # complementary arcs (short way / long way) of the very same
        # ellipse, not the same arc. Using phi_abs here would silently
        # collapse both to the same (short-way) arc and report the wrong
        # time of flight whenever delta_theta_deg > 180 -- so the signed
        # `phi` (not phi_abs) is used from this point on, exactly as
        # `geometry.independent_check` itself already expects (it derives
        # delta_theta = pi - phi internally, needing phi's sign to get
        # the real swept angle rather than always the <=180 one).
        geom_check = (
            independent_check(p_star, phi, R, e, S1, S2) if validate else None
        )
        c, s = math.cos(phi), math.sin(phi)
        x = p_star - 1.0
        y = p_star / R - 1.0
        e_sin_nu1 = -(x * c + y) / s
        nu1 = math.atan2(e_sin_nu1, x)
        nu2 = nu1 + (math.pi - phi)  # == nu1 + delta_theta_rad, the real sweep

    vc1 = circular_speed_km_s(r1_km, mu_km3_s2)
    delta_v1_km_s = (S1 ** 0.5) * vc1
    delta_v2_km_s = (S2 ** 0.5) * vc1
    delta_v_total_km_s = delta_v1_km_s + delta_v2_km_s

    hohmann_ref = hohmann.solve_hohmann(R)
    hohmann_total_km_s = hohmann_ref.delta_v_total * vc1
    penalty_pct = 100.0 * (delta_v_total_km_s / hohmann_total_km_s - 1.0)

    # --- Unconditional derived quantities (no extra input needed) -------
    p_star_km = p_star * r1_km
    a_km = a * r1_km
    r_peri_km = derived.periapsis_radius_km(p_star_km, e)
    r_apo_km = derived.apoapsis_radius_km(p_star_km, e)
    period_s = derived.orbital_period_s(a_km, mu_km3_s2)
    tof_s = derived.time_of_flight_s(a_km, e, nu1, nu2, mu_km3_s2)

    # --- Optional derived quantities (only if the extra input was given) -
    alt_info = (
        derived.altitude_info(r_peri_km, r_apo_km, body_radius_km)
        if body_radius_km is not None
        else None
    )
    prop_departure = (
        derived.propellant_info(delta_v1_km_s, isp_s) if isp_s is not None else None
    )
    prop_arrival = (
        derived.propellant_info(delta_v2_km_s, isp_s) if isp_s is not None else None
    )

    # --- Per-output numerical-conditioning warnings (not a domain error:
    # e is validated by this point, this only flags that a FEW downstream
    # quantities built from 1/(1-e) or 1/(1-e^2) can carry a much larger
    # relative floating-point error than the rest of the output, at this
    # specific eccentricity). See derived.eccentricity_conditioning_warning
    # for exactly which quantities are affected and why.
    numerical_warnings = []
    e_cond_warning = derived.eccentricity_conditioning_warning(e)
    if e_cond_warning is not None:
        numerical_warnings.append(e_cond_warning)
    numerical_warnings = tuple(numerical_warnings)

    lambert_result = None
    if lambert_check:
        if is_hohmann:
            # Exactly dnu=180 deg is a genuine degeneracy of THIS
            # solver's own A-formula (A = sqrt(2 r1 r2) cos(dnu/2) -> 0),
            # not of the underlying geometry -- lambert.solve_lambert
            # correctly raises ConvergenceError here rather than return a
            # numerically unstable answer. This package already has an
            # exact closed form for Hohmann (hohmann.py) and never needed
            # Lambert to confirm that case; skip it here the same way
            # geometry_check is skipped for Hohmann above, rather than
            # letting the exception propagate out of a successful call.
            lambert_result = None
        else:
            try:
                lambert_result = derived.lambert_cross_check(
                    r1_km=r1_km,
                    r2_km=r2_km,
                    delta_theta_rad=delta_theta_rad,
                    tof_s=tof_s,
                    mu_km3_s2=mu_km3_s2,
                    p_star_km=p_star_km,
                    e=e,
                    delta_v1_km_s=delta_v1_km_s,
                    delta_v2_km_s=delta_v2_km_s,
                )
            except ConvergenceError:
                # The Lambert cross-check is a secondary, optional
                # verification tool -- an EXPECTED numerical failure in
                # it (no valid bracket for some extreme geometry, a
                # degenerate angle, an internal consistency check
                # tripping) must never turn a successful primary solve
                # into a hard error for the caller. None signals "not
                # available for this case", exactly like
                # geometry_check=None for Hohmann. Only ConvergenceError
                # is caught here, deliberately NOT a bare `except
                # Exception` (flagged during independent review as able
                # to silently swallow a real programming bug) -- any
                # other exception (a TypeError, an AttributeError, ...)
                # is left to propagate, since those indicate an actual
                # defect that must not be hidden behind a clean-looking
                # lambert_check=None.
                lambert_result = None

    rk4_result = None
    if rk4_check:
        # Unlike lambert_check, the exact Hohmann transfer (nu1=0,
        # nu2=pi) is not a degeneracy for RK4 -- it is simply propagated
        # like any other case, no special-casing needed.
        try:
            rk4_result = rk4_validate(
                p_star_km=p_star_km,
                e=e,
                nu1_rad=nu1,
                delta_theta_rad=delta_theta_rad,
                tof_s=tof_s,
                mu_km3_s2=mu_km3_s2,
                r2_km=r2_km,
                n_steps=rk4_n_steps,
            )
        except ConvergenceError:
            # Same philosophy as lambert_check: an expected numerical
            # failure in this optional, secondary cross-check (the
            # convergence study not showing the expected order for this
            # specific geometry/step count) must never turn a successful
            # primary solve into a hard error. None signals "not
            # available for this case". Only ConvergenceError is caught,
            # deliberately not a bare `except Exception`, for the same
            # reason given above for lambert_check.
            rk4_result = None

    return TransferResult(
        r1_km=r1_km,
        r2_km=r2_km,
        delta_theta_deg=delta_theta_deg,
        mu_km3_s2=mu_km3_s2,
        R=R,
        phi_rad=phi,
        is_hohmann=is_hohmann,
        p_star_km=p_star_km,
        semi_major_axis_km=a_km,
        eccentricity=e,
        delta_v1_km_s=delta_v1_km_s,
        delta_v2_km_s=delta_v2_km_s,
        delta_v_total_km_s=delta_v_total_km_s,
        hohmann_delta_v_total_km_s=hohmann_total_km_s,
        penalty_vs_hohmann_pct=penalty_pct,
        geometry_check=geom_check,
        nu1_rad=nu1,
        nu2_rad=nu2,
        r_periapsis_km=r_peri_km,
        r_apoapsis_km=r_apo_km,
        orbital_period_s=period_s,
        time_of_flight_s=tof_s,
        altitude_info=alt_info,
        propellant_info_departure=prop_departure,
        propellant_info_arrival=prop_arrival,
        lambert_check=lambert_result,
        rk4_check=rk4_result,
        numerical_warnings=numerical_warnings,
    )
