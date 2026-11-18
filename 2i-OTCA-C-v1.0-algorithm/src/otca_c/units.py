"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: none (physical constants and unit conversions
only; standard astronomical data, not a result proven in this
repository).

Physical units and gravitational-parameter constants.

The underlying mathematical model (see `core.py` and the proof documents in
`circular/` at the repository root) is stated in normalized units: distances
in units of r1 (the departure, inner circular orbit's radius), speeds in
units of the departure orbit's own circular speed `vc1 = sqrt(mu/r1)`.

This module converts between those normalized units and real physical units
(kilometres, seconds), so that a caller can pass real orbital radii and read
back a real delta-v directly in km/s.

Standard gravitational parameters (mu = G*M) below are IAU/JPL values,
reproduced here as plain constants (not re-derived) -- they are standard
astronomical data, not a result proven elsewhere in this repository.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Bodies:
    """Standard gravitational parameters mu = G*M, in km^3/s^2.

    Source: JPL Planetary and Lunar Ephemerides / IAU current best estimates
    (values accurate to the digits shown; sufficient for the two-impulse
    cost model here, which does not itself claim higher precision than its
    own idealizations -- point-mass, two-body, coplanar, impulsive burns).
    """

    SUN: float = 1.32712440018e11
    MERCURY: float = 2.2032e4
    VENUS: float = 3.24859e5
    EARTH: float = 3.986004418e5
    MOON: float = 4.9048695e3
    MARS: float = 4.282837e4
    JUPITER: float = 1.26686534e8
    SATURN: float = 3.7931187e7
    URANUS: float = 5.793939e6
    NEPTUNE: float = 6.836529e6


@dataclass(frozen=True)
class BodyRadii:
    """Mean equatorial radii, in km (IAU current best estimates), one entry
    per body listed in `Bodies` above so the two stay consistent. Kept
    separate from `Bodies` (mu) because a periapsis/apoapsis *altitude*
    needs the body's physical radius, which an orbit computation never
    implies on its own -- it is a second, independent physical constant,
    not derived from mu. `compute_transfer` only reports an altitude when
    this value is explicitly supplied (here or as a custom float of the
    caller's own choosing -- any body not listed here, e.g. a specific
    asteroid or space station, can simply pass its own radius in km),
    rather than guessing one."""

    SUN: float = 695700.0
    MERCURY: float = 2439.7
    VENUS: float = 6051.8
    EARTH: float = 6378.137
    MOON: float = 1737.4
    MARS: float = 3389.5
    JUPITER: float = 69911.0
    SATURN: float = 58232.0
    URANUS: float = 25362.0
    NEPTUNE: float = 24622.0


# Convenience plain constant, used throughout this package's defaults.
EARTH_MU_KM3_S2 = Bodies.EARTH


def circular_speed_km_s(r_km: float, mu_km3_s2: float) -> float:
    """v_c = sqrt(mu/r), the circular orbital speed at radius r_km."""
    if r_km <= 0:
        raise ValueError(f"radius must be positive, got {r_km!r} km")
    if mu_km3_s2 <= 0:
        raise ValueError(f"mu must be positive, got {mu_km3_s2!r} km^3/s^2")
    return (mu_km3_s2 / r_km) ** 0.5


def to_normalized_speed(delta_v_km_s: float, r1_km: float, mu_km3_s2: float) -> float:
    """Convert a real delta-v (km/s) into the model's normalized speed unit
    (multiples of vc1 = sqrt(mu/r1))."""
    return delta_v_km_s / circular_speed_km_s(r1_km, mu_km3_s2)


def to_physical_speed_km_s(delta_v_normalized: float, r1_km: float, mu_km3_s2: float) -> float:
    """Convert a normalized speed (multiples of vc1) into a real delta-v
    (km/s)."""
    return delta_v_normalized * circular_speed_km_s(r1_km, mu_km3_s2)
