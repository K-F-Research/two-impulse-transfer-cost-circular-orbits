"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: Theorem 2 of the manuscript (strict local
minimum of the cost at the Hohmann transfer, phi=0 -- see below for
the known numbering note on this result's proof folder).

Exact closed-form Hohmann transfer (the phi=0 boundary case).

`core.solve` evaluates the general-phi primitive formulas, which contain a
`1/sin(phi)` term (see `core.py`'s `_e2`) and are therefore not directly
evaluable at exactly `phi=0`. The Hohmann transfer is the single classical,
closed-form special case this component's general solver approaches as a
limit (Theorem 2 of the manuscript -- proof content filed internally
under `circular/hohmann_local_minimum/`, labeled there as "Theorem 1";
see STATUS.md's numbering note -- proves it is in fact a
strict local minimum of the cost); it is reproduced here directly from the
standard textbook formula (not re-derived from this repository's own
algebra), both to serve requests at exactly phi=0 and as an independent
reference point to validate the general solver's phi->0 limit against
(see `tests/test_hohmann_limit.py`).
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class HohmannResult:
    """Normalized units (r1=1, mu=1), consistent with core.py."""

    p: float
    e: float
    a: float
    delta_v1: float
    delta_v2: float
    delta_v_total: float


def solve_hohmann(R: float) -> HohmannResult:
    """Exact classical Hohmann transfer between circular orbits of radius
    1 and R (R>1), normalized units (r1=1, mu=1)."""
    if not R > 1:
        raise ValueError(f"Hohmann transfer requires R>1, got R={R!r}")
    a = (1.0 + R) / 2.0
    e = (R - 1.0) / (R + 1.0)
    p = a * (1 - e**2)
    vc1 = 1.0
    vc2 = 1.0 / (R ** 0.5)
    v_peri = (2.0 / 1.0 - 1.0 / a) ** 0.5  # vis-viva at periapsis r=1
    v_apo = (2.0 / R - 1.0 / a) ** 0.5  # vis-viva at apoapsis r=R
    delta_v1 = abs(v_peri - vc1)
    delta_v2 = abs(vc2 - v_apo)
    return HohmannResult(
        p=p, e=e, a=a,
        delta_v1=delta_v1, delta_v2=delta_v2,
        delta_v_total=delta_v1 + delta_v2,
    )
