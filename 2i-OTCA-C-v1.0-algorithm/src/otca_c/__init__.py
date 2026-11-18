"""
2i-OTCA-C -- Two-Impulse Orbital Transfer Cost Algorithm -- Circular Case

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Theorems relied upon: none directly (this module only re-exports the
public names defined elsewhere; see each module's own docstring, and
`docs_atbd/module_theorem_correspondence.md`, for the exact citations).

Public API:

    from otca_c import compute_transfer, Bodies

    result = compute_transfer(
        r1_km=6778.0,           # 400 km altitude LEO (Earth radius 6378 km)
        r2_km=13556.0,          # R = 2 case
        delta_theta_deg=170.0,  # transfer angle actually traversed
        mu_km3_s2=Bodies.EARTH,
    )
    print(result.delta_v_total_km_s, result.p_star_km)

For a fast, closed-form estimate near the Hohmann configuration (Theorem
2's quadratic surrogate -- cheap, but only validated within 30 deg of
Hohmann; see `otca_c.surrogate`'s module docstring for the honest limits):

    from otca_c import estimate_cost, compare_to_exact

    fast = estimate_cost(r1_km=6678.0, r2_km=10017.0, delta_theta_deg=170.0)
    print(fast.delta_v_total_estimate_km_s, fast.within_measured_domain)

To find the largest angular deviation a given delta-v budget affords
(Theorem 1 for existence/uniqueness, Theorem 3's convexity for the
solver's own structure; see `otca_c.budget_inversion`'s module docstring):

    from otca_c import invert_phi_for_budget

    budget_result = invert_phi_for_budget(
        r1_km=6678.0, r2_km=13356.0, delta_v_budget_km_s=4.5,
    )
    print(budget_result.delta_theta_min_deg, budget_result.phi_max_rad)

To classify a R=2 transfer relative to Theorem 4's certified curvature
sign-change point (R=2 only; see `otca_c.curvature_regime`'s module
docstring):

    from otca_c import classify_curvature_regime

    regime = classify_curvature_regime(
        r1_km=6678.0, r2_km=13356.0, delta_theta_deg=20.0,
    )
    print(regime.regime, regime.distance_to_wc)

See `README.md` in this folder for the full mathematical contract (domain
of validity, angle convention, units, and references to the theorems this
implementation relies on), and `docs/mathematical_conventions.md` at the
repository root for the shared notation.
"""
from .api import compute_transfer, TransferResult
from .derived import AltitudeInfo, LambertCrossCheck, PropellantInfo
from .rk4 import RK4Validation
from .units import Bodies, BodyRadii
from .errors import DomainError, ConvergenceError
from .surrogate import (
    MAX_VALIDATED_PHI_DEG,
    SurrogateComparison,
    SurrogateEstimate,
    compare_to_exact,
    curvature_coefficient_km_s,
    estimate_cost,
)
from .budget_inversion import BudgetInversionResult, invert_phi_for_budget
from .curvature_regime import (
    PHI_C_CERTIFIED_DEG,
    W_C_CERTIFIED,
    CurvatureRegime,
    classify_curvature_regime,
)

__all__ = [
    "compute_transfer",
    "TransferResult",
    "Bodies",
    "BodyRadii",
    "AltitudeInfo",
    "PropellantInfo",
    "LambertCrossCheck",
    "RK4Validation",
    "DomainError",
    "ConvergenceError",
    "estimate_cost",
    "compare_to_exact",
    "curvature_coefficient_km_s",
    "SurrogateEstimate",
    "SurrogateComparison",
    "MAX_VALIDATED_PHI_DEG",
    "invert_phi_for_budget",
    "BudgetInversionResult",
    "classify_curvature_regime",
    "CurvatureRegime",
    "W_C_CERTIFIED",
    "PHI_C_CERTIFIED_DEG",
]

__version__ = "0.1.0"
