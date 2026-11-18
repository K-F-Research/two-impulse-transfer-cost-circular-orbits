"""
Theorem 2 -- Global Monotonicity of the Optimized Transfer Cost
Independent Verification Script

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Independent re-implementation, written directly from the definitions in
`theorem2_global_monotonicity_proof.md` (Sec. 1), that performs two checks:

  1. Symbolic re-derivation (SymPy) of the proof's two central closed-form
     identities -- D'(p) = 2(1+c)/(p^3(1-c)) (Sec. 6, Lemma 6.1) and the
     factorization g(p) = A(p-p1*)(p-p2*) (Sec. 8, Lemma 8.3) -- confirmed
     against the definitions with zero symbolic residual.
  2. Direct numerical confirmation of the monotonicity claim itself: for
     several radius ratios R, V_R(phi) := min_p F(p,phi,R) is evaluated at a
     grid of angles and checked to be strictly increasing.

Both checks are confirmatory only, exactly as the methodology stated
throughout the proof document requires: they do not replace the analytic
proof (Sec. 1-13 of the companion document) and are not used to derive any
of its conclusions.

Status: EXPLORATORY / confirmatory. Reproducibility: REPRODUCED -- both
checks below have been executed with the results stated in this docstring:
zero symbolic residual, zero monotonicity violations over 472
consecutive-pair comparisons (8 ratios R in [1.2, 50], 60 angles each from
0.5 deg to 179.5 deg).
"""

import numpy as np
import sympy as sp
from scipy.optimize import minimize_scalar


def check_central_identities():
    """Symbolically re-derive the two central closed-form identities of the
    proof (Lemma 6.1 and Lemma 8.3) directly from the Sec. 1 definitions,
    and confirm them with zero residual."""
    p, c, R = sp.symbols('p c R', positive=True)
    x = p - 1
    y = p / R - 1
    s2 = 1 - c**2
    e2 = x**2 + (x * c + y)**2 / s2

    # D(p) := [d/dp(e^2)*p - (e^2-1)] / p^2  (Sec. 6 definition)
    dE2 = sp.diff(e2, p)
    D = (dE2 * p - (e2 - 1)) / p**2
    Dp = sp.simplify(sp.diff(D, p))
    claimed_Dp = 2 * (1 + c) / (p**3 * (1 - c))
    residual_D = sp.simplify(Dp - claimed_Dp)

    # g(p) := (x*c+y)(x+y*c)  (Sec. 8 definition), claimed factorization
    g = (x * c + y) * (x + y * c)
    p1star = R * (1 + c) / (1 + R * c)
    p2star = R * (1 + c) / (R + c)
    A = (R + c) * (1 + R * c) / R**2
    factored = sp.expand(A * (p - p1star) * (p - p2star))
    residual_g = sp.simplify(sp.expand(g) - factored)

    return residual_D, residual_g


def transfer_cost(p, phi, R):
    """F(p, phi, R): total two-impulse transfer cost, r1=1, r2=R (Sec. 1)."""
    c = np.cos(phi)
    x = p - 1.0
    y = p / R - 1.0
    e2 = x**2 + (x * c + y)**2 / (1.0 - c**2)
    v1_sq = 2.0 + (e2 - 1.0) / p
    v2_sq = 2.0 / R + (e2 - 1.0) / p
    v_theta1 = np.sqrt(p)
    v_theta2 = np.sqrt(p) / R
    v_c1 = 1.0
    v_c2 = 1.0 / np.sqrt(R)
    S1 = v1_sq + v_c1**2 - 2 * v_c1 * v_theta1
    S2 = v2_sq + v_c2**2 - 2 * v_c2 * v_theta2
    return np.sqrt(max(S1, 1e-300)) + np.sqrt(max(S2, 1e-300))


def optimized_cost(phi, R):
    """V_R(phi) := min_p F(p, phi, R)."""
    result = minimize_scalar(
        lambda p: transfer_cost(p, phi, R),
        bounds=(1e-6, 1e6), method='bounded',
        options={'xatol': 1e-12},
    )
    return result.fun


def check_monotonicity(
    ratios=(1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 25.0, 50.0),
    angles_deg=np.linspace(0.5, 179.5, 60),
):
    """Confirms V_R(phi_1) < V_R(phi_2) for all consecutive sampled pairs
    phi_1 < phi_2, across several ratios R. Returns the number of pairs
    checked and the number of violations found."""
    checked = 0
    violations = []
    for R in ratios:
        values = [optimized_cost(np.radians(a), R) for a in angles_deg]
        for i in range(len(values) - 1):
            checked += 1
            if not (values[i] < values[i + 1] - 1e-9):
                violations.append((R, angles_deg[i], angles_deg[i + 1]))
    return checked, violations


if __name__ == "__main__":
    print("Theorem 2 -- independent verification")
    print("=" * 60)

    print("\n[1] Re-deriving the two central closed-form identities...")
    residual_D, residual_g = check_central_identities()
    print(f"    D'(p) residual (expect 0): {residual_D}")
    print(f"    g(p) factorization residual (expect 0): {residual_g}")
    assert residual_D == 0, "D'(p) identity FAILED to reproduce"
    assert residual_g == 0, "g(p) factorization FAILED to reproduce"
    print("    Both identities confirmed with zero residual.")

    print("\n[2] Independently checking global monotonicity numerically...")
    checked, violations = check_monotonicity()
    print(f"    Consecutive-pair comparisons checked: {checked}")
    print(f"    Violations found: {len(violations)}")
    if violations:
        for v in violations:
            print(f"      VIOLATION at R={v[0]}, phi1={v[1]} deg, phi2={v[2]} deg")
    assert not violations, "Monotonicity violation found"
    print("    No violations: monotonicity confirmed on this grid.")

    print("\nAll checks passed.")
