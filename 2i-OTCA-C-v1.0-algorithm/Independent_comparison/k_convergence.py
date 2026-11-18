"""
2i-OTCA-C v1.0 -- Independent numerical cross-check (Lambert p-iteration method)
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See the project LICENSE file / ATBD Annex A for the full license text.

Independent numerical verification of the curvature coefficient K(r1,r2) of
Theorem 2, by finite differences of the EXACT solver's own cost function
(otca_c.api.compute_transfer), evaluated around the Hohmann point (phi=0).

Simple quotient:
  K_num(phi) = [dV_exact(phi) - dV_exact(0)] / phi^2

Symmetric (centered) quotient, cancels odd-order terms in the expansion:
  K_sym(h) = [dV_exact(+h) + dV_exact(-h) - 2*dV_exact(0)] / (2*h^2)

Both are compared to the closed-form K(1,R)*sqrt(mu/r1) of Theorem 2
(otca_c.surrogate._K_normalized), across a convergence sweep of phi,
to find the regime where O(phi^2) truncation error dominates vs. where
the exact solver's own numerical error dominates (expected to degrade
again as phi -> 0 because of subtractive cancellation of nearly-equal
costs).

Everything here is evaluated at 50-digit precision throughout: the EXACT
solver itself is double precision (otca_c.core), so this script cannot
remove that floor, but using mpmath for the quotient arithmetic and
differencing at least isolates the solver's own double-precision floor
from any additional error introduced by the finite-difference bookkeeping.
"""
import os
import sys
# Adjust this path to the otca_c package's src/ directory if run standalone.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import mpmath as mp
mp.mp.dps = 50
from otca_c import api, surrogate

mu = 398600.4418

def dv_exact(r1, r2, phi_rad):
    """dV_exact(phi), phi signed deviation from Hohmann in radians."""
    delta_theta_deg = 180.0 - mp.degrees(phi_rad)
    res = api.compute_transfer(r1, r2, float(delta_theta_deg), mu_km3_s2=mu, validate=False)
    return mp.mpf(res.delta_v_total_km_s)

def K_theoretical(r1, r2):
    Rn = r2 / r1
    K_norm = surrogate._K_normalized(Rn)  # dimensionless K(1,R)
    vc1 = (mu / r1) ** 0.5
    return mp.mpf(K_norm) * mp.mpf(vc1)   # physical K(r1,r2) in km/s

def convergence_study(r1, r2, phis_rad):
    K_th = K_theoretical(r1, r2)
    dv0 = dv_exact(r1, r2, mp.mpf(0))
    print(f"  K_theoretical = {float(K_th):.12f} km/s   (dv0 = {float(dv0):.10f} km/s)")
    print(f"  {'phi(rad)':>12} {'K_num(simple)':>16} {'rel_err_simple':>14} {'K_sym(centered)':>18} {'rel_err_sym':>12}")
    rows = []
    for phi in phis_rad:
        phi_mp = mp.mpf(phi)
        dvp = dv_exact(r1, r2, phi_mp)
        dvm = dv_exact(r1, r2, -phi_mp)
        K_simple = (dvp - dv0) / phi_mp**2
        K_sym = (dvp + dvm - 2*dv0) / (2*phi_mp**2)
        rel_simple = abs(K_simple - K_th) / K_th
        rel_sym = abs(K_sym - K_th) / K_th
        rows.append((float(phi), float(K_simple), float(rel_simple), float(K_sym), float(rel_sym)))
        print(f"  {float(phi):>12.6f} {float(K_simple):>16.8f} {float(rel_simple):>14.3e} {float(K_sym):>18.8f} {float(rel_sym):>12.3e}")
    return rows

if __name__ == "__main__":
    import json
    phis = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001]
    results = {}
    for R in [1.2, 1.5, 2.0, 5.0]:
        r1 = 6678.0
        r2 = r1 * R
        print(f"\n=== R={R} (r1={r1}, r2={r2}) ===")
        rows = convergence_study(r1, r2, phis)
        results[str(R)] = rows
    out_dir = os.path.dirname(os.path.abspath(__file__))
    json.dump(results, open(os.path.join(out_dir, 'k_convergence_results.json'), 'w'), indent=2)
