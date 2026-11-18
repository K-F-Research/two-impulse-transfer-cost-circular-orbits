"""
2i-OTCA-C v1.0 -- Independent numerical cross-check (Lambert p-iteration method)
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See the project LICENSE file / ATBD Annex A for the full license text.

Control C: physical scaling properties, tested against the INDEPENDENT
script itself (not otca_c), so these checks do not depend on core.py at all.
  - Length scaling: r1,r2 -> lambda*r1, lambda*r2 (mu fixed) should scale
    dV by 1/sqrt(lambda).
  - Gravitational-parameter scaling: mu -> q*mu (geometry fixed) should
    scale dV by sqrt(q).
"""
from lambert_independent_mp import min_cost_transfer_mp

mu0 = 398600.4418
r1, r2, theta = 6678.0, 10017.0, 170.0

_, dv0, _ = min_cost_transfer_mp(r1, r2, theta, mu0)

print("Length scaling (dV should scale as 1/sqrt(lambda)):")
for lam in [0.25, 0.5, 2.0, 4.0, 10.0]:
    _, dv_s, _ = min_cost_transfer_mp(r1*lam, r2*lam, theta, mu0)
    predicted = dv0 / lam**0.5
    rel = abs(dv_s - predicted) / predicted
    print(f"  lambda={lam:<6} dV_scaled={dv_s:.10f}  predicted={predicted:.10f}  rel_diff={rel:.3e}")

print("\nGravitational-parameter scaling (dV should scale as sqrt(q)):")
for q in [0.25, 0.5, 2.0, 4.0, 10.0]:
    _, dv_s, _ = min_cost_transfer_mp(r1, r2, theta, mu0*q)
    predicted = dv0 * q**0.5
    rel = abs(dv_s - predicted) / predicted
    print(f"  q={q:<6} dV_scaled={dv_s:.10f}  predicted={predicted:.10f}  rel_diff={rel:.3e}")
