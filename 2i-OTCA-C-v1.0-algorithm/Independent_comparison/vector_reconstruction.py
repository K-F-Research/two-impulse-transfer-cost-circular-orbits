"""
2i-OTCA-C v1.0 -- Independent numerical cross-check (Lambert p-iteration method)
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See the project LICENSE file / ATBD Annex A for the full license text.

Control A: reconstruct the two impulses as explicit 2D vectors (radial +
tangential components in the local orbital-plane frame at each point),
rather than trusting the pre-simplified law-of-cosines scalar formula
(Shen & Tsiotras Eq. 4) used in lambert_independent_mp.py.

This is an algebraic cross-check of that formula, not a new physical
method -- vr^2 + (vt-vc)^2 is algebraically identical to
v^2+vc^2-2*v*vc*cos(gamma) when cos(gamma)=vt/v. Running both paths
independently in code catches any transcription/coding error in deriving
or typing the scalar formula.
"""
import mpmath as mp
mp.mp.dps = 50
from lambert_independent_mp import _e_of_p, _bracket_admissible_interval, _golden_section_min

def dv_vector_form(p, r1, r2, theta, mu):
    ct, st = mp.cos(theta), mp.sin(theta)
    ecf1 = p / r1 - 1
    esf1 = ((p / r1 - 1) * ct - (p / r2 - 1)) / st
    e = mp.sqrt(ecf1**2 + esf1**2)
    if e >= 1:
        return mp.nan
    f1 = mp.atan2(esf1, ecf1)
    f2 = f1 + theta
    h = mp.sqrt(mu * p)

    def dv_at(r, f):
        vr = (mu / h) * e * mp.sin(f)   # radial component, transfer orbit
        vt = h / r                       # tangential component, transfer orbit
        vc = mp.sqrt(mu / r)             # circular speed (purely tangential, vr=0)
        # vector subtraction in the (radial, tangential) local frame:
        dvr = vr - 0
        dvt = vt - vc
        return mp.sqrt(dvr**2 + dvt**2)

    return dv_at(r1, f1) + dv_at(r2, f2)

if __name__ == "__main__":
    from lambert_independent_mp import transfer_cost_mp, min_cost_transfer_mp
    mu = 398600.4418
    cases = [(6678.0, 6678.0*1.2, 175.0), (6678.0, 6678.0*1.5, 170.0),
             (6678.0, 6678.0*5.0, 155.0), (6678.0, 6678.0*1.02, 178.0)]
    print("Cross-check: law-of-cosines formula vs explicit vector-component reconstruction")
    for r1, r2, theta_deg in cases:
        p_opt, dv_scalar, bracket = min_cost_transfer_mp(r1, r2, theta_deg, mu)
        theta = mp.radians(mp.mpf(theta_deg))
        dv_vec = float(dv_vector_form(mp.mpf(p_opt), mp.mpf(r1), mp.mpf(r2), theta, mp.mpf(mu)))
        rel = abs(dv_vec - dv_scalar) / abs(dv_scalar)
        print(f"  r1={r1:.1f} r2={r2:.1f} theta={theta_deg}:  scalar={dv_scalar:.14f}  vector={dv_vec:.14f}  rel_diff={rel:.3e}")
