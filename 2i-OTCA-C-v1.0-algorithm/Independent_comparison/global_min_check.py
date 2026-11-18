"""
2i-OTCA-C v1.0 -- Independent numerical cross-check (Lambert p-iteration method)
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See the project LICENSE file / ATBD Annex A for the full license text.

Control B: verify the golden-section result is the GLOBAL minimum of dV(p)
on the admissible bracket, not a local minimum or boundary artifact.
Does this by dense sampling across the full admissible (p_lo, p_hi) interval
and confirming the sampled minimum matches the golden-section result, and
by checking the curve shape (single interior minimum, no other stationary
points) for each representative case.
"""
import mpmath as mp
mp.mp.dps = 50
from lambert_independent_mp import transfer_cost_mp, min_cost_transfer_mp, _bracket_admissible_interval

def scan_and_verify(r1, r2, theta_deg, mu, n=4000):
    theta = mp.radians(mp.mpf(theta_deg))
    ct, st = mp.cos(theta), mp.sin(theta)
    r1m, r2m, mum = mp.mpf(r1), mp.mpf(r2), mp.mpf(mu)
    d = mp.sqrt(r1m**2 + r2m**2 - 2*r1m*r2m*ct)
    s = (r1m + r2m + d) / 2
    p_lo, p_hi = _bracket_admissible_interval(r1m, r2m, ct, st, s)
    eps = (p_hi - p_lo) * mp.mpf('1e-6')
    ps = [p_lo + eps + (p_hi - p_lo - 2*eps) * k / (n-1) for k in range(n)]
    vals = [float(transfer_cost_mp(p, r1m, r2m, theta, mum)) for p in ps]
    i_min = min(range(n), key=lambda i: vals[i])
    p_scan_min, dv_scan_min = float(ps[i_min]), vals[i_min]

    p_opt, dv_opt, _ = min_cost_transfer_mp(r1, r2, theta_deg, mu)

    # count local minima (sign changes in discrete derivative) to confirm unimodality
    n_local_min = 0
    for i in range(1, n-1):
        if vals[i] < vals[i-1] and vals[i] < vals[i+1]:
            n_local_min += 1

    rel_gap = abs(dv_opt - dv_scan_min) / abs(dv_scan_min)
    near_lo = (p_scan_min - p_lo) < (p_hi - p_lo) * 0.02
    near_hi = (p_hi - p_scan_min) < (p_hi - p_lo) * 0.02
    print(f"  r1={r1:.2f} r2={r2:.2f} theta={theta_deg}: bracket=({float(p_lo):.3f},{float(p_hi):.3f}) "
          f"n_local_minima={n_local_min} scan_min_p={p_scan_min:.4f} golden_p={p_opt:.4f} "
          f"rel_gap_golden_vs_scan={rel_gap:.2e} min_near_boundary={'YES-SUSPECT' if (near_lo or near_hi) else 'no'}")
    return n_local_min, rel_gap, (near_lo or near_hi)

if __name__ == "__main__":
    mu = 398600.4418
    cases = [
        (6678.0, 6678.0*1.2, 179.0),   # near Hohmann
        (6678.0, 6678.0*1.5, 170.0),   # worked example case
        (6678.0, 6678.0*5.0, 150.0),   # high R, large deviation
        (6678.0, 6678.0*1.02, 179.0),  # near-degenerate R
        (6678.0, 6678.0*20.0, 150.0),  # extreme R from the 230-case grid
    ]
    print("Dense scan verification of global minimum (4000 samples/case):")
    all_ok = True
    for r1, r2, theta in cases:
        n_local, gap, suspect = scan_and_verify(r1, r2, theta, mu)
        if n_local != 1 or gap > 1e-6 or suspect:
            all_ok = False
    print("\nALL CASES CONFIRMED SINGLE GLOBAL MINIMUM, NOT NEAR BOUNDARY:", all_ok)
