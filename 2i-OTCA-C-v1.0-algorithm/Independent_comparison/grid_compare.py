"""
2i-OTCA-C v1.0 -- Independent numerical cross-check (Lambert p-iteration method)
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See the project LICENSE file / ATBD Annex A for the full license text.

Run A (36 cases): compares the independent p-iteration method against
otca_c.api.compute_transfer on a 6 x 6 grid of radius ratios R and angular
deviations from Hohmann. See grid_compare_full.py for the larger, 230-case
Run B that matches ATBD Annex B.3 exactly.
"""
import os
import sys
import json
import time
# Adjust this path to the otca_c package's src/ directory if run standalone.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from otca_c import api
from lambert_independent_mp import min_cost_transfer_mp

mu = 398600.4418
r1 = 6678.0

R_list = [1.2, 1.5, 2.0, 3.0, 5.0, 10.0]
dev_list = [5, 10, 15, 20, 25, 30]

rows = []
t0 = time.time()
for R in R_list:
    r2 = R * r1
    for dev in dev_list:
        theta = 180.0 - dev
        p_opt, dv_ind, bracket = min_cost_transfer_mp(r1, r2, theta, mu)
        exact = api.compute_transfer(r1, r2, theta, mu_km3_s2=mu, validate=False)
        dv_core = exact.delta_v_total_km_s
        rel_err = abs(dv_ind - dv_core) / abs(dv_core)
        rows.append(dict(R=R, dev=dev, p_opt=p_opt, dv_independent=dv_ind,
                          dv_core=dv_core, rel_diff=rel_err))
        print(f"R={R:<6} dev={dev:<4} dV_ind={dv_ind:.10f}  dV_core={dv_core:.10f}  rel_diff={rel_err:.3e}")
print(f"elapsed: {time.time()-t0:.1f}s")
out_dir = os.path.dirname(os.path.abspath(__file__))
json.dump(rows, open(os.path.join(out_dir, 'grid_compare_results.json'), 'w'), indent=2)
max_err = max(r['rel_diff'] for r in rows)
print(f"\nMax relative difference across {len(rows)} cases: {max_err:.3e}")
