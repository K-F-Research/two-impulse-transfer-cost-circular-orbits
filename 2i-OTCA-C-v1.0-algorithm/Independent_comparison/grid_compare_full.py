"""
2i-OTCA-C v1.0 -- Independent numerical cross-check (Lambert p-iteration method)
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See the project LICENSE file / ATBD Annex A for the full license text.

Run B (230 cases, aligned with ATBD Annex B.3): compares the independent
p-iteration method against otca_c.api.compute_transfer on a 23 x 10 grid of
radius ratios R and angular deviations from Hohmann. See grid_compare.py
for the smaller, 36-case Run A.
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

R_list = [1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.6,1.75,1.9,2.0,2.25,2.5,2.75,3.0,3.5,4.0,4.5,5.0,7.0,10.0,20.0]
dev_list = [1,2,3,5,7,10,15,20,25,30]

rows = []
t0 = time.time()
n_fail = 0
for R in R_list:
    r2 = R * r1
    for dev in dev_list:
        theta = 180.0 - dev
        try:
            p_opt, dv_ind, bracket = min_cost_transfer_mp(r1, r2, theta, mu)
        except Exception as e:
            n_fail += 1
            print(f"FAIL R={R} dev={dev}: {e}")
            continue
        exact = api.compute_transfer(r1, r2, theta, mu_km3_s2=mu, validate=False)
        dv_core = exact.delta_v_total_km_s
        rel_err = abs(dv_ind - dv_core) / abs(dv_core)
        rows.append(dict(R=R, dev=dev, p_opt=p_opt, dv_independent=dv_ind,
                          dv_core=dv_core, rel_diff=rel_err))
print(f"elapsed: {time.time()-t0:.1f}s, cases={len(rows)}, failures={n_fail}")
out_dir = os.path.dirname(os.path.abspath(__file__))
json.dump(rows, open(os.path.join(out_dir, 'grid_compare_full_results.json'), 'w'), indent=2)
diffs = [r['rel_diff'] for r in rows]
print(f"max rel diff: {max(diffs):.3e}")
print(f"min rel diff: {min(diffs):.3e}")
print(f"mean rel diff: {sum(diffs)/len(diffs):.3e}")
worst = sorted(rows, key=lambda r: -r['rel_diff'])[:5]
print("worst 5 cases:")
for w in worst:
    print(w)
