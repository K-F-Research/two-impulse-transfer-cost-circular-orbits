"""
2i-OTCA-C v1.0 -- Independent numerical cross-check (Lambert p-iteration method)
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See the project LICENSE file / ATBD Annex A for the full license text.

Control B (figure): plots dV(p) across the admissible bracket for four
representative cases, marking the golden-section minimum and the bracket
boundaries, to show visually that each case has a single interior minimum
(not a boundary artifact). Saves the figure as PNG, SVG and PDF.
"""
import os
import mpmath as mp
mp.mp.dps = 30
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]
plt.rcParams["mathtext.fontset"] = "stix"
from lambert_independent_mp import transfer_cost_mp, _bracket_admissible_interval, min_cost_transfer_mp

mu = 398600.4418
cases = [
    (6678.0, 6678.0*1.2, 179.0, "R=1.2, theta=179 deg (near Hohmann)"),
    (6678.0, 6678.0*1.5, 170.0, "R=1.5, theta=170 deg (worked example, Sec. 9.1)"),
    (6678.0, 6678.0*5.0, 150.0, "R=5, theta=150 deg (large deviation)"),
    (6678.0, 6678.0*1.02, 179.0, "R=1.02, theta=179 deg (near-degenerate)"),
]

fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
for ax, (r1, r2, theta_deg, label) in zip(axes.flat, cases):
    theta = mp.radians(mp.mpf(theta_deg))
    ct, st = mp.cos(theta), mp.sin(theta)
    r1m, r2m = mp.mpf(r1), mp.mpf(r2)
    d = mp.sqrt(r1m**2 + r2m**2 - 2*r1m*r2m*ct)
    s = (r1m + r2m + d) / 2
    p_lo, p_hi = _bracket_admissible_interval(r1m, r2m, ct, st, s)
    eps = (p_hi - p_lo) * mp.mpf('1e-6')
    n = 600
    ps = [float(p_lo + eps + (p_hi - p_lo - 2*eps) * k / (n-1)) for k in range(n)]
    vals = [float(transfer_cost_mp(mp.mpf(p), r1m, r2m, theta, mp.mpf(mu))) for p in ps]
    p_opt, dv_opt, _ = min_cost_transfer_mp(r1, r2, theta_deg, mu)
    ax.plot(ps, vals, lw=1.3, color="#2b6cb0")
    ax.plot([p_opt], [dv_opt], 'o', color="#c0392b", zorder=5, ms=6,
            label=f"minimum found\np*={p_opt:.2f} km\ndV*={dv_opt:.6f} km/s")
    ax.axvline(p_lo, color="gray", ls="--", lw=0.8)
    ax.axvline(p_hi, color="gray", ls="--", lw=0.8)
    ax.set_title(label, fontsize=10)
    ax.set_xlabel("p (semi-latus rectum), km")
    ax.set_ylabel("dV(p), km/s")
    ax.legend(fontsize=7, loc="upper center")
    ax.grid(alpha=0.3)

fig.suptitle("Control B - Verification of the global minimum of dV(p) on the admissible interval\n"
              "(independent method, p-iteration / Shen & Tsiotras 1999)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])

out_dir = os.path.dirname(os.path.abspath(__file__))
for ext in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(out_dir, f"dv_curves_global_min_check.{ext}"),
                dpi=300 if ext == "png" else None)
print("saved: dv_curves_global_min_check.{png,svg,pdf}")
