"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Example: visualize a two-impulse circular-to-circular transfer solved by
otca_c, and the cost-vs-angle curve this repository's Theorem 1 proves is
strictly monotonic on either side of the Hohmann angle.

Run with: python3 visualize_transfer.py
Requires: matplotlib (in addition to otca_c's own numpy/scipy/sympy).
Writes: transfer_visualization.png in the current directory.
"""
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from otca_c import Bodies, compute_transfer
from otca_c.core import find_p_star, cost as normalized_cost


def orbit_xy(p, e, omega, nu_array):
    """Classical conic r(nu) = p/(1+e cos(nu)), rotated by the orbit's own
    argument of periapsis omega, for plotting in a fixed inertial frame."""
    r = p / (1 + e * np.cos(nu_array))
    theta = nu_array + omega
    return r * np.cos(theta), r * np.sin(theta)


def main():
    r1_km = 6778.0  # ~400 km altitude LEO (Earth radius 6378 km)
    r2_km = 13556.0  # R = 2
    delta_theta_deg = 170.0

    result = compute_transfer(r1_km, r2_km, delta_theta_deg, mu_km3_s2=Bodies.EARTH)
    print(result.summary())

    R = result.R
    phi = result.phi_rad
    p_norm = result.p_star_km / r1_km
    e = result.eccentricity

    # Departure true anomaly nu1 on the transfer orbit (same derivation as
    # otca_c.geometry.independent_check), used only to orient the plot.
    c, s = math.cos(phi), math.sin(phi)
    x = p_norm - 1.0
    y = p_norm / R - 1.0
    e_sin_nu1 = -(x * c + y) / s
    nu1 = math.atan2(e_sin_nu1, x)
    delta_theta_rad = math.pi - abs(phi)
    nu2 = nu1 + delta_theta_rad
    omega = -nu1  # so the departure point sits at inertial angle 0

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    # --- Left panel: transfer geometry -------------------------------
    ax = axes[0]
    theta_full = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(theta_full), np.sin(theta_full), "--", color="0.6", label="departure orbit (r1)")
    ax.plot(R * np.cos(theta_full), R * np.sin(theta_full), "--", color="0.4", label="arrival orbit (r2)")

    nu_arc = np.linspace(nu1, nu2, 200)
    tx, ty = orbit_xy(p_norm, e, omega, nu_arc)
    ax.plot(tx, ty, color="tab:orange", linewidth=2, label="transfer orbit arc")

    dep_x, dep_y = orbit_xy(p_norm, e, omega, np.array([nu1]))
    arr_x, arr_y = orbit_xy(p_norm, e, omega, np.array([nu2]))
    ax.scatter(dep_x, dep_y, color="tab:green", zorder=5, label="departure burn")
    ax.scatter(arr_x, arr_y, color="tab:red", zorder=5, label="arrival burn")
    ax.scatter([0], [0], color="black", marker="+", s=100)

    ax.set_aspect("equal")
    ax.set_title(
        f"Transfer geometry (R={R:.3f}, "
        f"delta_theta={delta_theta_deg:.0f} deg)\n"
        f"p*={result.p_star_km:.0f} km, e={e:.4f}"
    )
    ax.legend(loc="upper left", fontsize=8)

    # --- Right panel: cost vs transfer angle --------------------------
    ax = axes[1]
    angles_deg = np.linspace(2, 178, 180)
    costs = []
    for a_deg in angles_deg:
        p_star = find_p_star(math.radians(a_deg), R)
        costs.append(normalized_cost(p_star, math.radians(a_deg), R))
    costs = np.array(costs)
    hohmann_cost = normalized_cost(find_p_star(math.radians(0.5), R), math.radians(0.5), R)

    ax.plot(180 - angles_deg, costs, color="tab:blue")
    ax.axvline(180.0, color="0.7", linestyle=":", label="Hohmann (delta_theta=180 deg)")
    ax.axvline(delta_theta_deg, color="tab:orange", linestyle="--", label="this transfer")
    ax.set_xlabel("delta_theta (deg)")
    ax.set_ylabel("optimized cost V_R(phi)  [units of vc1]")
    ax.set_title(
        "Cost vs. transfer angle (Theorem 1: strictly\n"
        "monotonic on each side of the Hohmann angle)"
    )
    ax.legend(loc="upper center", fontsize=8)

    fig.tight_layout()
    fig.savefig("transfer_visualization.png", dpi=150)
    print("\nSaved transfer_visualization.png")


if __name__ == "__main__":
    main()
