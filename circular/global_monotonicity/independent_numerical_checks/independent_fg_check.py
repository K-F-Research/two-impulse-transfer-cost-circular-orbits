"""
Theorem 1 -- Global Monotonicity of the Optimized Transfer Cost
Independent Cross-Check of the Cost Function via Lagrange f,g Coefficients

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Re-derives the two-impulse Delta-v total via the classical Lagrange f,g
coefficients (vector Lambert formulation, Vallado/Curtis textbook route)
instead of the vis-viva + radial/tangential decomposition used throughout
the proof. Agreement between the two independently-coded formulas rules out
a shared bug risk in everything built on F(p,phi,R).
"""

import numpy as np

def F_original(p, phi, R):
    r1, r2 = 1.0, R
    vc1, vc2 = 1.0, 1/np.sqrt(R)
    dtheta = np.pi - phi
    x = p/r1 - 1
    y = p/r2 - 1
    s = np.sin(dtheta)
    bnum = x*np.cos(dtheta) - y
    e2 = x**2 + (bnum/s)**2
    total = 0.0
    for rk, vck in ((r1, vc1), (r2, vc2)):
        vthk = np.sqrt(p)/rk
        vk2 = 2/rk + (e2-1)/p
        Sk = vk2 + vck**2 - 2*vck*vthk
        total += np.sqrt(max(Sk, 0.0))
    return total

def F_lagrange_fg(p, phi, R):
    """Independent route: Lagrange f,g coefficients, vector velocities."""
    mu = 1.0
    r1, r2 = 1.0, R
    dtheta = np.pi - phi
    r1v = np.array([r1, 0.0])
    r2v = np.array([r2*np.cos(dtheta), r2*np.sin(dtheta)])
    h = np.sqrt(mu*p)
    f = 1 - (r2/p)*(1-np.cos(dtheta))
    g = r1*r2*np.sin(dtheta)/h
    gdot = 1 - (r1/p)*(1-np.cos(dtheta))
    v1v = (r2v - f*r1v)/g
    v2v = (gdot*r2v - r1v)/g
    # circular velocities, tangential (prograde, +90 deg from radius vector)
    def tangential(rv, vc):
        rhat = rv/np.linalg.norm(rv)
        that = np.array([-rhat[1], rhat[0]])
        return vc*that
    vc1v = tangential(r1v, 1.0)
    vc2v = tangential(r2v, 1/np.sqrt(R))
    dv1 = np.linalg.norm(v1v - vc1v)
    dv2 = np.linalg.norm(v2v - vc2v)
    return dv1+dv2

print(f"{'R':>6} {'deg':>6} {'p/pH':>6} {'F_original':>14} {'F_lagrange_fg':>14} {'|diff|':>12}")
maxdiff = 0.0
n=0
rng = np.random.default_rng(42)
for _ in range(500):
    R = rng.uniform(1.05, 50)
    deg = rng.uniform(1, 179)
    phi = np.radians(deg)
    pH = 2*R/(R+1)
    pfac = rng.uniform(0.3, 3.0)
    p = pfac*pH
    f1 = F_original(p, phi, R)
    f2 = F_lagrange_fg(p, phi, R)
    d = abs(f1-f2)
    maxdiff = max(maxdiff, d)
    n += 1
    if d > 1e-8:
        print(f"{R:6.2f} {deg:6.1f} {pfac:6.2f} {f1:14.10f} {f2:14.10f} {d:12.2e}")
print(f"\n{n} random draws tested. Maximum discrepancy |F_original - F_lagrange_fg| = {maxdiff:.3e}")
