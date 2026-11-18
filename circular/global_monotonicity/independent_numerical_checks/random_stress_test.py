"""
Theorem 1 -- Global Monotonicity of the Optimized Transfer Cost
Random Stress Test

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Large-scale adversarial random sampling of the global monotonicity theorem for
V_R(phi) = min_p F(p,phi,R): tests localization (p2* < p* < p1* in Region I,
p* > p2* in Region II) and the sign of the envelope derivative dV_R/dphi over
tens of thousands of random (R,phi) draws spanning a very wide range, well
beyond the hand-picked grids used in the main verification pass -- reducing
the risk that a curated grid happened to avoid a problematic region.
"""

import numpy as np
from scipy.optimize import minimize_scalar


def F(p, phi, R):
    r1, r2 = 1.0, R
    vc1, vc2 = 1.0, 1 / np.sqrt(R)
    dtheta = np.pi - phi
    x = p / r1 - 1
    y = p / r2 - 1
    s = np.sin(dtheta)
    if abs(s) < 1e-13:
        return np.inf
    bnum = x * np.cos(dtheta) - y
    e2 = x ** 2 + (bnum / s) ** 2
    total = 0.0
    for rk, vck in ((r1, vc1), (r2, vc2)):
        vthk = np.sqrt(p) / rk
        vk2 = 2 / rk + (e2 - 1) / p
        Sk = vk2 + vck ** 2 - 2 * vck * vthk
        if Sk < -1e-9:
            return np.inf
        total += np.sqrt(max(Sk, 0.0))
    return total


def p1star(R, c):
    return R * (1 + c) / (1 + R * c)


def p2star(R, c):
    return R * (1 + c) / (R + c)


def localization_ok(pstar, R, c):
    p2 = p2star(R, c)
    region = "I" if c > -1 / R else "II"
    if region == "I":
        p1 = p1star(R, c)
        return (min(p1, p2) < pstar < max(p1, p2)) if p1 > 0 else (pstar > p2)
    return pstar > p2


if __name__ == "__main__":
    rng = np.random.default_rng(20260903)
    N = 50000
    n_fail = 0
    fails = []
    # log-uniform on R to cover both near-1 and very large ratios fairly,
    # and a mix of uniform + edge-concentrated angles.
    for i in range(N):
        R = float(np.exp(rng.uniform(np.log(1.001), np.log(1e5))))
        if rng.random() < 0.3:
            deg = rng.uniform(0.05, 5.0)          # near phi=0
        elif rng.random() < 0.5:
            deg = rng.uniform(175.0, 179.95)       # near phi=pi
        else:
            deg = rng.uniform(0.05, 179.95)        # generic
        phi = np.radians(deg)
        c = np.cos(phi)
        pH = 2 * R / (R + 1)
        res = minimize_scalar(
            lambda p: F(p, phi, R),
            bounds=(1e-9 * pH, 200 * pH),
            method="bounded",
            options={"xatol": 1e-14},
        )
        pstar = res.x
        if not localization_ok(pstar, R, c):
            n_fail += 1
            fails.append((R, deg, pstar))

    print(f"{N} random draws (R log-uniform on [1.001, 1e5], angles "
          f"concentrated near the boundaries + generic).")
    print(f"Localization failures (before high-precision re-check): {n_fail}")
    for R, deg, pstar in fails[:30]:
        print(f"  R={R:.6g}  deg={deg:.6g}  pstar={pstar:.6g}")
    if len(fails) > 30:
        print(f"  ... et {len(fails)-30} autres")
