"""
Theorem 1 -- Global Monotonicity of the Optimized Transfer Cost
Stress Test (Double Precision)

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Searches for counterexamples or pathological regimes in four families:
R -> 1+, R -> infinity, multiple local minima in p at fixed (phi,R), and
phi -> 0+ / phi -> pi- (outside the theorem's stated domain, tested for
caution). Any apparent failure found here must be re-checked at arbitrary
precision (see stress_test_high_precision.py) before being treated as a
real anomaly: double precision is known to lose resolution when p1*(phi,R)
and p2*(phi,R) become close (near R=1 or near the angular endpoints).
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
    if abs(s) < 1e-14:
        return np.inf
    bnum = x * np.cos(dtheta) - y
    e2 = x ** 2 + (bnum / s) ** 2
    total = 0.0
    for rk, vck in ((r1, vc1), (r2, vc2)):
        vthk = np.sqrt(p) / rk
        vk2 = 2 / rk + (e2 - 1) / p
        Sk = vk2 + vck ** 2 - 2 * vck * vthk
        if Sk < -1e-10:
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


def find_pstar(phi, R, pH):
    res = minimize_scalar(
        lambda p: F(p, phi, R),
        bounds=(1e-8 * pH, 50 * pH),
        method="bounded",
        options={"xatol": 1e-14},
    )
    return res.x


if __name__ == "__main__":
    print("=== Test 1: R -> 1+ ===")
    for R in [1.5, 1.1, 1.01, 1.001, 1.0001]:
        pH = 2 * R / (R + 1)
        ok = True
        for deg in [1, 45, 90, 135, 179]:
            phi = np.radians(deg)
            c = np.cos(phi)
            pstar = find_pstar(phi, R, pH)
            if not localization_ok(pstar, R, c):
                ok = False
        print(f"R={R:9.5f}  correct localization: {ok}")

    print("\n=== Test 2: R -> infinity ===")
    for R in [10, 100, 1e3, 1e4, 1e5, 1e6]:
        pH = 2 * R / (R + 1)
        ok = True
        for deg in [1, 45, 90, 135, 179]:
            phi = np.radians(deg)
            c = np.cos(phi)
            pstar = find_pstar(phi, R, pH)
            if not localization_ok(pstar, R, c):
                ok = False
        print(f"R={R:12.1f}  correct localization: {ok}")

    print("\n=== Test 3: multiple local minima (fine scan of F(p)) ===")
    tested, bad_multi = 0, 0
    for R in [1.2, 2, 5, 20, 100]:
        pH = 2 * R / (R + 1)
        for deg in [5, 30, 60, 90, 120, 150, 175]:
            phi = np.radians(deg)
            pgrid = np.geomspace(1e-6 * pH, 100 * pH, 4000)
            vals = np.array([F(p, phi, R) for p in pgrid])
            finite = np.isfinite(vals)
            v = vals[finite]
            dv = np.diff(v)
            n_local_min = len(np.where(np.diff(np.sign(dv)) > 0)[0])
            tested += 1
            if n_local_min > 1:
                bad_multi += 1
                print(f"  R={R}, deg={deg}: {n_local_min} local minima detected")
    print(f"Configurations tested: {tested}, with >1 local minimum: {bad_multi}")

    print("\n=== Test 4: phi -> 0+ and phi -> pi- (outside theorem scope) ===")
    for R in [1.5, 5, 50]:
        pH = 2 * R / (R + 1)
        for deg in [0.01, 0.001, 179.99, 179.999]:
            phi = np.radians(deg)
            c = np.cos(phi)
            pstar = find_pstar(phi, R, pH)
            print(f"  R={R}, deg={deg}: loc_ok={localization_ok(pstar, R, c)}, pstar={pstar:.6g}")
