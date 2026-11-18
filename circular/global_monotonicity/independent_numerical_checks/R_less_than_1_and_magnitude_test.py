"""
Theorem 1 -- Global Monotonicity of the Optimized Transfer Cost
Supplementary Numerical Checks: R<1 and Magnitude Consistency

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Two further checks, beyond the main localization grid:
(a) R<1 (r2<r1): the localization/monotonicity mechanism should carry over
    unchanged under simple relabeling of which body is "1" and which is "2"
    (the physics is symmetric in that choice). Tested here directly, without
    relabeling in the code, to catch any hidden asymmetry.
(b) Magnitude check: V_R(phi2)-V_R(phi1) should scale roughly linearly with
    (phi2-phi1) for small gaps (consistent with a locally bounded-below
    derivative), not just be positive in sign.
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


def VR(phi, R):
    pH = 2 * R / (R + 1)
    res = minimize_scalar(
        lambda p: F(p, phi, R), bounds=(1e-9 * pH, 200 * pH),
        method="bounded", options={"xatol": 1e-14},
    )
    return res.fun


if __name__ == "__main__":
    print("=== Test (a): R<1 (r2<r1), outside the range for which the "
          "theorem's code was written and proved (R>1) ===")
    ok = True
    for R in [0.9, 0.5, 0.2, 0.05, 0.01]:
        vals = [VR(np.radians(d), R) for d in [1, 30, 60, 90, 120, 150, 179]]
        monotone = all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
        print(f"  R={R:6.3f}  V(phi) increasing on the sample: {monotone}  "
              f"values={[round(v,5) for v in vals]}")
        if not monotone:
            ok = False
    print(f"  -> Monotonicity for R<1 (outside the proven domain, tested out of curiosity): {ok}")

    print("\n=== Test (b): magnitude of V_R(phi2)-V_R(phi1) for very small "
          "phi2-phi1 gaps ===")
    for R in [1.5, 10, 100]:
        for deg0 in [30, 90, 150]:
            phi0 = np.radians(deg0)
            for dh_deg in [1.0, 0.1, 0.01, 0.001]:
                dh = np.radians(dh_deg)
                dV = VR(phi0 + dh, R) - VR(phi0, R)
                rate = dV / dh
                print(f"  R={R:6.1f} phi0={deg0:4d}deg  dphi={dh_deg:7.4f}deg  "
                      f"dV={dV:.3e}  dV/dphi={rate:.6f}")
            print()
