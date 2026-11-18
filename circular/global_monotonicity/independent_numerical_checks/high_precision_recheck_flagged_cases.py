"""
Theorem 1 -- Global Monotonicity of the Optimized Transfer Cost
Arbitrary-Precision Re-Check of Flagged Boundary Cases

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

High-precision (mpmath, 50 digits) re-check of the apparent localization
anomalies found by random_stress_test.py, all in the extreme-R regime
combined with phi very close to 0 or pi -- the regime already known to
exceed double-precision resolution (see stress_test_high_precision.py).
"""

import mpmath as mp

mp.mp.dps = 50


def Sk_terms(p, phi, R):
    r1, r2 = mp.mpf(1), R
    vc1, vc2 = mp.mpf(1), 1 / mp.sqrt(R)
    dtheta = mp.pi - phi
    x = p / r1 - 1
    y = p / r2 - 1
    s = mp.sin(dtheta)
    bnum = x * mp.cos(dtheta) - y
    e2 = x ** 2 + (bnum / s) ** 2
    out = []
    for rk, vck in ((r1, vc1), (r2, vc2)):
        vthk = mp.sqrt(p) / rk
        vk2 = 2 / rk + (e2 - 1) / p
        out.append(vk2 + vck ** 2 - 2 * vck * vthk)
    return out


def feasible(p, phi, R):
    if p <= 0:
        return False
    S1, S2 = Sk_terms(p, phi, R)
    return S1 >= 0 and S2 >= 0


def F(p, phi, R):
    S1, S2 = Sk_terms(p, phi, R)
    return mp.sqrt(S1) + mp.sqrt(S2)


def dF_dp(p, phi, R, h=None):
    if h is None:
        h = p * mp.mpf(10) ** (-22)
    return (F(p + h, phi, R) - F(p - h, phi, R)) / (2 * h)


def find_pstar(phi, R, pH):
    mid0 = pH
    while not feasible(mid0, phi, R):
        mid0 = mid0 * mp.mpf("0.7")

    def fv(p):
        if not feasible(p, phi, R):
            return None
        return dF_dp(p, phi, R)

    hi, factor = mid0, mp.mpf("1.5")
    fhi = fv(hi)
    tries = 0
    while (fhi is None or fhi <= 0) and tries < 400:
        hi *= factor
        fhi = fv(hi)
        tries += 1
    lo = mid0
    flo = fv(lo)
    tries = 0
    while (flo is None or flo >= 0) and tries < 800:
        lo /= factor
        flo = fv(lo)
        tries += 1
    for _ in range(500):
        m = (lo + hi) / 2
        fm = fv(m)
        if fm is None:
            m = (m + mid0) / 2
            fm = fv(m)
            if fm is None:
                hi = m
                continue
        if fm < 0:
            lo = m
        else:
            hi = m
        if hi - lo < lo * mp.mpf(10) ** (-40):
            break
    return (lo + hi) / 2


def p1star(R, c):
    return R * (1 + c) / (1 + R * c)


def p2star(R, c):
    return R * (1 + c) / (R + c)


cases = [
    (60408.4, 179.935), (46063.5, 179.95), (21184.1, 0.054153),
    (94096.7, 179.81), (2650.04, 179.921), (31136.6, 179.947),
]

for R, deg in cases:
    R = mp.mpf(str(R)); deg = mp.mpf(str(deg))
    pH = 2 * R / (R + 1)
    phi = mp.radians(deg)
    c = mp.cos(phi)
    pstar = find_pstar(phi, R, pH)
    p2 = p2star(R, c)
    region = "I" if c > -1 / R else "II"
    if region == "I":
        p1 = p1star(R, c)
        loc_ok = (min(p1, p2) < pstar < max(p1, p2)) if p1 > 0 else (pstar > p2)
        extra = f"p1={float(p1):.6g} p2={float(p2):.6g}"
    else:
        loc_ok = pstar > p2
        extra = f"p2={float(p2):.6g}"
    print(f"R={float(R):10.3f} deg={float(deg):9.5f} region={region} pstar={float(pstar):.6e} {extra} loc_ok={loc_ok}")
