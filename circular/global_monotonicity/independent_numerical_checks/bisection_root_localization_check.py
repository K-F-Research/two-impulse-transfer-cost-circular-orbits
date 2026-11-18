"""
Theorem 1 -- Global Monotonicity of the Optimized Transfer Cost
Independent High-Precision Validation (bisection root-finding)

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Independent, higher-rigor numerical validation of the global monotonicity
theorem for V_R(phi): arbitrary precision (mpmath, 50 significant digits),
minimizer found via bisection on dF/dp=0 (robust version, avoiding the
complex-sqrt excursions that Newton's method can hit near degenerate
geometry as phi approaches pi), with an explicit check dF/dp(p*)~0 and
measurement of the signed distance of p* to p1*, p2*. Run over 9 values of
R in [1.2,50] and 19 angles from 1 deg to 179 deg (171 configurations).
"""

import mpmath as mp
mp.mp.dps = 50

def Sk_terms(p, phi, R):
    r1, r2 = mp.mpf(1), R
    vc1, vc2 = mp.mpf(1), 1/mp.sqrt(R)
    dtheta = mp.pi - phi
    x = p/r1 - 1
    y = p/r2 - 1
    s = mp.sin(dtheta)
    bnum = x*mp.cos(dtheta) - y
    e2 = x**2 + (bnum/s)**2
    out = []
    for rk, vck in ((r1, vc1), (r2, vc2)):
        vthk = mp.sqrt(p)/rk
        vk2 = 2/rk + (e2-1)/p
        Sk = vk2 + vck**2 - 2*vck*vthk
        out.append(Sk)
    return out, e2

def feasible(p, phi, R):
    (S1, S2), _ = Sk_terms(p, phi, R)
    return S1 >= 0 and S2 >= 0

def F(p, phi, R):
    (S1, S2), _ = Sk_terms(p, phi, R)
    return mp.sqrt(S1) + mp.sqrt(S2)

def dF_dp(p, phi, R, h=None):
    if h is None:
        h = p * mp.mpf(10)**(-25)
    return (F(p+h, phi, R) - F(p-h, phi, R)) / (2*h)

def dF_dphi(p, phi, R, h=None):
    if h is None:
        h = mp.mpf(10)**(-25)
    return (F(p, phi+h, R) - F(p, phi-h, R)) / (2*h)

def p1_star(R, c):
    return R*(1+c)/(1+R*c)

def p2_star(R, c):
    return R*(1+c)/(R+c)

def find_pstar_bisect(phi, R, pH):
    """Bisection on dF/dp's sign. dF/dp <0 below minimum, >0 above (F convex-ish
    in the feasible p-range around the minimum). Bracket search expands outward
    from pH while staying in the feasible region for both endpoints."""
    lo = pH * mp.mpf('0.05')
    hi = pH * mp.mpf('20')
    # ensure feasibility of the bracket endpoints; shrink hi if needed
    while not feasible(hi, phi, R) and hi > pH:
        hi = (hi + pH) / 2
    while not feasible(lo, phi, R) and lo < pH:
        lo = (lo + pH) / 2
    flo = dF_dp(lo, phi, R)
    fhi = dF_dp(hi, phi, R)
    if flo > 0 or fhi < 0:
        # fallback: scan for sign change
        n = 200
        pts = [lo + (hi-lo)*i/n for i in range(n+1)]
        vals = [dF_dp(pt, phi, R) if feasible(pt, phi, R) else None for pt in pts]
        for i in range(n):
            if vals[i] is not None and vals[i+1] is not None and vals[i] < 0 <= vals[i+1]:
                lo, hi = pts[i], pts[i+1]
                break
    for _ in range(200):
        mid = (lo+hi)/2
        if not feasible(mid, phi, R):
            # nudge toward pH
            mid = (mid + pH) / 2
        fm = dF_dp(mid, phi, R)
        if fm < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < pH * mp.mpf(10)**(-40):
            break
    return (lo+hi)/2

R_values = [mp.mpf(x) for x in ['1.2','1.5','2','3','5','8','11','20','50']]
deg_list = [1,5,10,20,30,45,60,75,89,91,105,120,135,150,160,170,175,178,179]

results = []
print(f"{'R':>6} {'deg':>5} {'p*':>14} {'dF/dp(p*)':>10} {'dF/dphi(p*)':>14} {'p1*':>13} {'p2*':>12} {'p*-p2*':>12} {'p1*-p*':>12} {'reg':>4} {'g sign':>7}")
for R in R_values:
    pH = 2*R/(R+1)
    for deg in deg_list:
        phi = mp.radians(deg)
        c = mp.cos(phi)
        pstar = find_pstar_bisect(phi, R, pH)
        dp = dF_dp(pstar, phi, R)
        dphi = dF_dphi(pstar, phi, R)
        p1 = p1_star(R, c)
        p2 = p2_star(R, c)
        region = "I" if c > -1/R else "II"
        x = pstar - 1
        y = pstar/R - 1
        g = (x*c + y)*(x + y*c)
        gs = "-" if g < 0 else ("+" if g > 0 else "0")
        ok_loc = None
        if region == "I":
            ok_loc = (p2 < pstar < p1)
        else:
            ok_loc = (pstar > p2)
        results.append((float(R), deg, region, gs, bool(dphi>0), ok_loc))
        print(f"{float(R):6.2f} {deg:5d} {float(pstar):14.8f} {float(dp):10.1e} {float(dphi):14.6e} {float(p1):13.6f} {float(p2):12.6f} {float(pstar-p2):12.6f} {float(p1-pstar):12.6f} {region:>4} {gs:>7}")

n_bad = sum(1 for r in results if r[3] != '-' or not r[4] or not r[5])
print(f"\nTotal configurations: {len(results)}; anomalies (g>=0 or dphi<=0 or wrong localization): {n_bad}")
