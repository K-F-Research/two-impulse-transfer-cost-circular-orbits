"""
Theorem 1 -- Global Monotonicity of the Optimized Transfer Cost
Independent High-Precision Validation (near phi -> pi refinement)

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Refinement of bisection_root_localization_check.py's bisection-based
validation, extending the lower search bracket so it can shrink well below
0.05*pH: this is needed because the true minimizer p*(phi) shrinks quickly
as phi -> pi (near-degenerate transfer geometry), which otherwise corrupts
results in that regime. Run over 9 values of R in [1.2,50] and 8 angles
from 150 deg to 179.9 deg (72 configurations). Combined with
bisection_root_localization_check.py, this gives the 243 configurations
(angles from 1 deg to 179.9 deg) cited as the arbitrary-precision
localization check in the proof document.
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
    return out

def feasible(p, phi, R):
    if p <= 0:
        return False
    (S1, S2) = Sk_terms(p, phi, R)
    return S1 >= 0 and S2 >= 0

def F(p, phi, R):
    (S1, S2) = Sk_terms(p, phi, R)
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
    # locate a feasible starting point near pH, then find a bracket [lo,hi]
    # with dF_dp(lo)<0<dF_dp(hi) by geometric expansion in both directions,
    # shrinking toward feasibility as needed.
    mid0 = pH
    while not feasible(mid0, phi, R):
        mid0 = mid0 * mp.mpf('0.7')
        if mid0 < pH*mp.mpf(10)**(-20):
            raise RuntimeError("no feasible start found")

    def feas_val(p):
        if not feasible(p, phi, R):
            return None
        return dF_dp(p, phi, R)

    # expand upward for hi with dF_dp>0
    hi = mid0
    fhi = feas_val(hi)
    factor = mp.mpf('1.5')
    tries = 0
    while (fhi is None or fhi <= 0) and tries < 200:
        hi = hi * factor
        fhi = feas_val(hi)
        tries += 1
    # expand downward for lo with dF_dp<0
    lo = mid0
    flo = feas_val(lo)
    tries = 0
    while (flo is None or flo >= 0) and tries < 400:
        lo = lo / factor
        flo = feas_val(lo)
        tries += 1
    if flo is None or fhi is None or flo >= 0 or fhi <= 0:
        raise RuntimeError(f"bracket failed phi={phi} R={R} lo={lo} hi={hi} flo={flo} fhi={fhi}")
    for _ in range(300):
        m = (lo+hi)/2
        fm = feas_val(m)
        if fm is None:
            # infeasible midpoint (can happen if feasible region is disconnected
            # in a weird way) -- nudge toward mid0
            m = (m+mid0)/2
            fm = feas_val(m)
            if fm is None:
                hi = m  # shrink conservatively
                continue
        if fm < 0:
            lo = m
        else:
            hi = m
        if hi-lo < (lo)*mp.mpf(10)**(-40):
            break
    return (lo+hi)/2

R_values = [mp.mpf(x) for x in ['1.2','1.5','2','3','5','8','11','20','50']]
deg_list = [150,160,170,175,178,179,179.5,179.9]

print(f"{'R':>6} {'deg':>7} {'p*':>16} {'dF/dp(p*)':>10} {'dF/dphi(p*)':>14} {'p1*':>13} {'p2*':>12} {'p*-p2*':>14} {'reg':>4} {'g sign':>7}")
anomalies = 0
total = 0
for R in R_values:
    pH = 2*R/(R+1)
    for deg in deg_list:
        phi = mp.radians(mp.mpf(str(deg)))
        c = mp.cos(phi)
        try:
            pstar = find_pstar_bisect(phi, R, pH)
        except RuntimeError as e:
            print("FAIL", R, deg, e)
            continue
        dp = dF_dp(pstar, phi, R)
        dphi = dF_dphi(pstar, phi, R)
        p1 = p1_star(R, c)
        p2 = p2_star(R, c)
        region = "I" if c > -1/R else "II"
        x = pstar - 1
        y = pstar/R - 1
        g = (x*c + y)*(x + y*c)
        gs = "-" if g < 0 else ("+" if g > 0 else "0")
        total += 1
        if gs != "-" or not (dphi > 0):
            anomalies += 1
        print(f"{float(R):6.2f} {float(deg):7.2f} {float(pstar):16.10f} {float(dp):10.1e} {float(dphi):14.6e} {float(p1):13.6f} {float(p2):12.6f} {float(pstar-p2):14.8f} {region:>4} {gs:>7}")
print(f"\ntotal={total} anomalies={anomalies}")
