"""
Theorem 1 -- Global Monotonicity of the Optimized Transfer Cost
Bracketing Confirmation: Sign of dF/dp at p1*, p2* (200 configurations)

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Confirms, at arbitrary precision (mpmath, 50 significant digits), that the
sign of dF/dp at the two closed-form points p1*, p2* alone suffices to
localize the true minimizer p* without ever solving for it explicitly: in
Region I, dF/dp(p2*) and dF/dp(p1*) have opposite signs (bracketing p* via
the intermediate value theorem, since F is convex in p near the minimum);
in Region II, dF/dp(p2*)<0 places the minimum to the right of p2*. Run over
200 configurations (10 values of R in [1.2,100] x 20 angles from 1 deg to
179.9 deg), consistent with Lemma 7.3 of the proof document.
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

def F(p, phi, R):
    S1, S2 = Sk_terms(p, phi, R)
    return mp.sqrt(S1) + mp.sqrt(S2)

def dF_dp(p, phi, R, h=None):
    if h is None:
        h = p * mp.mpf(10)**(-25)
    return (F(p+h, phi, R) - F(p-h, phi, R)) / (2*h)

def p1_star(R, c):
    return R*(1+c)/(1+R*c)

def p2_star(R, c):
    return R*(1+c)/(R+c)

R_values = [mp.mpf(x) for x in ['1.2','1.5','2','3','5','8','11','20','50','100']]
deg_list = [1,5,10,20,30,45,60,75,89,91,105,120,135,150,160,170,175,178,179,179.9]

print(f"{'R':>6} {'deg':>7} {'reg':>4} {'dF/dp(p2*)':>14} {'dF/dp(p1*)':>14}  bracket_ok")
n_bad = 0
n_tot = 0
for R in R_values:
    for deg in deg_list:
        phi = mp.radians(mp.mpf(str(deg)))
        c = mp.cos(phi)
        p2 = p2_star(R, c)
        region = "I" if c > -1/R else "II"
        dFp2 = dF_dp(p2, phi, R)
        if region == "I":
            p1 = p1_star(R, c)
            dFp1 = dF_dp(p1, phi, R)
            ok = (dFp2 < 0) and (dFp1 > 0)
            print(f"{float(R):6.1f} {float(deg):7.2f} {region:>4} {float(dFp2):14.6e} {float(dFp1):14.6e}  {ok}")
        else:
            dFp1 = None
            ok = (dFp2 < 0)
            print(f"{float(R):6.1f} {float(deg):7.2f} {region:>4} {float(dFp2):14.6e} {'  (p1<0, n/a)':>14}  {ok}")
        n_tot += 1
        if not ok:
            n_bad += 1
print(f"\ntotal={n_tot}  failures={n_bad}")
