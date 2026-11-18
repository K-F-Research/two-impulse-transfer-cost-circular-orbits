"""
Theorem 1 -- Global Monotonicity of the Optimized Transfer Cost
Independent High-Precision Validation (Newton root-finding)

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Independent, higher-rigor numerical validation of the global monotonicity
theorem for V_R(phi): arbitrary precision (mpmath, 50 significant digits),
minimizer found via root-finding on dF/dp=0 (Newton's method) rather than
scipy.optimize.minimize_scalar, with an explicit check dF/dp(p*)~0 and
measurement of the signed distance of p* to p1*, p2*.
"""

import mpmath as mp
mp.mp.dps = 50  # 50 decimal digits

def F(p, phi, R):
    """Total Δv(p, phi; R) with r1=1, r2=R, vc1=1, vc2=1/sqrt(R)."""
    r1, r2 = mp.mpf(1), R
    vc1, vc2 = mp.mpf(1), 1/mp.sqrt(R)
    dtheta = mp.pi - phi
    x = p/r1 - 1
    y = p/r2 - 1
    s = mp.sin(dtheta)
    bnum = x*mp.cos(dtheta) - y
    e2 = x**2 + (bnum/s)**2
    total = mp.mpf(0)
    for rk, vck in ((r1, vc1), (r2, vc2)):
        vthk = mp.sqrt(p)/rk
        vk2 = 2/rk + (e2-1)/p
        Sk = vk2 + vck**2 - 2*vck*vthk
        total += mp.sqrt(Sk)
    return total

def dF_dp(p, phi, R, h=None):
    if h is None:
        h = mp.mpf(10)**(-30)
    return (F(p+h, phi, R) - F(p-h, phi, R)) / (2*h)

def dF_dphi(p, phi, R, h=None):
    if h is None:
        h = mp.mpf(10)**(-30)
    return (F(p, phi+h, R) - F(p, phi-h, R)) / (2*h)

def p1_star(R, c):
    return R*(1+c)/(1+R*c)

def p2_star(R, c):
    return R*(1+c)/(R+c)

def find_pstar(phi, R, pH):
    """Root-find dF/dp = 0 via mpmath.findroot (secant/Newton via mp.findroot),
    seeded near Hohmann value pH -- independent of scipy.optimize.minimize_scalar."""
    f = lambda p: dF_dp(p, phi, R)
    p0 = pH
    try:
        root = mp.findroot(f, p0, solver='newton', tol=mp.mpf(10)**(-35))
    except Exception:
        root = mp.findroot(f, (p0*mp.mpf('0.5'), p0*mp.mpf('2')), solver='anderson')
    return root

R_values = [mp.mpf(x) for x in ['1.2','1.5','2','3','5','8','11','20','50']]
deg_list = [1,5,10,20,30,45,60,75,89,91,105,120,135,150,160,170,175,179]

print(f"{'R':>6} {'deg':>5} {'p*':>14} {'dF/dp(p*)':>14} {'dF/dphi(p*)':>14} {'p1*':>12} {'p2*':>12} {'p*-p2*':>12} {'p1*-p*':>12} {'region':>8} {'g(p*) sign':>10}")
for R in R_values:
    pH = 2*R/(R+1)
    for deg in deg_list:
        phi = mp.radians(deg)
        c = mp.cos(phi)
        pstar = find_pstar(phi, R, pH)
        dp = dF_dp(pstar, phi, R)
        dphi = dF_dphi(pstar, phi, R)
        p1 = p1_star(R, c)
        p2 = p2_star(R, c)
        region = "I" if c > -1/R else "II"
        # g(p) = (x c + y)(x + y c), x = p-1, y = p/R -1
        x = pstar - 1
        y = pstar/R - 1
        g = (x*c + y)*(x + y*c)
        gs = "-" if g < 0 else ("+" if g > 0 else "0")
        print(f"{float(R):6.1f} {deg:5d} {float(pstar):14.8f} {float(dp):14.2e} {float(dphi):14.6e} {float(p1):12.6f} {float(p2):12.6f} {float(pstar-p2):12.6f} {float(p1-pstar):12.6f} {region:>8} {gs:>10}")
