# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
"""
certify_Ppp_v2.py -- corrected version, following an audit finding.

Correction relative to certify_Ppp.py: the earlier version bounded G_Y,
Y', Y'', P'' on a box [Ylo,Yhi] chosen heuristically (from the known
branch slope), WITHOUT ever verifying that the true branch Y(w) stays
inside that box for every w in the tested W-box.

This version adds the missing check: before accepting the P'' certificate,
it verifies, by sign-bracketing of G itself (not only of its derivatives),
that:
    (i)  G(w, Y_lo) > 0   for every w in W-box
    (ii) G(w, Y_hi) < 0   for every w in W-box
    (iii) G_Y(w,Y) < 0 strictly for every (w,Y) in W-box x [Ylo,Yhi]
(i)+(iii) => existence of a unique zero of G(w,.) in (Ylo,Yhi) for each
fixed w (G continuous, strictly decreasing in Y, changes sign); (ii)
confirms the sign change on the other side. Combined, (i)-(iii) certify
that the true branch Y(w) (the root of G(w,.)=0) is indeed contained in
[Ylo,Yhi] for every w in the W-box -- which was the unverified hypothesis
of the previous version.

Architecture kept identical (same construction of G, Q_ww, Q_wY, Q_YY,
Q_Y via the chain rule, same certify() function); only the certify()
function is extended with the missing check, BEFORE computing/accepting
Ppp_box.
"""
import sympy as sp
import mpmath as mp
import pickle, time

w, Y, u1, u2 = sp.symbols('w Y u1 u2', positive=True)
with open('/tmp/Q_second_partials.pkl','rb') as f:
    d = pickle.load(f)
Q_ww, Q_wY, Q_YY, Q_Y = d['Q_ww'], d['Q_wY'], d['Q_YY'], d['Q_Y']

S1 = (8*Y**2*w**3 - 12*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 16*Y*w + 12*Y + 1) / (4*Y*(2 - Y*w**2))
S2 = (2*sp.sqrt(2)*Y**2*w**3 - 6*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 4*sp.sqrt(2)*Y*w + 1) / (4*Y*(2 - Y*w**2))

p, x = sp.symbols('p x')
x1 = p - 1; x2 = p/2 - 1; s2 = 1 - x**2
e2 = x1**2 + (x1*x + x2)**2 / s2
v1sq = 2 + (e2 - 1)/p; v2sq = 1 + (e2 - 1)/p
S1_orig = v1sq + 1 - 2*sp.sqrt(p)
S2_orig = v2sq + sp.Rational(1,2) - sp.sqrt(p/2)
S1_p_orig = sp.diff(S1_orig, p); S2_p_orig = sp.diff(S2_orig, p)
subs = {p: w**2, x: -1 + w**2*Y}
S1_p_reg = sp.simplify(S1_p_orig.subs(subs, simultaneous=True))
S2_p_reg = sp.simplify(S2_p_orig.subs(subs, simultaneous=True))
Cp_reg = S1_p_reg/(2*sp.sqrt(S1)) + S2_p_reg/(2*sp.sqrt(S2))
G = sp.together(w**2 * Cp_reg)

# G itself needs to be evaluable in interval mode with sqrt(S1), sqrt(S2)
# explicitly positive -> re-expressed via u1=sqrt(S1), u2=sqrt(S2) to reuse
# exactly the same radicals as Q_ww etc. (strict consistency).
G_u = sp.together(w**2 * (S1_p_reg/(2*u1) + S2_p_reg/(2*u2)))
# check: G_u with u1->sqrt(S1), u2->sqrt(S2) must reproduce G
check = sp.simplify(G_u.subs({u1: sp.sqrt(S1), u2: sp.sqrt(S2)}) - G)
assert check == 0, f"G_u substitution mismatch: {check}"

G_w = sp.diff(G,w); G_Y = sp.diff(G,Y)
G_ww = sp.diff(G,w,w); G_wY = sp.diff(G,w,Y); G_YY = sp.diff(G,Y,Y)

mp.mp.dps = 50
mp.iv.dps = 50
modules_iv = [{'sqrt': mp.iv.sqrt, 'mpf': mp.iv.mpf}, 'mpmath']

Q_ww_f = sp.lambdify((w,Y,u1,u2), Q_ww, modules=modules_iv)
Q_wY_f = sp.lambdify((w,Y,u1,u2), Q_wY, modules=modules_iv)
Q_YY_f = sp.lambdify((w,Y,u1,u2), Q_YY, modules=modules_iv)
Q_Y_f  = sp.lambdify((w,Y,u1,u2), Q_Y, modules=modules_iv)
Gw_f  = sp.lambdify((w,Y), G_w, modules=modules_iv)
GY_f  = sp.lambdify((w,Y), G_Y, modules=modules_iv)
Gww_f = sp.lambdify((w,Y), G_ww, modules=modules_iv)
GwY_f = sp.lambdify((w,Y), G_wY, modules=modules_iv)
GYY_f = sp.lambdify((w,Y), G_YY, modules=modules_iv)
S1_f  = sp.lambdify((w,Y), S1, modules=modules_iv)
S2_f  = sp.lambdify((w,Y), S2, modules=modules_iv)
# G_u evaluated directly with explicit u1,u2 (same sqrt-branch convention
# as everywhere else) -- used ONLY for the sign-bracketing check at
# Y=Ylo/Yhi.
Gu_f  = sp.lambdify((w,Y,u1,u2), G_u, modules=modules_iv)


def iv_pos_sqrt(box):
    """sqrt of an interval box known/required to be strictly positive."""
    if box.a <= 0:
        raise ValueError(f"radicand not certified positive: {box}")
    return mp.iv.sqrt(box)


def certify(w_lo, w_hi, Y_lo, Y_hi, label, verbose=True):
    Wbox = mp.iv.mpf([w_lo, w_hi])
    Ybox = mp.iv.mpf([Y_lo, Y_hi])

    # ---- (0) sanity: S1,S2 > 0 on the full box (needed for all sqrt calls) ----
    S1_box = S1_f(Wbox, Ybox); S2_box = S2_f(Wbox, Ybox)
    if S1_box.a <= 0 or S2_box.a <= 0:
        if verbose: print(f"[{label}] ABORT: S1 or S2 not certified positive on full box")
        return None
    u1_box = mp.iv.sqrt(S1_box); u2_box = mp.iv.sqrt(S2_box)

    # ---- (1) NEW: branch-inclusion check -----------------------------------
    # G_Y strictly negative over the WHOLE box (w,Y) in Wbox x Ybox
    GY_box_full = GY_f(Wbox, Ybox)
    if not (GY_box_full.b < 0):
        if verbose: print(f"[{label}] ABORT: G_Y not certified strictly negative on full box: {GY_box_full}")
        return None

    # G(w, Y_lo) > 0 for all w in Wbox  (Y_lo is the point Y=Y_lo, zero-width interval)
    Ylo_pt = mp.iv.mpf(Y_lo)
    S1_lo = S1_f(Wbox, Ylo_pt); S2_lo = S2_f(Wbox, Ylo_pt)
    if S1_lo.a <= 0 or S2_lo.a <= 0:
        if verbose: print(f"[{label}] ABORT: S1/S2 not positive at Y=Ylo edge")
        return None
    u1_lo = mp.iv.sqrt(S1_lo); u2_lo = mp.iv.sqrt(S2_lo)
    G_at_Ylo = Gu_f(Wbox, Ylo_pt, u1_lo, u2_lo)
    if not (G_at_Ylo.a > 0):
        if verbose: print(f"[{label}] ABORT: G(w,Ylo) not certified > 0 for all w in box: {G_at_Ylo}")
        return None

    # G(w, Y_hi) < 0 for all w in Wbox
    Yhi_pt = mp.iv.mpf(Y_hi)
    S1_hi = S1_f(Wbox, Yhi_pt); S2_hi = S2_f(Wbox, Yhi_pt)
    if S1_hi.a <= 0 or S2_hi.a <= 0:
        if verbose: print(f"[{label}] ABORT: S1/S2 not positive at Y=Yhi edge")
        return None
    u1_hi = mp.iv.sqrt(S1_hi); u2_hi = mp.iv.sqrt(S2_hi)
    G_at_Yhi = Gu_f(Wbox, Yhi_pt, u1_hi, u2_hi)
    if not (G_at_Yhi.b < 0):
        if verbose: print(f"[{label}] ABORT: G(w,Yhi) not certified < 0 for all w in box: {G_at_Yhi}")
        return None

    if verbose:
        print(f"[{label}] branch-inclusion OK: G(w,Ylo)={G_at_Ylo} > 0, "
              f"G(w,Yhi)={G_at_Yhi} < 0, G_Y<=0 on box: {GY_box_full}")

    # ---- (2) now (and only now) proceed to the P'' bound, as before --------
    GY_box = GY_box_full
    Gw_box = Gw_f(Wbox, Ybox)
    Gww_box = Gww_f(Wbox, Ybox)
    GwY_box = GwY_f(Wbox, Ybox)
    GYY_box = GYY_f(Wbox, Ybox)

    Yp_box = -Gw_box/GY_box
    Ypp_box = -(Gww_box + 2*GwY_box*Yp_box + GYY_box*Yp_box**2)/GY_box

    Qww_box = Q_ww_f(Wbox, Ybox, u1_box, u2_box)
    QwY_box = Q_wY_f(Wbox, Ybox, u1_box, u2_box)
    QYY_box = Q_YY_f(Wbox, Ybox, u1_box, u2_box)
    QY_box  = Q_Y_f(Wbox, Ybox, u1_box, u2_box)
    Ppp_box = Qww_box + 2*QwY_box*Yp_box + QYY_box*Yp_box**2 + QY_box*Ypp_box

    if verbose:
        print(f"[{label}] CERTIFIED (branch-inclusion proved): P'' in {Ppp_box}  "
              f"(sup P'' = {Ppp_box.b})")
    return Ppp_box


if __name__ == "__main__":
    t0 = time.time()
    # Same box as the previous version, now with the branch-inclusion
    # check added before accepting the result.
    b_coef = mp.mpf('0.57')
    wstar = '0.00015'
    wv = mp.mpf(wstar)
    margin = b_coef*wv*mp.mpf('1.5')
    Ylo = mp.mpf('0.5') - margin
    Yhi = mp.mpf('0.5') + margin*mp.mpf('0.3')
    print(f"Testing W=[0,{wstar}], Y=[{mp.nstr(Ylo,15)},{mp.nstr(Yhi,15)}]")
    result = certify(0, wstar, mp.nstr(Ylo,15), mp.nstr(Yhi,15), wstar)
    print(f"Result: {result}")
    print(f"Elapsed: {time.time()-t0:.2f}s")
