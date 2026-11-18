# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
# Sign certificate for the 04_hohmann_endpoint region: t in [0, 1e-4],
# fixed box a in [-0.25, 0.25]. Certifies S1 > 0, S2 > 0, t*Cpp > 0,
# opposite signs of Ghat at the a-box boundaries (existence+uniqueness
# of the stationary root), and t*H > 0 (via a subdivision search for a
# comfortable margin above a threshold of 1.0).
import sympy as sp
import mpmath as mp
import pickle

mp.iv.dps = 50
mp.mp.dps = 60

with open('/tmp/hohmann_regular_exprs.pkl','rb') as f:
    d = pickle.load(f)

t, a = sp.symbols('t a')
Ghat = d['Ghat']; tCpphat = d['tCpphat']; tHhat = d['tHhat']
S1_hat = d['S1_hat']; S2_hat = d['S2_hat']

f_Ghat = sp.lambdify((t,a), Ghat, mp.iv)
f_tCpphat = sp.lambdify((t,a), tCpphat, mp.iv)
f_tHhat = sp.lambdify((t,a), tHhat, mp.iv)
f_S1hat = sp.lambdify((t,a), S1_hat, mp.iv)
f_S2hat = sp.lambdify((t,a), S2_hat, mp.iv)

t0 = mp.mpf('1e-4')
a_lo, a_hi = mp.mpf('-0.25'), mp.mpf('0.25')
Bt_full = mp.iv.mpf([0, t0])

print("Step 1+2: S1,S2>0 on the full box")
S1_box = f_S1hat(Bt_full, mp.iv.mpf([a_lo,a_hi]))
S2_box = f_S2hat(Bt_full, mp.iv.mpf([a_lo,a_hi]))
print("S1:", S1_box, "> 0:", mp.mpf(S1_box.a)>0)
print("S2:", S2_box, "> 0:", mp.mpf(S2_box.a)>0)

print("\nStep 3: t*Cpp>0 on the full box")
tCpp_box = f_tCpphat(Bt_full, mp.iv.mpf([a_lo,a_hi]))
print(tCpp_box, " >0:", mp.mpf(tCpp_box.a)>0)

print("\nStep 4: opposite signs at the edges")
Gm = f_Ghat(Bt_full, mp.iv.mpf([a_lo,a_lo]))
Gp = f_Ghat(Bt_full, mp.iv.mpf([a_hi,a_hi]))
print("G(-0.25):", Gm, " <0:", mp.mpf(Gm.b)<0)
print("G(+0.25):", Gp, " >0:", mp.mpf(Gp.a)>0)

print("\nStep 5: t*H>1.0 -- full box first, then subdivided boxes")
box_full = f_tHhat(Bt_full, mp.iv.mpf([a_lo,a_hi]))
print("full box:", box_full, " lower>1.0:", mp.mpf(box_full.a)>1)

def min_lower_bound(n_a, n_t, thresh):
    worst = None
    fails = 0
    for j in range(n_t):
        tl = t0*j/n_t; th = t0*(j+1)/n_t
        Bt_sub = mp.iv.mpf([tl, th])
        for i in range(n_a):
            al = a_lo + (a_hi-a_lo)*i/n_a
            ah = a_lo + (a_hi-a_lo)*(i+1)/n_a
            box = f_tHhat(Bt_sub, mp.iv.mpf([al, ah]))
            lo = mp.mpf(box.a)
            if lo <= thresh:
                fails += 1
            if worst is None or lo < worst[0]:
                worst = (lo, al, ah, tl, th)
    return worst, fails

for n_a, n_t in [(50,1), (100,1), (200,1)]:
    w, fails = min_lower_bound(n_a, n_t, mp.mpf('1.0'))
    print(f"n_a={n_a} n_t={n_t}: worst={mp.nstr(w[0],12)} at a in [{mp.nstr(w[1],6)},{mp.nstr(w[2],6)}]  fails(<=1.0): {fails}/{n_a*n_t}")
