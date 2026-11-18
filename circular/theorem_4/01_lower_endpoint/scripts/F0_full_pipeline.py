# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
import sympy as sp
import mpmath as mp
import time

p, x, w, Y = sp.symbols('p x w Y', positive=True)

A = p - 1
B = p/2 - 1
N = A*x + B
E = A**2 + N**2/(1 - x**2)
S1 = 3 + (E - 1)/p - 2*sp.sqrt(p)
S2 = sp.Rational(3, 2) + (E - 1)/p - sp.sqrt(p/2)
S = [S1, S2]

def C_derivative(q, r=None):
    result = 0
    for Si in S:
        if r is None:
            result += sp.diff(Si, q)/(2*sp.sqrt(Si))
        else:
            Si_qr = sp.diff(Si, q, r)
            Si_q = sp.diff(Si, q)
            Si_r = sp.diff(Si, r)
            result += (Si_qr/(2*sp.sqrt(Si)) - Si_q*Si_r/(4*Si**sp.Rational(3, 2)))
    return result

t0=time.time()
Cp  = C_derivative(p)
Cpp = C_derivative(p, p)
Cx  = C_derivative(x)
Cxx = C_derivative(x, x)
Cpx = C_derivative(p, x)
H = ((1-x**2)*Cxx - x*Cx)*Cpp - (1-x**2)*Cpx**2
print(f"[{time.time()-t0:.1f}s] H built (pre-substitution)")

sub = {p: w**2, x: -1 + w**2*Y}
F0_raw = w**4 * H.subs(sub, simultaneous=True)
F0 = sp.together(F0_raw)
print(f"[{time.time()-t0:.1f}s] F0 substituted+together, count_ops={sp.count_ops(F0)}")

# --- lambdify directly, no factor/cancel (not needed for numeric/interval eval) ---
mp.mp.dps = 50
mp.iv.dps = 50
modules_mp = [{'sqrt': mp.sqrt, 'mpf': mp.mpf}, 'mpmath']
modules_iv = [{'sqrt': mp.iv.sqrt, 'mpf': mp.iv.mpf}, 'mpmath']

print("Lambdifying F0 (point + interval)...")
t1=time.time()
F0_f_mp = sp.lambdify((w,Y), F0, modules=modules_mp)
F0_f_iv = sp.lambdify((w,Y), F0, modules=modules_iv)
print(f"[{time.time()-t1:.1f}s] lambdify done")

# Test B: limit at w->0 (use tiny w, Y=0.5)
L0 = -(9+4*mp.sqrt(2))/32
print("\n--- Test B: F0(w->0, Y=0.5) vs L0 ---")
for wv in [mp.mpf('1e-3'), mp.mpf('1e-6'), mp.mpf('1e-9')]:
    val = F0_f_mp(wv, mp.mpf('0.5'))
    print(f"w={wv}: F0 = {mp.nstr(val,20)}")
print("L0 exact =", mp.nstr(L0,20))

print("\n--- Test B corrected: evaluate F0 at the TRUE branch Y(w), not fixed Y=0.5 ---")
b = -mp.mpf(1)/3 - mp.sqrt(2)/6
c = (2+mp.sqrt(2))/4
for wv in [mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-5'), mp.mpf('1e-6')]:
    Yv = mp.mpf('0.5') + b*wv + c*wv**2
    val = F0_f_mp(wv, Yv)
    print(f"w={mp.nstr(wv,6)}  Y(w)={mp.nstr(Yv,15)}  F0={mp.nstr(val,20)}")
print("L0 exact =", mp.nstr(L0,20))

print("\n--- Test A: identity check (substituted S1,S2 match previously validated regularized forms) ---")
S1_confirmed = (8*Y**2*w**3 - 12*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 16*Y*w + 12*Y + 1) / (4*Y*(2 - Y*w**2))
S2_confirmed = (2*sp.sqrt(2)*Y**2*w**3 - 6*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 4*sp.sqrt(2)*Y*w + 1) / (4*Y*(2 - Y*w**2))
S1_sub = sp.simplify(S1.subs(sub))
S2_sub = sp.simplify(S2.subs(sub))
print("S1 diff:", sp.simplify(S1_sub - S1_confirmed))
print("S2 diff:", sp.simplify(S2_sub - S2_confirmed))

print("\n--- Test C: cross-check vs fully independent H (F0_veryhighprec.py route), same w ---")
# F0_veryhighprec used symbolic diff of C(p,x) directly (different code path), branch solved via bisection
# Compare at w=1e-6 equivalent eps
