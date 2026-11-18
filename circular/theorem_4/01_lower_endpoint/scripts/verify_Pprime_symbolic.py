# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
import sympy as sp
import time, pickle

w, Y = sp.symbols('w Y', positive=True)
with open('/tmp/P_built.pkl','rb') as f:
    P = pickle.load(f)
u1, u2 = sp.symbols('u1 u2', positive=True)

S1 = (8*Y**2*w**3 - 12*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 16*Y*w + 12*Y + 1) / (4*Y*(2 - Y*w**2))
S2 = (2*sp.sqrt(2)*Y**2*w**3 - 6*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 4*sp.sqrt(2)*Y*w + 1) / (4*Y*(2 - Y*w**2))

print("Substituting u1=sqrt(S1), u2=sqrt(S2) into P (keeping w,Y symbolic)...")
t0=time.time()
P_wY = P.subs({u1: sp.sqrt(S1), u2: sp.sqrt(S2)}, simultaneous=True)
print(f"[{time.time()-t0:.1f}s] count_ops(P_wY) = {sp.count_ops(P_wY)}")

print("Computing partial derivatives dP/dw, dP/dY (symbolic, no further simplify)...")
t0=time.time()
P_w = sp.diff(P_wY, w)
print(f"[{time.time()-t0:.1f}s] P_w done, count_ops={sp.count_ops(P_w)}")
t0=time.time()
P_Y = sp.diff(P_wY, Y)
print(f"[{time.time()-t0:.1f}s] P_Y done, count_ops={sp.count_ops(P_Y)}")

b = -sp.Rational(1,3) - sp.sqrt(2)/6
point = {w: 0, Y: sp.Rational(1,2)}

t0=time.time()
Pw0 = sp.simplify(P_w.subs(point))
print(f"[{time.time()-t0:.1f}s] P_w(0,1/2) = {Pw0}")
t0=time.time()
PY0 = sp.simplify(P_Y.subs(point))
print(f"[{time.time()-t0:.1f}s] P_Y(0,1/2) = {PY0}")

Pprime0 = Pw0 + PY0*b
print("\nP'(0) along branch = P_w(0,1/2) + P_Y(0,1/2)*b =", sp.simplify(Pprime0))

print("\n--- Second derivative ---")
t0=time.time()
P_ww = sp.diff(P_w, w)
print(f"[{time.time()-t0:.1f}s] P_ww done")
t0=time.time()
P_wY2 = sp.diff(P_w, Y)
print(f"[{time.time()-t0:.1f}s] P_wY done")
t0=time.time()
P_YY = sp.diff(P_Y, Y)
print(f"[{time.time()-t0:.1f}s] P_YY done")

Pww0 = sp.simplify(P_ww.subs(point))
PwY0 = sp.simplify(P_wY2.subs(point))
PYY0 = sp.simplify(P_YY.subs(point))
print("\nP_ww(0,1/2) =", Pww0)
print("P_wY(0,1/2) =", PwY0)
print("P_YY(0,1/2) =", PYY0)

c = (2+sp.sqrt(2))/4
Y2prime0 = 2*c  # Y''(0) = 2*c since Y(w)=1/2+b*w+c*w^2+...

Pprimeprime0 = Pww0 + 2*PwY0*b + PYY0*b**2 + PY0*Y2prime0
Pprimeprime0_s = sp.simplify(Pprimeprime0)
print("\nP''(0) along branch =", Pprimeprime0_s)
print("P''(0)/2 =", sp.simplify(Pprimeprime0_s/2))

L0 = -(9+4*sp.sqrt(2))/32
print("L0 =", L0, " numeric:", sp.N(L0,20))
print("Match P''(0)/2 == L0 ?", sp.simplify(Pprimeprime0_s/2 - L0) == 0)

print("\n--- Third derivatives (for P'''(0) point value as first estimate) ---")
t0=time.time()
P_www = sp.diff(P_ww, w)
print(f"[{time.time()-t0:.1f}s] P_www done")
t0=time.time()
P_wwY = sp.diff(P_ww, Y)
print(f"[{time.time()-t0:.1f}s] P_wwY done")
t0=time.time()
P_wYY = sp.diff(P_wY2, Y)
print(f"[{time.time()-t0:.1f}s] P_wYY done")
t0=time.time()
P_YYY = sp.diff(P_YY, Y)
print(f"[{time.time()-t0:.1f}s] P_YYY done")

Pwww0 = sp.simplify(P_www.subs(point))
PwwY0 = sp.simplify(P_wwY.subs(point))
PwYY0 = sp.simplify(P_wYY.subs(point))
PYYY0 = sp.simplify(P_YYY.subs(point))
print("\nP_www(0,1/2)=",Pwww0," P_wwY(0,1/2)=",PwwY0," P_wYY(0,1/2)=",PwYY0," P_YYY(0,1/2)=",PYYY0)

Yp3 = -(sp.Rational(59,36) + 131*sp.sqrt(2)/144)  # Y'''(0) just computed
Ypp2 = sp.sqrt(2)/2 + 1  # Y''(0)
bb = -sp.Rational(1,3) - sp.sqrt(2)/6  # Y'(0)

Ppppp0 = (Pwww0 + 3*PwwY0*bb + 3*PwYY0*bb**2 + PYYY0*bb**3
          + 3*PwY0*Ypp2 + 3*PYY0*bb*Ypp2 + PY0*Yp3)
Ppppp0_s = sp.simplify(Ppppp0)
print("\nP'''(0) along branch =", Ppppp0_s, " numeric:", sp.N(Ppppp0_s,15))
