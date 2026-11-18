# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
import sympy as sp
import pickle, time

w, Y, u1, u2 = sp.symbols('w Y u1 u2', positive=True)
with open('/tmp/P_built.pkl','rb') as f:
    P = pickle.load(f)

t0=time.time()
P_at_w0 = P.subs(w, 0)
print(f"[{time.time()-t0:.1f}s] P(0,Y,u1,u2) count_ops={sp.count_ops(P_at_w0)}")
print(P_at_w0)

print("\nSubstituting Y=1/2, u1=sqrt(2), u2=sqrt(1/2):")
val = P_at_w0.subs({Y: sp.Rational(1,2), u1: sp.sqrt(2), u2: sp.sqrt(sp.Rational(1,2))})
print("raw:", val)
print("simplified:", sp.nsimplify(sp.simplify(val)))
print("numeric:", sp.N(val, 20))
