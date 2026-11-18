# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
import sympy as sp
import pickle

w, Y, u1, u2 = sp.symbols('w Y u1 u2', positive=True)
with open('/tmp/F0_dummy_cancelled.pkl','rb') as f:
    F0_cancelled = pickle.load(f)

num, den = sp.fraction(F0_cancelled)
print("den =", sp.factor(den))

# P := w^2 * F0_cancelled -- should cancel the w^2 in the denominator exactly
P = sp.cancel(w**2 * F0_cancelled)
num_P, den_P = sp.fraction(P)
print("\nden_P (should have NO w factor left) =", sp.factor(den_P))
print("count_ops(P) =", sp.count_ops(P))

with open('/tmp/P_built.pkl','wb') as f:
    pickle.dump(P, f)

print("\n--- Verify P(0, Y) at Y=1/2, u1=sqrt(2), u2=sqrt(1/2) ---")
P_at_w0 = P.subs(w, 0)
P_at_w0_Y05 = P_at_w0.subs({Y: sp.Rational(1,2), u1: sp.sqrt(2), u2: sp.sqrt(sp.Rational(1,2))})
print("P(0, 1/2) =", sp.simplify(P_at_w0_Y05))

print("\n--- P(0, Y) as a function of Y alone (with u1=sqrt(S1(0,Y)), u2=sqrt(S2(0,Y))) ---")
S1 = (8*Y**2*w**3 - 12*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 16*Y*w + 12*Y + 1) / (4*Y*(2 - Y*w**2))
S2 = (2*sp.sqrt(2)*Y**2*w**3 - 6*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 4*sp.sqrt(2)*Y*w + 1) / (4*Y*(2 - Y*w**2))
S1_at0 = S1.subs(w,0)
S2_at0 = S2.subs(w,0)
print("S1(0,Y) =", S1_at0, " S2(0,Y) =", S2_at0)
P_at_w0_genY = P_at_w0.subs({u1: sp.sqrt(S1_at0), u2: sp.sqrt(S2_at0)})
P_at_w0_genY_s = sp.simplify(P_at_w0_genY)
print("P(0,Y) simplified =", P_at_w0_genY_s)
# The claim used in lower_endpoint_sign_conclusion.md is P(0,1/2)=0,
# printed above, not an enumeration of every root of P(0,.) in Y.
