# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
"""
crosscheck_Ppp_methods.py

Confronts, term by term, the TWO independent methods for computing the
total second derivative of P along the branch, at the point (w=0, Y=1/2):

  Method A (verify_Pprime_symbolic.py): direct substitution
    u1 = sqrt(S1(w,Y)), u2 = sqrt(S2(w,Y)) into P(w,Y,u1,u2) -> P(w,Y),
    then sp.diff directly in (w,Y). The dependence of u1,u2 on (w,Y) is
    therefore already absorbed BEFORE differentiation.

  Method B (build_Ppp_chainrule.py): FORMAL differentiation of
    P(w,Y,u1,u2) with respect to its 4 independent arguments (fast),
    then manual recomposition via the chain rule using the known
    derivatives of u1(w,Y)=sqrt(S1(w,Y)), u2(w,Y)=sqrt(S2(w,Y))
    (Q_ww, Q_wY, Q_YY, Q_Y -- the "total" second derivatives of
    P(w,Y,u1(w,Y),u2(w,Y)) with respect to (w,Y), but WITHOUT yet
    composing with Y(w) -- that last composition happens AFTER, in
    certify_Ppp_v2.py, via Y'(w), Y''(w)).

This script verifies that the two methods give EXACTLY the same
P_ww(0,1/2), P_wY(0,1/2), P_YY(0,1/2), P_Y(0,1/2) -- before even going as
far as the composition with Y(w). This is the requested consistency
check: confronting Q_ww/Q_wY/Q_YY/Q_Y (method B) against
P_ww/P_wY/P_YY/P_Y (method A), both evaluated at the same point
(w,Y,u1,u2) = (0, 1/2, sqrt(2), 1/sqrt(2)).
"""
import sympy as sp
import pickle, time

w, Y, u1, u2 = sp.symbols('w Y u1 u2', positive=True)

# ---------------------------------------------------------------------
# Method A: repeat the computation of verify_Pprime_symbolic.py (direct
# substitution before differentiation), keeping only what is needed for
# the check point (0,1/2) -- third derivatives are not needed here.
# ---------------------------------------------------------------------
with open('/tmp/P_built.pkl', 'rb') as f:
    P = pickle.load(f)

S1 = (8*Y**2*w**3 - 12*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 16*Y*w + 12*Y + 1) / (4*Y*(2 - Y*w**2))
S2 = (2*sp.sqrt(2)*Y**2*w**3 - 6*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 4*sp.sqrt(2)*Y*w + 1) / (4*Y*(2 - Y*w**2))

print("[Method A] direct substitution u1=sqrt(S1(w,Y)), u2=sqrt(S2(w,Y))...")
t0 = time.time()
P_wY = P.subs({u1: sp.sqrt(S1), u2: sp.sqrt(S2)}, simultaneous=True)
print(f"  [{time.time()-t0:.1f}s] count_ops(P_wY) = {sp.count_ops(P_wY)}")

P_w_A = sp.diff(P_wY, w)
P_Y_A = sp.diff(P_wY, Y)
P_ww_A = sp.diff(P_w_A, w)
P_wY_A = sp.diff(P_w_A, Y)
P_YY_A = sp.diff(P_Y_A, Y)

point = {w: 0, Y: sp.Rational(1, 2)}
Pw0_A = sp.simplify(P_w_A.subs(point))
PY0_A = sp.simplify(P_Y_A.subs(point))
Pww0_A = sp.simplify(P_ww_A.subs(point))
PwY0_A = sp.simplify(P_wY_A.subs(point))
PYY0_A = sp.simplify(P_YY_A.subs(point))

print(f"  P_w (0,1/2)  = {Pw0_A}")
print(f"  P_Y (0,1/2)  = {PY0_A}")
print(f"  P_ww(0,1/2)  = {Pww0_A}")
print(f"  P_wY(0,1/2)  = {PwY0_A}")
print(f"  P_YY(0,1/2)  = {PYY0_A}")

# ---------------------------------------------------------------------
# Method B: reload Q_ww, Q_wY, Q_YY, Q_Y (already built by
# build_Ppp_chainrule.py -- the TOTAL second derivatives of
# P(w,Y,u1(w,Y),u2(w,Y)) with respect to (w,Y), via formal chain rule),
# and evaluate them at the same numeric point by supplying u1,u2
# explicitly.
# ---------------------------------------------------------------------
print("\n[Method B] reloading Q_ww, Q_wY, Q_YY, Q_Y (formal chain rule)...")
with open('/tmp/Q_second_partials.pkl', 'rb') as f:
    d = pickle.load(f)
Q_ww, Q_wY, Q_YY, Q_Y = d['Q_ww'], d['Q_wY'], d['Q_YY'], d['Q_Y']

point_B = {w: 0, Y: sp.Rational(1, 2), u1: sp.sqrt(2), u2: 1/sp.sqrt(2)}
Qww0_B = sp.simplify(Q_ww.subs(point_B))
QwY0_B = sp.simplify(Q_wY.subs(point_B))
QYY0_B = sp.simplify(Q_YY.subs(point_B))
QY0_B = sp.simplify(Q_Y.subs(point_B))

print(f"  Q_ww(0,1/2,sqrt2,1/sqrt2) = {Qww0_B}")
print(f"  Q_wY(0,1/2,sqrt2,1/sqrt2) = {QwY0_B}")
print(f"  Q_YY(0,1/2,sqrt2,1/sqrt2) = {QYY0_B}")
print(f"  Q_Y (0,1/2,sqrt2,1/sqrt2) = {QY0_B}")

# ---------------------------------------------------------------------
# Term-by-term confrontation
# ---------------------------------------------------------------------
print("\n=== CONFRONTATION ===")
checks = [
    ("P_Y  (A)  vs  Q_Y  (B)", PY0_A, QY0_B),
    ("P_ww (A)  vs  Q_ww (B)", Pww0_A, Qww0_B),
    ("P_wY (A)  vs  Q_wY (B)", PwY0_A, QwY0_B),
    ("P_YY (A)  vs  Q_YY (B)", PYY0_A, QYY0_B),
]
all_ok = True
for label, a, b in checks:
    diff = sp.simplify(a - b)
    ok = (diff == 0)
    all_ok &= ok
    print(f"  {label}: A={a}  B={b}  diff={diff}  {'OK' if ok else 'MISMATCH'}")

print(f"\nAll quantities coincide term by term: {all_ok}")

# ---------------------------------------------------------------------
# If everything agrees, recompose P''(0) along the branch via method B
# (chain rule + Y', Y'') and compare against L0, cross-checked against
# method A (already done in verify_Pprime_symbolic.py: P''(0)/2=L0).
# ---------------------------------------------------------------------
b_branch = -sp.Rational(1, 3) - sp.sqrt(2)/6   # Y'(0)
c_branch = (2 + sp.sqrt(2))/4                   # Y(w)=1/2+b w + c w^2 + ...
Ypp0 = 2*c_branch                               # Y''(0) = 2c

Pprimeprime0_B = Qww0_B + 2*QwY0_B*b_branch + QYY0_B*b_branch**2 + QY0_B*Ypp0
Pprimeprime0_B_s = sp.simplify(Pprimeprime0_B)
L0 = -(9 + 4*sp.sqrt(2))/32
print(f"\nP''(0) recomposed via method B = {Pprimeprime0_B_s}")
print(f"P''(0)/2 (method B) = {sp.simplify(Pprimeprime0_B_s/2)}")
print(f"L0 = {L0}")
print(f"Does method B match L0? {sp.simplify(Pprimeprime0_B_s/2 - L0) == 0}")
