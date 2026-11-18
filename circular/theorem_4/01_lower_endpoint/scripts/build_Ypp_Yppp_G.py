# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
import sympy as sp

p, x = sp.symbols('p x')
w, Y = sp.symbols('w Y', positive=True)

x1 = p - 1
x2 = p/2 - 1
s2 = 1 - x**2
e2 = x1**2 + (x1*x + x2)**2 / s2
v1sq = 2 + (e2 - 1)/p
v2sq = 1 + (e2 - 1)/p
S1_orig = v1sq + 1 - 2*sp.sqrt(p)
S2_orig = v2sq + sp.Rational(1,2) - sp.sqrt(p/2)
S1_p_orig = sp.diff(S1_orig, p)
S2_p_orig = sp.diff(S2_orig, p)
subs = {p: w**2, x: -1 + w**2*Y}
S1_p_reg = sp.simplify(S1_p_orig.subs(subs, simultaneous=True))
S2_p_reg = sp.simplify(S2_p_orig.subs(subs, simultaneous=True))

S1 = (8*Y**2*w**3 - 12*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 16*Y*w + 12*Y + 1) / (4*Y*(2 - Y*w**2))
S2 = (2*sp.sqrt(2)*Y**2*w**3 - 6*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 4*sp.sqrt(2)*Y*w + 1) / (4*Y*(2 - Y*w**2))
Cp_reg = S1_p_reg/(2*sp.sqrt(S1)) + S2_p_reg/(2*sp.sqrt(S2))
G = sp.together(w**2 * Cp_reg)

print("Building G derivatives up to order 3 (needed for Y''', via implicit function thm)...")
G_w = sp.diff(G, w); G_Y = sp.diff(G, Y)
G_ww = sp.diff(G, w, w); G_wY = sp.diff(G, w, Y); G_YY = sp.diff(G, Y, Y)
G_www = sp.diff(G, w, w, w); G_wwY = sp.diff(G, w, w, Y); G_wYY = sp.diff(G, w, Y, Y); G_YYY = sp.diff(G, Y, Y, Y)

point = {w: 0, Y: sp.Rational(1,2)}
vals = {}
for name, expr in [("G_w",G_w),("G_Y",G_Y),("G_ww",G_ww),("G_wY",G_wY),("G_YY",G_YY),
                    ("G_www",G_www),("G_wwY",G_wwY),("G_wYY",G_wYY),("G_YYY",G_YYY)]:
    vals[name] = sp.simplify(expr.subs(point))
    print(f"{name}(0,1/2) = {vals[name]}")

Yp = -vals["G_w"]/vals["G_Y"]
print("\nY'(0) =", sp.simplify(Yp))

Ypp = -(vals["G_ww"] + 2*vals["G_wY"]*Yp + vals["G_YY"]*Yp**2)/vals["G_Y"]
Ypp = sp.simplify(Ypp)
print("Y''(0) =", Ypp, " vs 2c known =", sp.simplify(2*(2+sp.sqrt(2))/4))

Yppp = -(vals["G_www"] + 3*vals["G_wwY"]*Yp + 3*vals["G_wYY"]*Yp**2 + vals["G_YYY"]*Yp**3
         + 3*(vals["G_wY"]+vals["G_YY"]*Yp)*Ypp) / vals["G_Y"]
Yppp = sp.simplify(Yppp)
print("Y'''(0) =", Yppp, " numeric:", sp.N(Yppp,15))
