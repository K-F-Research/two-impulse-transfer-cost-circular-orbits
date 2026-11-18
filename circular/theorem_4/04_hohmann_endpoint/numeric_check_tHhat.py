# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
"""Numerical check of tHhat(t,a) near the Hohmann endpoint (Theorem 4, R=2)."""
import sympy as sp
import mpmath as mp
import random

mp.mp.dps = 60

p, x, t, a = sp.symbols('p x t a', real=True)
u1, u2 = sp.symbols('u1 u2', positive=True)

A0 = p - 1
B0 = p/2 - 1
N = A0*x + B0
E = A0**2 + N**2/(1-x**2)

S1x = sp.cancel(sp.together(3 + (E-1)/p - 2*sp.sqrt(p)))
S2x = sp.cancel(sp.together(sp.Rational(3,2) + (E-1)/p - sp.sqrt(p/2)))

def _cdiff(expr, *vs):
    return sp.cancel(sp.together(sp.diff(expr, *vs)))

S1t = S1x.subs(x, 1-2*t)
S2t = S2x.subs(x, 1-2*t)

pH = sp.Rational(4, 3)
p_of_t_a = pH + a*t
def reg(expr):
    return sp.cancel(expr.subs(p, p_of_t_a))

S1h = reg(S1t); S2h = reg(S2t)
S1p_h = reg(_cdiff(S1t, p)); S2p_h = reg(_cdiff(S2t, p))
S1t_h = reg(_cdiff(S1t, t)); S2t_h = reg(_cdiff(S2t, t))
S1pp_h = reg(_cdiff(S1t, p, 2)); S2pp_h = reg(_cdiff(S2t, p, 2))
S1pt_h = reg(_cdiff(_cdiff(S1t, p), t)); S2pt_h = reg(_cdiff(_cdiff(S2t, p), t))
S1tt_h = reg(_cdiff(S1t, t, 2)); S2tt_h = reg(_cdiff(S2t, t, 2))

Ghat = S1p_h/(2*u1) + S2p_h/(2*u2)
Cthat = S1t_h/(2*u1) + S2t_h/(2*u2)
tCpphat = (reg(t*S1pp_h)/(2*u1) - t*S1p_h**2/(4*S1h*u1)
         + reg(t*S2pp_h)/(2*u2) - t*S2p_h**2/(4*S2h*u2))
tCpthat = (reg(t*S1pt_h)/(2*u1) - t*S1p_h*S1t_h/(4*S1h*u1)
         + reg(t*S2pt_h)/(2*u2) - t*S2p_h*S2t_h/(4*S2h*u2))
tCtthat = (reg(t*S1tt_h)/(2*u1) - t*S1t_h**2/(4*S1h*u1)
         + reg(t*S2tt_h)/(2*u2) - t*S2t_h**2/(4*S2h*u2))

tHhat_constructed = (1-t)*tCtthat*tCpphat + (1-2*t)*Cthat*tCpphat/2 - (1-t)*tCpthat**2

# Direct H(pH+at,t), recomputed via une route TOTALEMENT independante :
# on differencie S1t,S2t par rapport a (p,t) EN PREMIER, avec p libre
# (exactement la discipline "derivees partielles avant substitution sur la
# branche" -- substituer p=pH+a*t avant de differencier par rapport a p
# donnerait trivialement 0, bug deja piege et corrige ici), PUIS seulement
# on substitue p=pH+a*t dans les derivees deja calculees.
Ct_direct_free = sp.diff(S1t, t)/(2*sp.sqrt(S1t)) + sp.diff(S2t, t)/(2*sp.sqrt(S2t))

def Csecond_free(S1, S2, q, r):
    res = 0
    for S in (S1, S2):
        Sq = sp.diff(S, q); Sr = sp.diff(S, r); Sqr = sp.diff(S, q, r)
        res += Sqr/(2*sp.sqrt(S)) - Sq*Sr/(4*S**sp.Rational(3,2))
    return res

Cpp_direct_free = Csecond_free(S1t, S2t, p, p)
Ctt_direct_free = Csecond_free(S1t, S2t, t, t)
Cpt_direct_free = Csecond_free(S1t, S2t, p, t)

# Substitution p=pH+a*t APRES differentiation (seule la substitution est
# faite ici, pas de nouvelle derivation).
Ct_direct_expr = Ct_direct_free.subs(p, p_of_t_a)
Cpp_direct_expr = Cpp_direct_free.subs(p, p_of_t_a)
Ctt_direct_expr = Ctt_direct_free.subs(p, p_of_t_a)
Cpt_direct_expr = Cpt_direct_free.subs(p, p_of_t_a)

# lambdify numerically (mpmath), evaluate pointwise -- never call cancel/
# simplify on these, so no OOM risk; sp.diff alone (no simplification) on
# a MODERATE expression evaluated numerically is cheap.
f_Ct = sp.lambdify((t, a), Ct_direct_expr, 'mpmath')
f_Cpp = sp.lambdify((t, a), Cpp_direct_expr, 'mpmath')
f_Ctt = sp.lambdify((t, a), Ctt_direct_expr, 'mpmath')
f_Cpt = sp.lambdify((t, a), Cpt_direct_expr, 'mpmath')

f_tHhat = sp.lambdify((t, a, u1, u2), tHhat_constructed, 'mpmath')
f_S1h = sp.lambdify((t, a), S1h, 'mpmath')
f_S2h = sp.lambdify((t, a), S2h, 'mpmath')

random.seed(12345)
n_points = 40
max_rel_err = mp.mpf(0)
n_ok = 0
n_tested = 0

print("Controle numerique haute precision (dps=60) : tHhat_construit vs t*H direct")
print("(route INDEPENDANTE : H direct differencie S1t,S2t substituees en p=pH+at,")
print(" sans passer par S1h/S2h/tCpphat/tCpthat/tCtthat deja construits)\n")

attempts = 0
while n_tested < n_points and attempts < 2000:
    attempts += 1
    tv = mp.mpf(random.uniform(0.001, 0.49))
    av = mp.mpf(random.uniform(-0.3, 0.3))
    try:
        S1v = f_S1h(tv, av); S2v = f_S2h(tv, av)
        if S1v <= 0 or S2v <= 0:
            continue
        u1v = mp.sqrt(S1v); u2v = mp.sqrt(S2v)
        lhs = f_tHhat(tv, av, u1v, u2v)

        Ctv = f_Ct(tv, av); Cppv = f_Cpp(tv, av); Cttv = f_Ctt(tv, av); Cptv = f_Cpt(tv, av)
        Hv = (tv*(1-tv)*Cttv + (1-2*tv)/2*Ctv)*Cppv - tv*(1-tv)*Cptv**2
        rhs = tv*Hv

        if rhs == 0:
            continue
        rel_err = abs(lhs - rhs) / abs(rhs)
        n_tested += 1
        if rel_err > max_rel_err:
            max_rel_err = rel_err
        if rel_err < mp.mpf(10)**(-45):
            n_ok += 1
        else:
            print(f"  ECART au point t={tv}, a={av}: lhs={lhs}, rhs={rhs}, rel_err={rel_err}")
    except (ValueError, ZeroDivisionError):
        continue

print(f"\nPoints testes: {n_tested}  (tentatives: {attempts})")
print(f"Points avec accord a 1e-45 relatif ou mieux: {n_ok}/{n_tested}")
print(f"Pire ecart relatif observe: {max_rel_err}")
print(f"\nRESULTAT: {'PASS (controle numerique, pas une preuve symbolique)' if n_ok == n_tested else 'ECHEC -- DIVERGENCE DETECTEE'}")
