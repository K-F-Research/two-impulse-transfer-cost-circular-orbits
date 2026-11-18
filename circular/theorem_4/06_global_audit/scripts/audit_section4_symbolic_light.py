"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Memory-light symbolic verification of section 4 (tHhat = t*H). Instead of
calling sp.cancel() on the entire combined expression -- which requires
factoring/GCD of enormous multivariate polynomials and is memory-intensive
-- this script goes through sp.together() + sp.fraction() to isolate the
numerator, then sp.expand() alone (never factor/cancel on the whole
expression) to verify it is identically zero: expand() distributes without
ever computing a polynomial GCD, avoiding cancel()'s memory peak.
"""
import sympy as sp
import time, resource

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

def logmem(label):
    rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"[{label}] peak RSS so far: {rss_mb:.0f} MB", flush=True)

t0 = time.time()
S1h = reg(S1t); S2h = reg(S2t)
S1p_h = reg(_cdiff(S1t, p)); S2p_h = reg(_cdiff(S2t, p))
S1t_h = reg(_cdiff(S1t, t)); S2t_h = reg(_cdiff(S2t, t))
S1pp_h = reg(_cdiff(S1t, p, 2)); S2pp_h = reg(_cdiff(S2t, p, 2))
S1pt_h = reg(_cdiff(_cdiff(S1t, p), t)); S2pt_h = reg(_cdiff(_cdiff(S2t, p), t))
S1tt_h = reg(_cdiff(S1t, t, 2)); S2tt_h = reg(_cdiff(S2t, t, 2))
print(f"pieces: {time.time()-t0:.1f}s", flush=True); logmem("pieces"); t0=time.time()

Ghat = S1p_h/(2*u1) + S2p_h/(2*u2)
Cthat = S1t_h/(2*u1) + S2t_h/(2*u2)
tCpphat = (reg(t*S1pp_h)/(2*u1) - t*S1p_h**2/(4*S1h*u1)
         + reg(t*S2pp_h)/(2*u2) - t*S2p_h**2/(4*S2h*u2))
tCpthat = (reg(t*S1pt_h)/(2*u1) - t*S1p_h*S1t_h/(4*S1h*u1)
         + reg(t*S2pt_h)/(2*u2) - t*S2p_h*S2t_h/(4*S2h*u2))
tCtthat = (reg(t*S1tt_h)/(2*u1) - t*S1t_h**2/(4*S1h*u1)
         + reg(t*S2tt_h)/(2*u2) - t*S2t_h**2/(4*S2h*u2))
print(f"Ghat etc: {time.time()-t0:.1f}s", flush=True); logmem("Ghat"); t0=time.time()

tHhat_constructed = (1-t)*tCtthat*tCpphat + (1-2*t)*Cthat*tCpphat/2 - (1-t)*tCpthat**2

Cpp_hat_full = tCpphat/t
Cpt_hat_full = tCpthat/t
Ctt_hat_full = tCtthat/t
Cphiphi_hat = t*(1-t)*Ctt_hat_full + (1-2*t)/2*Cthat

# PAS de cancel() sur le produit complet -- on construit la difference
# brute, on la "together"-ise (mise au meme denominateur, sans factoriser),
# puis on isole le numerateur et on l'expand (distribue, pas de GCD).
diff_expr = tHhat_constructed - t*(Cphiphi_hat*Cpp_hat_full - t*(1-t)*Cpt_hat_full**2)
print(f"diff_expr built: {time.time()-t0:.1f}s, ops={sp.count_ops(diff_expr)}", flush=True)
logmem("diff_expr"); t0=time.time()

tog = sp.together(diff_expr)
print(f"together: {time.time()-t0:.1f}s, ops={sp.count_ops(tog)}", flush=True)
logmem("together"); t0=time.time()

num, den = sp.fraction(tog)
print(f"fraction split: {time.time()-t0:.1f}s, num_ops={sp.count_ops(num)}", flush=True)
logmem("fraction"); t0=time.time()

num_expanded = sp.expand(num)
print(f"expand numerator: {time.time()-t0:.1f}s, ops={sp.count_ops(num_expanded)}", flush=True)
logmem("expand"); t0=time.time()

print(f"\nRESULTAT: numerateur == 0 : {num_expanded == 0}", flush=True)
if num_expanded != 0:
    print("Residu (apres expand, PAS factorise) -- a examiner :")
    print(num_expanded)
