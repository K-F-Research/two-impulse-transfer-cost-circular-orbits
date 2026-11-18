"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Independent audit of the Hohmann endpoint certificate (Theorem 4, R=2).

This script does not reload any pickle or symbolic object from the main
construction script (hohmann_build_Rta.py): it rebuilds the regularized
expressions Ghat(t,a), tHhat(t,a), S1_hat, S2_hat from scratch, through its
own computation path, then redoes the interval evaluation over the 100
sub-boxes in `a`, at two different precisions (50 and 80 digits), as an
independent reproducibility check against the main script's own
implementation.
"""
import sympy as sp
import mpmath as mp

p, t, a = sp.symbols('p t a')
R = sp.Integer(2); r1 = sp.Integer(1); r2 = R
pH = sp.Rational(4,3)
cc = 1 - 2*t
s2 = 4*t*(1-t)
x = p/r1 - 1; y = p/r2 - 1
e2 = x**2 + (x*cc+y)**2/s2
vtheta1 = sp.sqrt(p)/r1; vtheta2 = sp.sqrt(p)/r2
vc1 = 1/sp.sqrt(r1); vc2 = 1/sp.sqrt(r2)
v1_2 = 2/r1 + (e2-1)/p
v2_2 = 2/r2 + (e2-1)/p
S1 = sp.together(v1_2 + vc1**2 - 2*vc1*vtheta1)
S2 = sp.together(v2_2 + vc2**2 - 2*vc2*vtheta2)

S1p = sp.diff(S1,p); S2p = sp.diff(S2,p)
S1t = sp.diff(S1,t); S2t = sp.diff(S2,t)
S1pp = sp.diff(S1p,p); S2pp = sp.diff(S2p,p)
S1pt = sp.diff(S1p,t); S2pt = sp.diff(S2p,t)
S1tt = sp.diff(S1t,t); S2tt = sp.diff(S2t,t)

p_of_t_a = pH + a*t
def reg(expr, mult=1):
    return sp.cancel((mult*expr).subs(p, p_of_t_a))

S1_hat = reg(S1); S2_hat = reg(S2)
S1p_hat = reg(S1p); S2p_hat = reg(S2p)
S1t_hat = reg(S1t); S2t_hat = reg(S2t)
tS1pp = reg(S1pp, t); tS2pp = reg(S2pp, t)
tS1pt = reg(S1pt, t); tS2pt = reg(S2pt, t)
tS1tt = reg(S1tt, t); tS2tt = reg(S2tt, t)

sqrtS1 = sp.sqrt(S1_hat); sqrtS2 = sp.sqrt(S2_hat)

Ghat = S1p_hat/(2*sqrtS1) + S2p_hat/(2*sqrtS2)
Cthat = S1t_hat/(2*sqrtS1) + S2t_hat/(2*sqrtS2)
tCpphat = tS1pp/(2*sqrtS1) - t*S1p_hat**2/(4*S1_hat**sp.Rational(3,2)) \
        + tS2pp/(2*sqrtS2) - t*S2p_hat**2/(4*S2_hat**sp.Rational(3,2))
tCpthat = tS1pt/(2*sqrtS1) - t*S1p_hat*S1t_hat/(4*S1_hat**sp.Rational(3,2)) \
        + tS2pt/(2*sqrtS2) - t*S2p_hat*S2t_hat/(4*S2_hat**sp.Rational(3,2))
tCtthat = tS1tt/(2*sqrtS1) - t*S1t_hat**2/(4*S1_hat**sp.Rational(3,2)) \
        + tS2tt/(2*sqrtS2) - t*S2t_hat**2/(4*S2_hat**sp.Rational(3,2))

tHhat = (1-t)*tCtthat*tCpphat + (1-2*t)/2*Cthat*tCpphat - (1-t)*tCpthat**2

print("Expressions reconstruites independamment (chemin de calcul propre, sans pickle).")

def run_certificate(dps_iv, dps_mp, n_a=100):
    mp.iv.dps = dps_iv
    mp.mp.dps = dps_mp
    f_Ghat = sp.lambdify((t,a), Ghat, mp.iv)
    f_tCpphat = sp.lambdify((t,a), tCpphat, mp.iv)
    f_tHhat = sp.lambdify((t,a), tHhat, mp.iv)
    f_S1hat = sp.lambdify((t,a), S1_hat, mp.iv)
    f_S2hat = sp.lambdify((t,a), S2_hat, mp.iv)

    t0 = mp.mpf('1e-4')
    a_lo, a_hi = mp.mpf('-0.25'), mp.mpf('0.25')
    Bt_full = mp.iv.mpf([0, t0])

    S1b = f_S1hat(Bt_full, mp.iv.mpf([a_lo,a_hi]))
    S2b = f_S2hat(Bt_full, mp.iv.mpf([a_lo,a_hi]))
    tCppb = f_tCpphat(Bt_full, mp.iv.mpf([a_lo,a_hi]))
    Gm = f_Ghat(Bt_full, mp.iv.mpf([a_lo,a_lo]))
    Gp = f_Ghat(Bt_full, mp.iv.mpf([a_hi,a_hi]))

    rows = []
    worst = None
    for i in range(n_a):
        al = a_lo + (a_hi-a_lo)*i/n_a
        ah = a_lo + (a_hi-a_lo)*(i+1)/n_a
        Ba_sub = mp.iv.mpf([al, ah])
        box = f_tHhat(Bt_full, Ba_sub)
        S1sub = f_S1hat(Bt_full, Ba_sub)
        S2sub = f_S2hat(Bt_full, Ba_sub)
        lo = mp.mpf(box.a); hi = mp.mpf(box.b)
        rows.append((al, ah, lo, hi, mp.mpf(S1sub.a), mp.mpf(S2sub.a)))
        if worst is None or lo < worst[0]:
            worst = (lo, al, ah)

    return {
        'dps_iv': dps_iv, 'dps_mp': dps_mp,
        'S1_full': S1b, 'S2_full': S2b, 'tCpp_full': tCppb,
        'G_minus': Gm, 'G_plus': Gp,
        'rows': rows, 'worst': worst,
    }

results = {}
for dps_iv, dps_mp in [(50,60), (80,90)]:
    print(f"\n=== Certificat independant a dps_iv={dps_iv}, dps_mp={dps_mp} ===")
    res = run_certificate(dps_iv, dps_mp)
    results[dps_iv] = res
    print(f"S1 (boite pleine) : {res['S1_full']}")
    print(f"S2 (boite pleine) : {res['S2_full']}")
    print(f"t*Cpp (boite pleine) : {res['tCpp_full']}")
    print(f"G(t,-0.25) : {res['G_minus']}   G(t,+0.25) : {res['G_plus']}")
    print(f"Pire sous-boite : a in [{mp.nstr(res['worst'][1],20)}, {mp.nstr(res['worst'][2],20)}]")
    print(f"  t*Hhat lower bound (pire cas) : {mp.nstr(res['worst'][0], 40)}")
    all_pos = all(r[4]>0 and r[5]>0 for r in res['rows'])
    all_above = all(r[2] > mp.mpf('1.0') for r in res['rows'])
    print(f"S1,S2>0 sur les 100 sous-boites : {all_pos}")
    print(f"t*Hhat > 1.0 sur les 100 sous-boites : {all_above}")

print("\n=== Comparaison de reproductibilite dps=50 vs dps=80 ===")
w50 = results[50]['worst'][0]
w80 = results[80]['worst'][0]
print("Pire borne a dps=50:", mp.nstr(w50, 30))
print("Pire borne a dps=80:", mp.nstr(w80, 30))
print("Difference:", mp.nstr(abs(w50-w80), 15))

# Write full table (100 rows) to a file with high precision
with open('/tmp/hohmann_audit_table.csv', 'w') as f:
    f.write("a_lo,a_hi,tHhat_lo,tHhat_hi,S1_lo,S2_lo\n")
    for r in results[80]['rows']:
        f.write(",".join(mp.nstr(x, 35) for x in r) + "\n")
print("\nTable complete (dps=80) ecrite dans /tmp/hohmann_audit_table.csv")
