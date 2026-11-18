"""
Copyright (c) 2026 - Frederic G. Speyser & Kseniya Vyatkina

Note:
Continues the exact reduction chain: builds J from H and Q, computes RJ = Res_rho(Q, J),
factors out the common part G = gcd(RH, RJ), and verifies that the resulting A = RH/G
and B = RJ/G are coprime -- the algebraic fact underlying the sign-constancy argument
used to prove strict convexity of V(z) for every R > 1.
"""

import time, re, pickle
import flint

def parse_poly_t_w_rho(path):
    with open(path) as f:
        s = f.read().strip()
    if s[0] not in '+-':
        s = '+' + s
    terms = {}
    raw_terms = re.findall(r'[+-][^+-]+', s)
    for rt in raw_terms:
        sign = 1 if rt[0] == '+' else -1
        body = rt[1:]
        a=b=c=0
        coef=None
        for factor in body.split('*'):
            factor=factor.strip()
            if factor=='':continue
            if factor.startswith('t'):
                a = int(factor.split('^')[1]) if '^' in factor else 1
            elif factor.startswith('w'):
                b = int(factor.split('^')[1]) if '^' in factor else 1
            elif factor.startswith('rho'):
                rest=factor[3:]
                c = int(rest[1:]) if rest.startswith('^') else 1
            else:
                coef=int(factor)
        if coef is None: coef=1
        coef*=sign
        key=(a,b,c)
        terms[key]=terms.get(key,0)+coef
    return terms

Qterms = parse_poly_t_w_rho('Q_poly.txt')
Hterms = parse_poly_t_w_rho('H_poly.txt')

ctx3 = flint.fmpq_mpoly_ctx.get(('t','w','rho'), 'lex')
t_, w_, rho_ = ctx3.gens()

def to_mpoly3(terms, ctx):
    d = {k:flint.fmpq(v) for k,v in terms.items()}
    return ctx.from_dict(d)

print("building Q,H...", flush=True)
Q3 = to_mpoly3(Qterms, ctx3)
H3 = to_mpoly3(Hterms, ctx3)

print("computing J = Hw*Qr - Hr*Qw ...", flush=True)
start = time.time()
Hw = H3.derivative(1)  # derivative wrt w (index 1: t=0,w=1,rho=2)
Hr = H3.derivative(2)
Qw = Q3.derivative(1)
Qr = Q3.derivative(2)
J3 = Hw*Qr - Hr*Qw
print("J built, time:", time.time()-start, " terms:", len(list(J3.terms())), flush=True)

print("loading RH (already computed)...", flush=True)
with open('RH_exact_flint.pkl','rb') as f:
    RH_dict = pickle.load(f)
pass
RH = ctx3.from_dict(RH_dict)
print("RH loaded, terms:", len(list(RH.terms())), flush=True)

print("computing RJ = Res_rho(Q,J) ...", flush=True)
start = time.time()
RJ = Q3.resultant(J3, 'rho')
print("RJ time:", time.time()-start, " terms:", len(list(RJ.terms())), " total_degree:", RJ.total_degree(), flush=True)
with open('RJ_exact_flint.pkl','wb') as f:
    pickle.dump(RJ.to_dict(), f)

print("computing G = gcd(RH,RJ) ...", flush=True)
start = time.time()
G = RH.gcd(RJ)
print("G time:", time.time()-start, " terms:", len(list(G.terms())), " total_degree:", G.total_degree(), flush=True)
with open('G_exact_flint.pkl','wb') as f:
    pickle.dump(G.to_dict(), f)

print("computing A=RH/G, B=RJ/G ...", flush=True)
start = time.time()
A = RH.exact_division(G) if hasattr(RH,'exact_division') else RH / G
B = RJ.exact_division(G) if hasattr(RJ,'exact_division') else RJ / G
print("division time:", time.time()-start, flush=True)
print("A terms:", len(list(A.terms())), " B terms:", len(list(B.terms())), flush=True)
with open('A_exact_flint.pkl','wb') as f:
    pickle.dump(A.to_dict(), f)
with open('B_exact_flint.pkl','wb') as f:
    pickle.dump(B.to_dict(), f)

print("checking gcd(A,B)...", flush=True)
start=time.time()
gAB = A.gcd(B)
print("gcd(A,B) time:", time.time()-start, " terms:", len(list(gAB.terms())), flush=True)

print("ALL DONE CHAIN (before final big resultant)", flush=True)
