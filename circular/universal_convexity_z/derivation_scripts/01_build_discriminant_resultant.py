"""
Copyright (c) 2026 - Frederic G. Speyser & Kseniya Vyatkina

Note:
Builds the exact three-variable polynomials Q and H from their stored term lists and
computes RH = Res_rho(Q, H) over the rationals using FLINT, the first link in the
exact algebraic chain behind the universal-convexity theorem (variable z).
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

def to_mpoly3(terms, ctx):
    d = {k:flint.fmpq(v) for k,v in terms.items()}
    return ctx.from_dict(d)

print("building...", flush=True)
Q3 = to_mpoly3(Qterms, ctx3)
H3 = to_mpoly3(Hterms, ctx3)
print("built, computing resultant...", flush=True)

start2=time.time()
RH = Q3.resultant(H3, 'rho')
print("FULL symbolic resultant time:", time.time()-start2, flush=True)
print("RH num terms:", len(list(RH.terms())), flush=True)
print("RH total degree:", RH.total_degree(), flush=True)

with open('RH_exact_flint.pkl','wb') as f:
    pickle.dump(RH.to_dict(), f)
print("saved. ALL DONE", flush=True)
