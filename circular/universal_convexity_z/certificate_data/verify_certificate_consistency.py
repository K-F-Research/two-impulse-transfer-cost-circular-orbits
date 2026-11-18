# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
# Independent consistency check of the raw polynomial certificate data
# behind Theorem 3 (circular/universal_convexity_z/), Section 4. See
# certificate_data/README.md for what this does and does not establish.
#
# This script loads the five polynomial objects RH, RJ, G, A, B (saved as
# exact-rational multivariate polynomials in (t, w, rho) via python-flint)
# and checks the algebraic identities that the original derivation's final
# step depends on:
#
#   A * G == RH         (G divides RH exactly, quotient is A)
#   B * G == RJ         (G divides RJ exactly, quotient is B)
#   gcd(A, B) == 1      (A and B are coprime -- the fact used by the
#                        topological continuity argument of Sec. 4)
#   gcd(RH, RJ) == G    (G is genuinely the gcd of RH and RJ, not merely a
#                        common factor asserted without computing the gcd)
#
# All four checks are re-computed here, not re-read from a stored result.

import pickle
import flint

ctx3 = flint.fmpq_mpoly_ctx.get(('t', 'w', 'rho'), 'lex')


def load(name):
    with open(name, 'rb') as f:
        d = pickle.load(f)
    return ctx3.from_dict(d)


def main():
    RH = load('RH_exact_flint.pkl')
    RJ = load('RJ_exact_flint.pkl')
    G = load('G_exact_flint.pkl')
    A = load('A_exact_flint.pkl')
    B = load('B_exact_flint.pkl')

    print("RH: terms =", len(list(RH.terms())), " total_degree =", RH.total_degree())
    print("RJ: terms =", len(list(RJ.terms())), " total_degree =", RJ.total_degree())
    print("G:  terms =", len(list(G.terms())), " total_degree =", G.total_degree())
    print("A:  terms =", len(list(A.terms())), " total_degree =", A.total_degree())
    print("B:  terms =", len(list(B.terms())), " total_degree =", B.total_degree())

    ok1 = (A * G == RH)
    print("A*G == RH:", ok1)

    ok2 = (B * G == RJ)
    print("B*G == RJ:", ok2)

    gAB = A.gcd(B)
    ok3 = (gAB.total_degree() == 0)
    print("gcd(A,B) is a nonzero constant (A,B coprime):", ok3,
          " (gcd has", len(list(gAB.terms())), "term(s))")

    gRHRJ = RH.gcd(RJ)
    ok4 = (gRHRJ == G)
    print("gcd(RH,RJ) == G (recomputed, not merely re-read):", ok4)

    assert ok1 and ok2 and ok3 and ok4, "certificate consistency check FAILED"
    print("\nAll four consistency checks passed.")


if __name__ == '__main__':
    main()
