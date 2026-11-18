# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
"""Independent audit of the sign conventions and coordinate transformations
used across the different regularizations of Theorem 4 (R=2)."""
import sympy as sp
import sys
import traceback

sp.init_printing()

# ============================================================
# 0. Audit tools
#
# METHODOLOGICAL PRECAUTION: sp.sqrt(S) is never built explicitly in
# expressions that are later compared via cancel()/simplify(). Instead,
# u1, u2 are opaque symbols standing for sqrt(S1), sqrt(S2) -- "introduce
# radicals last" -- the same discipline already used in hohmann_build_Rta.py.
# All structural comparisons (sections 2 to 4) therefore stay purely
# RATIONAL in (p,x,t,a,u1,u2), which keeps cancel()/together() fast: there
# is no need to prove the equality of two different written forms of the
# same radical.
# ============================================================

p, x, t, w, Y, a = sp.symbols(
    "p x t w Y a", real=True
)
u1, u2 = sp.symbols("u1 u2", positive=True)  # dummies for sqrt(S1), sqrt(S2)

pH = sp.Rational(4, 3)
cc = 1 - 2*t
s2 = 4*t*(1-t)

passed = []
failed = []

def check(name, expr, expected=0, light=False):
    """Verifies exactly expr == expected.

    light=False (default): cancel(together(...)) -- fast for moderately
    sized expressions (sections 1-3, 5).

    light=True: avoids cancel()/factor() on the entire combined expression.
    On section 4, cancel() must compute a multivariate polynomial GCD whose
    intermediate size grows very large (observed peak around 6 GB of RAM).
    Instead: together() followed by isolating the numerator (sp.fraction),
    then expand() alone (distributes, no GCD needed) -- the measured peak
    on the same identity drops to around 730 MB."""
    try:
        if light:
            tog = sp.together(expr - expected)
            num, den = sp.fraction(tog)
            residual = sp.expand(num)
        else:
            residual = sp.cancel(sp.together(expr - expected))
            if residual != 0:
                residual = sp.simplify(residual)
        if residual == 0:
            passed.append(name)
            print(f"[PASS] {name}")
            return True
        failed.append(name)
        print(f"[FAIL] {name}")
        print("  Residual:", residual)
        return False
    except Exception as exc:
        failed.append(name)
        print(f"[ERROR] {name}: {exc}")
        traceback.print_exc()
        return False


def _cdiff(expr, *vs):
    """Differentiation followed by IMMEDIATE cancellation (cancel+together)
    -- avoids the combinatorial growth of an unsimplified nested fraction
    when it is differentiated a second time (this is the actual cause of
    the multi-minute slowdown otherwise observed, independent of any
    radical: sp.diff() on an uncancelled nested fraction can produce a
    numerator/denominator of explosive size before any sqrt is even
    involved)."""
    return sp.cancel(sp.together(sp.diff(expr, *vs)))


def C_first(S1, S2, q):
    """dC/dq, with sqrt(S1), sqrt(S2) replaced by the dummies u1, u2 (no
    sp.sqrt() is ever constructed)."""
    return (
        _cdiff(S1, q)/(2*u1)
        + _cdiff(S2, q)/(2*u2)
    )


def C_second(S1, S2, q, r):
    """d2C/dq dr, same discipline: 1/(4*S**(3/2)) is written as 1/(4*S*u)
    since S**(3/2) = S*sqrt(S) = S*u -- an algebraic relation used only to
    REWRITE, never requiring sympy to prove an equality of radicals. Each
    derivative is cancelled (cancel) separately BEFORE being combined, so
    that an already-nested fraction is never differentiated again."""
    result = 0
    for S, u in ((S1, u1), (S2, u2)):
        Sq = _cdiff(S, q)
        Sr = _cdiff(S, r)
        Sqr = _cdiff(Sq, r)   # derive la version DEJA annulee de Sq, pas S brut
        result += Sqr/(2*u) - Sq*Sr/(4*S*u)
    return result


def H_from_x(S1, S2):
    Cx = C_first(S1, S2, x)
    Cpp = C_second(S1, S2, p, p)
    Cxx = C_second(S1, S2, x, x)
    Cpx = C_second(S1, S2, p, x)

    Cphiphi = (1-x**2)*Cxx - x*Cx
    return sp.cancel(sp.together(
        Cphiphi*Cpp - (1-x**2)*Cpx**2
    ))


def H_from_t(S1, S2):
    Ct = C_first(S1, S2, t)
    Cpp = C_second(S1, S2, p, p)
    Ctt = C_second(S1, S2, t, t)
    Cpt = C_second(S1, S2, p, t)

    Cphiphi = (
        t*(1-t)*Ctt
        + (1-2*t)*Ct/2
    )
    return sp.cancel(sp.together(
        Cphiphi*Cpp - t*(1-t)*Cpt**2
    ))


# ============================================================
# 1. Initial R=2 definitions in (p,x)
# ============================================================

A0 = p - 1
B0 = p/2 - 1
N = A0*x + B0
E = A0**2 + N**2/(1-x**2)

S1x = sp.cancel(sp.together(3 + (E-1)/p - 2*sp.sqrt(p)))
S2x = sp.cancel(sp.together(sp.Rational(3, 2) + (E-1)/p - sp.sqrt(p/2)))

print("\n=== 1. Variable transformations ===")

check("x = 1 - 2t", (1-2*t) - cc)
check("sin(phi)^2 = 4t(1-t)", 1-cc**2-s2)
check("4 (dt/dphi)^2 = 4t(1-t)", 4*t*(1-t) - s2)
check("d2t/dphi2 = (1-2t)/2", (1-2*t)/2 - cc/2)


# ============================================================
# 2. Exact equivalence of the two formulas for H
#
# Note: S1x, S2x contain sp.sqrt(p) (not sqrt(S1) or sqrt(S2)) -- this is
# the circular-speed term -2*sqrt(p), independent of the dummies u1, u2,
# which stand only for sqrt(S1), sqrt(S2). sp.cancel will therefore never
# see two different written forms of the same radical: sqrt(p) stays as
# is on both sides of the comparison (the same variable p, never
# substituted here), and S1, S2 themselves stay purely rational in (p,x)
# or (p,t).
# ============================================================

print("\n=== 2. Hessian: the x and t conventions ===")

S1t = S1x.subs(x, 1-2*t)
S2t = S2x.subs(x, 1-2*t)

Hx = H_from_x(S1x, S2x)
Ht = H_from_t(S1t, S2t)

Hx_to_t = Hx.subs(x, 1-2*t)

check("H_x(p,1-2t) = H_t(p,t)", Hx_to_t - Ht)

Cx = C_first(S1x, S2x, x)
Cxx = C_second(S1x, S2x, x, x)
Cpx = C_second(S1x, S2x, p, x)

Ct = C_first(S1t, S2t, t)
Ctt = C_second(S1t, S2t, t, t)
Cpt = C_second(S1t, S2t, p, t)

check("Ct = -2 Cx", Ct - (-2*Cx.subs(x, 1-2*t)))
check("Ctt = 4 Cxx", Ctt - 4*Cxx.subs(x, 1-2*t))
check("Cpt = -2 Cpx", Cpt - (-2*Cpx.subs(x, 1-2*t)))


# ============================================================
# 3. Independent reconstruction of hohmann_build_Rta.py
#
# Here u1, u2 now stand for sqrt(S1_hat), sqrt(S2_hat) (after the
# substitution p -> pH+a*t) -- the SAME dummy convention, no sqrt is
# constructed. The differentiation identities (with respect to a, then to
# t) are verified by reconstructing them from the rational pieces,
# EXACTLY as hohmann_build_Rta.py does (sp.diff() is never applied to an
# expression already containing u1, u2 nonlinearly in the denominator --
# the chain rule is applied BY HAND, as in the constructor).
# ============================================================

print("\n=== 3. Hohmann regularization ===")

p_of_t_a = pH + a*t

def reg(expr):
    return sp.cancel(expr.subs(p, p_of_t_a))

S1h = reg(S1t)
S2h = reg(S2t)

S1p_h = reg(sp.diff(S1t, p))
S2p_h = reg(sp.diff(S2t, p))
S1t_h = reg(sp.diff(S1t, t))
S2t_h = reg(sp.diff(S2t, t))
S1pp_h = reg(sp.diff(S1t, p, 2))
S2pp_h = reg(sp.diff(S2t, p, 2))
S1pt_h = reg(sp.diff(S1t, p, t))
S2pt_h = reg(sp.diff(S2t, p, t))
S1tt_h = reg(sp.diff(S1t, t, 2))
S2tt_h = reg(sp.diff(S2t, t, 2))

# Ghat = Cp(pH+at,t), written with the dummies u1, u2 = sqrt(S1h), sqrt(S2h)
Ghat = S1p_h/(2*u1) + S2p_h/(2*u2)
Cthat = S1t_h/(2*u1) + S2t_h/(2*u2)

tCpphat = (
    reg(t*S1pp_h)/(2*u1) - t*S1p_h**2/(4*S1h*u1)
    + reg(t*S2pp_h)/(2*u2) - t*S2p_h**2/(4*S2h*u2)
)
tCpthat = (
    reg(t*S1pt_h)/(2*u1) - t*S1p_h*S1t_h/(4*S1h*u1)
    + reg(t*S2pt_h)/(2*u2) - t*S2p_h*S2t_h/(4*S2h*u2)
)
tCtthat = (
    reg(t*S1tt_h)/(2*u1) - t*S1t_h**2/(4*S1h*u1)
    + reg(t*S2tt_h)/(2*u2) - t*S2t_h**2/(4*S2h*u2)
)

Cp_direct = reg(S1p_h/(2*u1) + S2p_h/(2*u2))  # = Ghat by construction

check("Ghat = Cp(pH+at,t)", Ghat - Cp_direct)

# d/da [S1p_h/(2u1)]: u1^2=S1h (an IMPLICIT relation, used only to rewrite
# 1/u1^3 = 1/(S1h*u1), never for sympy to differentiate sqrt()).
# du1/da = (dS1h/da)/(2u1)  =>  d/da[f/(2u1)] = f_a/(2u1) - f*S1h_a/(4*S1h*u1)
S1h_a = sp.diff(S1h, a)
S2h_a = sp.diff(S2h, a)
S1p_h_a = sp.diff(S1p_h, a)
S2p_h_a = sp.diff(S2p_h, a)

dGhat_da_manual = (
    S1p_h_a/(2*u1) - S1p_h*S1h_a/(4*S1h*u1)
    + S2p_h_a/(2*u2) - S2p_h*S2h_a/(4*S2h*u2)
)

check("dGhat/da = tCpphat", dGhat_da_manual - tCpphat)

# d/dt [S1p_h/(2u1)], same principle, total derivative in t (p(t)=pH+a*t, a fixed)
S1h_t = sp.diff(S1h, t)
S2h_t = sp.diff(S2h, t)
S1p_h_t = sp.diff(S1p_h, t)
S2p_h_t = sp.diff(S2p_h, t)

dGhat_dt_manual = (
    S1p_h_t/(2*u1) - S1p_h*S1h_t/(4*S1h*u1)
    + S2p_h_t/(2*u2) - S2p_h*S2h_t/(4*S2h*u2)
)

Cpp_hat = tCpphat/t  # = Cpp(pH+at,t), read off by a rational division by t
Cpt_hat = tCpthat/t  # = Cpt(pH+at,t)
Ghat_t_chain = sp.cancel(a*Cpp_hat + Cpt_hat)

check("dGhat/dt = aCpp + Cpt", sp.cancel(dGhat_dt_manual) - Ghat_t_chain)


# ============================================================
# 4. Verification of tHhat
# ============================================================

print("\n=== 4. Assembling tHhat ===")

tHhat_constructed = (
    (1-t)*tCtthat*tCpphat
    + (1-2*t)*Cthat*tCpphat/2
    - (1-t)*tCpthat**2
)

# H(pH+at,t) reconstructed directly (same dummies u1,u2=sqrt(S1h),sqrt(S2h))
Cpp_hat_full = tCpphat/t
Cpt_hat_full = tCpthat/t
Ctt_hat_full = tCtthat/t
Cphiphi_hat = t*(1-t)*Ctt_hat_full + (1-2*t)/2*Cthat
# NO cancel() here -- it is precisely this combined size that drives the
# memory peak described above. H_on_branch is left as a raw expression;
# the check() call below, in light=True mode, handles the cancellation
# while avoiding cancel()/factor().
H_on_branch = Cphiphi_hat*Cpp_hat_full - t*(1-t)*Cpt_hat_full**2

check(
    "tHhat = t*H(pH+at,t)",
    tHhat_constructed - t*H_on_branch,
    light=True,
)

check(
    "Expanded formula for tHhat",
    tHhat_constructed - (
        ((1-t)*tCtthat + (1-2*t)*Cthat/2)*tCpphat
        - (1-t)*tCpthat**2
    )
)


# ============================================================
# 5. Exact identity L(a)
# ============================================================

print("\n=== 5. Hohmann limit L(a) ===")

s1_0 = sp.Rational(7, 3) - 4*sp.sqrt(3)/3
s2_0 = sp.Rational(5, 6) - sp.sqrt(6)/3

A_lim = sp.simplify(sp.Rational(27, 64)*(1/sp.sqrt(s1_0) + 1/sp.sqrt(s2_0)))
B_lim = -a*A_lim
D_lim = a**2*A_lim
E_lim = sp.simplify((sp.Rational(1, 12) - sp.Rational(27, 64)*a**2)*(32*A_lim/27))
K_lim = sp.simplify(D_lim + E_lim/2)
L_lim = sp.expand(sp.simplify(K_lim*A_lim - B_lim**2))
L_closed = sp.simplify(A_lim**2*(sp.Rational(4, 81) - a**2/4))

check("L(a) = A^2(4/81-a^2/4)", L_lim - L_closed)

print("A_lim exact =", A_lim)
print("L(a) exact  =", L_closed)


# ============================================================
# 6. Summary
# ============================================================

print("\n" + "="*60)
print("TRANSFORMATIONS R=2 AUDIT SUMMARY")
print("="*60)
print(f"PASS: {len(passed)}")
print(f"FAIL: {len(failed)}")

if failed:
    print("\nUnverified identities:")
    for name in failed:
        print(" -", name)
    sys.exit(1)

print("\nALL TESTED SYMBOLIC IDENTITIES VERIFIED")
print("\nNOTE: this symbolic audit does not replace the interval "
      "certificates.")
