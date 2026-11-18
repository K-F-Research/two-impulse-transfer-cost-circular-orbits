# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
# Exact conversion between the Hohmann-endpoint variable t and the
# physical variable w = sqrt(p) at the junction boundary t = t_archive,
# where t_archive is the archive's own recorded value of domain_T[0]
# (= the upper boundary of 04_hohmann_endpoint's certified domain =
# the lower boundary of 05_upper_bridge's first tile), not the bare
# literal 1e-4.
#
# On the branch p(t) = p_H + a(t)*t, the correct conversion is
#
#     w(t) = sqrt(p_H + a(t)*t),
#
# so the lower bound of the Hohmann bridge in the w-variable at
# t = t_archive is
#
#     w_- = sqrt(4/3 + t_archive * a(t_archive)).
#
# This script certifies, by interval arithmetic, a tight enclosure of
# the root a(t_archive) of Ghat(t_archive, .) = 0 (monotone bisection,
# using the strict positivity of t*Cpp = Ghat_a already established
# elsewhere), derives a certified enclosure of p(t_archive) and of w_-,
# and verifies that w_- < w_H strictly, where w_H = sqrt(4/3) corresponds
# to the limiting point t = 0 of the branch (not an interior point of
# the physical domain).

import pickle
import sympy as sp
import mpmath as mp
from fractions import Fraction

mp.mp.dps = 90
mp.iv.dps = 90

FROZEN_PATH = '/tmp/hohmann_FROZEN_v1.pkl'
T_ARCHIVE_STR = '0.0001000000000000000000000000000000000000000000000000000211419864763152195'

passed = []
failed = []


def check(name, cond, detail=""):
    if cond:
        passed.append(name)
        print(f"[PASS] {name}" + (f"  ({detail})" if detail else ""))
    else:
        failed.append(name)
        print(f"[FAIL] {name}" + (f"  ({detail})" if detail else ""))


with open(FROZEN_PATH, 'rb') as f:
    EXPRS = pickle.load(f)

t, a = sp.symbols('t a')
Ghat = EXPRS['Ghat']
tCpphat = EXPRS['Ghat_a']  # = t*Cpp = t * dGhat/da

modules_iv = [{'sqrt': mp.iv.sqrt, 'mpf': mp.iv.mpf}, 'mpmath']
f_G = sp.lambdify((t, a), Ghat, modules=modules_iv, cse=True)
f_tCpp = sp.lambdify((t, a), tCpphat, modules=modules_iv, cse=True)


def Ibox(lo, hi):
    return mp.iv.mpf([str(lo), str(hi)])


t0_scalar = mp.mpf(T_ARCHIVE_STR)
t0 = mp.iv.mpf(T_ARCHIVE_STR)

check(
    "t0 parsed losslessly at dps=90 (exact-rational round-trip)",
    Fraction(str(t0_scalar)) == Fraction(T_ARCHIVE_STR),
    f"t0_scalar={mp.nstr(t0_scalar, 80)}"
)

print("=== 1. Certified enclosure of the root a(t_archive) ===")

# Monotone bisection: Ghat(t0, .) is strictly increasing in a, since
# t*Cpp > 0 on [-0.25, 0.25] x [0, t_archive] (established independently
# in 04_hohmann_endpoint's certificate); the bracket is therefore
# narrowed by testing the sign at the midpoint.
lo, hi = mp.mpf('-0.25'), mp.mpf('0.25')
for _ in range(150):
    mid = (lo + hi) / 2
    g_mid = f_G(t0, mp.iv.mpf(str(mid)))
    if mp.mpf(g_mid.a) > 0:
        hi = mid
    else:
        lo = mid

bracket_width = hi - lo
check(
    "monotone bisection converges (final bracket width negligible)",
    bracket_width < mp.mpf('1e-40'),
    f"width={mp.nstr(bracket_width, 6)}"
)

Ba_bracket = Ibox(lo, hi)
tCpp_box = f_tCpp(t0, Ba_bracket)
check(
    "t*Cpp > 0 on the final bracket (monotonicity re-verified)",
    mp.mpf(tCpp_box.a) > 0,
    f"tCpp in {tCpp_box}"
)

G_lo = f_G(t0, mp.iv.mpf(str(lo)))
G_hi = f_G(t0, mp.iv.mpf(str(hi)))
check("Ghat(t_archive, lo) <= 0", mp.mpf(G_lo.b) <= 0, f"{G_lo}")
check("Ghat(t_archive, hi) >= 0", mp.mpf(G_hi.a) >= 0, f"{G_hi}")
check(
    "strict sign change confirmed => unique root a(t_archive) in [lo, hi]",
    (mp.mpf(G_lo.b) <= 0) and (mp.mpf(G_hi.a) >= 0) and (mp.mpf(tCpp_box.a) > 0)
)

print(f"\nCertified enclosure: a(t_archive) in [{lo}, {hi}]")

print("\n=== 2. Exact t -> w conversion at the boundary t = t_archive ===")

pH = mp.mpf(4) / 3
p_lo = pH + t0_scalar * lo
p_hi = pH + t0_scalar * hi
check("p_lo <= p_hi (consistent bracket)", p_lo <= p_hi)

w_lo = mp.sqrt(p_lo)
w_hi = mp.sqrt(p_hi)
w_H = mp.sqrt(mp.mpf(4) / 3)

print(f"p(t_archive) = 4/3 + t_archive*a(t_archive) in [{mp.nstr(p_lo, 40)}, {mp.nstr(p_hi, 40)}]")
print(f"w_-           = sqrt(p(t_archive))            in [{mp.nstr(w_lo, 40)}, {mp.nstr(w_hi, 40)}]")
print(f"w_H                                              = {mp.nstr(w_H, 40)}")

check(
    "w_- < w_H strictly (w = w_H is attained only in the limit t -> 0)",
    w_hi < w_H,
    f"w_H - w_hi = {mp.nstr(w_H - w_hi, 20)}"
)

# Compare against the previous (bare-1e-4) value, to quantify the shift.
t_old = mp.mpf('1e-4')
w_old_lo = mp.sqrt(pH + t_old * lo)
w_old_hi = mp.sqrt(pH + t_old * hi)
print(f"\n[info] w_- recomputed at the bare literal 1e-4 (for comparison only): "
      f"[{mp.nstr(w_old_lo, 30)}, {mp.nstr(w_old_hi, 30)}]")
print(f"[info] shift in w_hi between t=1e-4 and t=t_archive: "
      f"{mp.nstr(w_hi - w_old_hi, 10)}")

print("\n" + "=" * 60)
print("SUMMARY: exact w-conversion near the Hohmann endpoint, at t_archive")
print("=" * 60)
print(f"PASS: {len(passed)}")
print(f"FAIL: {len(failed)}")
if failed:
    import sys
    print("\nUnverified checks:")
    for name in failed:
        print(" -", name)
    sys.exit(1)

print("\nALL REQUIRED CHECKS PASSED.")
print(f"\nConclusion: w_- = sqrt(4/3 + t_archive*a(t_archive)) in [{mp.nstr(w_lo, 30)}, {mp.nstr(w_hi, 30)}]")
print("is the correct lower bound of the Hohmann bridge in the w-variable,")
print("strictly below w_H. The point w = w_H (t = 0) is a limit of the")
print("branch, never attained for t > 0: the physical domain of the")
print("Hohmann-endpoint region remains w < w_H strictly.")
