# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
"""
pilot_mesh_audit_level2.py

A second level of audit. The structural audit (pilot_mesh_audit.py)
checks the consistency of the numbers ALREADY ARCHIVED, but recomputes
nothing numerically -- by itself it would trust a corrupted archive, or
one produced by different code than the one it claims to be from.

This script RECONSTRUCTS the symbolic chain (the same exact construction
as full_mesh_bridge.py -- S1, S2, P, the regularized S1p/S2p, G, G_u,
Q_w, Q_Y) and, for EVERY tile in the archive, INDEPENDENTLY RE-EVALUATES
by interval arithmetic, starting only from the archived W, Y0, Y_B
bounds (not from the archived G_Y/S1/S2/Q/T values themselves):

  - G_Y, S1, S2, the denominator of G, G at the Y_lo/Y_hi edges,
  - Q0 (at the tile center), T (total derivative), the final Q
    (mean-value interval),

and checks that THESE RE-EVALUATED VALUES (not the archived ones)
satisfy every certification criterion. A discrepancy between an
archived value and its re-evaluation is also reported (it should only
arise from serialization rounding; any larger discrepancy is a serious
defect).

For the stitchings, this script also re-evaluates G_Y and G at the
edges of each boundary's archived hull.

It fails EXPLICITLY on the first defect found.
"""
import sympy as sp
import mpmath as mp
import pickle, json, sys, time

mp.mp.dps = 50
mp.iv.dps = 50

w, Y, u1, u2 = sp.symbols('w Y u1 u2', positive=True)

S1 = (8*Y**2*w**3 - 12*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 16*Y*w + 12*Y + 1) / (4*Y*(2 - Y*w**2))
S2 = (2*sp.sqrt(2)*Y**2*w**3 - 6*Y**2*w**2 + 4*Y**2 + 4*Y*w**2 - 4*sp.sqrt(2)*Y*w + 1) / (4*Y*(2 - Y*w**2))

with open('/tmp/P_built.pkl', 'rb') as f:
    P = pickle.load(f)

p, x = sp.symbols('p x')
A = p - 1; B = p/2 - 1; N = A*x + B; E = A**2 + N**2/(1 - x**2)
S1_orig = 3 + (E - 1)/p - 2*sp.sqrt(p)
S2_orig = sp.Rational(3, 2) + (E - 1)/p - sp.sqrt(p/2)
subs = {p: w**2, x: -1 + w**2*Y}
S1p = sp.simplify(sp.diff(S1_orig, p).subs(subs, simultaneous=True))
S2p = sp.simplify(sp.diff(S2_orig, p).subs(subs, simultaneous=True))

G_u = sp.together(w**2 * (S1p/(2*u1) + S2p/(2*u2)))
G = sp.together(G_u.subs({u1: sp.sqrt(S1), u2: sp.sqrt(S2)}, simultaneous=True))
assert sp.simplify(G - sp.together(w**2 * (S1p/(2*sp.sqrt(S1)) + S2p/(2*sp.sqrt(S2))))) == 0

G_rat = sp.cancel(sp.together(G))
num_G, den_G = sp.fraction(G_rat)

S1_w = sp.diff(S1, w); S1_Y = sp.diff(S1, Y)
S2_w = sp.diff(S2, w); S2_Y = sp.diff(S2, Y)
u1_w = S1_w/(2*u1); u1_Y = S1_Y/(2*u1)
u2_w = S2_w/(2*u2); u2_Y = S2_Y/(2*u2)
P_w = sp.diff(P, w); P_Y = sp.diff(P, Y)
P_u1 = sp.diff(P, u1); P_u2 = sp.diff(P, u2)
Q_w = sp.cancel(P_w + P_u1*u1_w + P_u2*u2_w)
Q_Y = sp.cancel(P_Y + P_u1*u1_Y + P_u2*u2_Y)
G_w = sp.diff(G, w)
G_Y = sp.diff(G, Y)

modules_iv = [{'sqrt': mp.iv.sqrt, 'mpf': mp.iv.mpf}, 'mpmath']
S1_iv = sp.lambdify((w, Y), S1, modules=modules_iv)
S2_iv = sp.lambdify((w, Y), S2, modules=modules_iv)
Gw_iv = sp.lambdify((w, Y), G_w, modules=modules_iv)
GY_iv = sp.lambdify((w, Y), G_Y, modules=modules_iv)
Gu_iv = sp.lambdify((w, Y, u1, u2), G_u, modules=modules_iv)
P_iv = sp.lambdify((w, Y, u1, u2), P, modules=modules_iv)
Qw_iv = sp.lambdify((w, Y, u1, u2), Q_w, modules=modules_iv)
QY_iv = sp.lambdify((w, Y, u1, u2), Q_Y, modules=modules_iv)
denG_iv = sp.lambdify((w, Y), den_G, modules=modules_iv)


def I(a, b=None):
    if b is None:
        return mp.iv.mpf(str(a))
    return mp.iv.mpf([str(a), str(b)])


def mpf_pair(x):
    if isinstance(x, list):
        return mp.mpf(x[0]), mp.mpf(x[1])
    v = mp.mpf(x)
    return v, v


def fail(msg):
    print(f"AUDIT L2 FAIL: {msg}")
    sys.exit(1)


def main(path):
    with open(path) as f:
        archive = json.load(f)

    MARGIN_TOL = mp.mpf(archive['params']['MARGIN_TOL'])
    YB_WIDTH_MAX = mp.mpf(archive['params']['YB_WIDTH_MAX'])
    DENG_SEP_MIN = mp.mpf(archive['params']['DENG_SEP_MIN'])

    domain_lo, domain_hi = (mp.mpf(v) for v in archive['domain'])
    tiles = sorted(archive['tiles'], key=lambda r: mp.mpf(r['w_lo']))
    n = len(tiles)
    print(f"Independent re-evaluation of {n} tiles on [{domain_lo},{domain_hi}]...")

    t0 = time.time()
    violations = []
    archive_mismatches = []
    margins = []

    for idx, r in enumerate(tiles):
        w_lo = mp.mpf(r['w_lo']); w_hi = mp.mpf(r['w_hi']); w0 = mp.mpf(r['w0'])
        Y0lo, Y0hi = mpf_pair(r['Y0'])
        Ylo, Yhi = mpf_pair(r['YB'])
        half_width = (w_hi - w_lo) / 2

        Wbox = I(w_lo, w_hi)
        YB = mp.iv.mpf([Ylo, Yhi])
        Y0box = mp.iv.mpf([Y0lo, Y0hi])

        # --- independent re-evaluation, from the archived W/Y0/YB bounds only ---
        GY_box = GY_iv(Wbox, YB)
        S1B = S1_iv(Wbox, YB); S2B = S2_iv(Wbox, YB)
        denG_box = denG_iv(Wbox, YB)

        if S1B.a <= 0 or S2B.a <= 0:
            violations.append((idx, "S1S2_not_positive_REEVAL"))
            continue
        U1B = mp.iv.sqrt(S1B); U2B = mp.iv.sqrt(S2B)

        Ylo_pt = I(Ylo); Yhi_pt = I(Yhi)
        S1lo = S1_iv(Wbox, Ylo_pt); S2lo = S2_iv(Wbox, Ylo_pt)
        S1hi = S1_iv(Wbox, Yhi_pt); S2hi = S2_iv(Wbox, Yhi_pt)
        if S1lo.a <= 0 or S2lo.a <= 0 or S1hi.a <= 0 or S2hi.a <= 0:
            violations.append((idx, "S1S2_not_positive_at_edges_REEVAL"))
            continue
        Glo = Gu_iv(Wbox, Ylo_pt, mp.iv.sqrt(S1lo), mp.iv.sqrt(S2lo))
        Ghi = Gu_iv(Wbox, Yhi_pt, mp.iv.sqrt(S1hi), mp.iv.sqrt(S2hi))

        S1_0 = S1_iv(I(w0), Y0box); S2_0 = S2_iv(I(w0), Y0box)
        if S1_0.a <= 0 or S2_0.a <= 0:
            violations.append((idx, "S1S2_not_positive_at_center_REEVAL"))
            continue
        U1_0 = mp.iv.sqrt(S1_0); U2_0 = mp.iv.sqrt(S2_0)
        Q0 = P_iv(I(w0), Y0box, U1_0, U2_0)

        Gw_box = Gw_iv(Wbox, YB)
        Yp_box = -Gw_box / GY_box
        Qw_box = Qw_iv(Wbox, YB, U1B, U2B)
        QY_box = QY_iv(Wbox, YB, U1B, U2B)
        T_box = Qw_box + QY_box * Yp_box

        Qfinal_REEVAL = Q0 + T_box * mp.iv.mpf([-half_width, half_width])
        margin_REEVAL = -Qfinal_REEVAL.b

        # --- criteria, applied to the RE-EVALUATED values (not the archived ones) ---
        if not (GY_box.b < 0):
            violations.append((idx, "GY_not_negative_REEVAL"))
        if not (Glo.a > 0):
            violations.append((idx, "G_lo_not_positive_REEVAL"))
        if not (Ghi.b < 0):
            violations.append((idx, "G_hi_not_negative_REEVAL"))
        den_sep_ok = (denG_box.a > 0 or denG_box.b < 0) and \
            min(abs(denG_box.a), abs(denG_box.b)) >= DENG_SEP_MIN
        if not den_sep_ok:
            violations.append((idx, "denG_not_separated_REEVAL"))
        if (Yhi - Ylo) > YB_WIDTH_MAX:
            violations.append((idx, "YB_too_wide_REEVAL"))
        if not (margin_REEVAL > MARGIN_TOL):
            violations.append((idx, f"margin_insufficient_REEVAL ({margin_REEVAL})"))
        margins.append(margin_REEVAL)

        # --- comparison with the archived values (should agree up to
        #     serialization rounding; a large discrepancy signals a
        #     corrupted archive or one produced by different code) ---
        Q_arch_lo, Q_arch_hi = mpf_pair(r['Q'])
        tol_compare = mp.mpf('1e-30')
        if abs(Qfinal_REEVAL.a - Q_arch_lo) > tol_compare or abs(Qfinal_REEVAL.b - Q_arch_hi) > tol_compare:
            archive_mismatches.append((idx, "Q_mismatch",
                                        (mp.nstr(Qfinal_REEVAL.a, 15), mp.nstr(Qfinal_REEVAL.b, 15)),
                                        (mp.nstr(Q_arch_lo, 15), mp.nstr(Q_arch_hi, 15))))

    if violations:
        for idx, reason in violations:
            print(f"  TILE {idx}: {reason}")
        fail(f"{len(violations)} violation(s) detected by independent RE-EVALUATION")

    if archive_mismatches:
        print(f"\n{len(archive_mismatches)} archive/re-evaluation discrepancy(ies) detected "
              f"(beyond the comparison tolerance {tol_compare}):")
        for idx, reason, reeval, arch in archive_mismatches[:10]:
            print(f"  TILE {idx}: reeval={reeval}  archive={arch}")
        fail(f"{len(archive_mismatches)} tile(s) where the archive does not match "
             f"the independent re-evaluation")

    print(f"Independent numerical re-evaluation: PASS for all {n} tiles "
          f"(no violation, no archive/re-evaluation discrepancy)")
    print(f"Minimum margin (re-evaluated): {min(margins)}")
    print(f"Maximum margin (re-evaluated): {max(margins)}")

    # --- stitchings: re-evaluation on the archived hull ---
    stitching = archive['stitching']
    stitch_violations = []
    for i, st in enumerate(stitching):
        if 'w_s' not in st:
            stitch_violations.append((i, "no_w_s_recorded"))
            continue
        w_s = mp.mpf(st['w_s'])
        Ylo_hull, Yhi_hull = mpf_pair(st['Y_hull'])
        Wpt = I(w_s)
        YBhull = mp.iv.mpf([Ylo_hull, Yhi_hull])
        GY_hull = GY_iv(Wpt, YBhull)
        if not (GY_hull.b < 0):
            stitch_violations.append((i, f"GY_not_negative_on_hull_REEVAL: {GY_hull}"))
            continue
        Ylo_pt = I(Ylo_hull); Yhi_pt = I(Yhi_hull)
        S1lo = S1_iv(Wpt, Ylo_pt); S2lo = S2_iv(Wpt, Ylo_pt)
        S1hi = S1_iv(Wpt, Yhi_pt); S2hi = S2_iv(Wpt, Yhi_pt)
        if S1lo.a <= 0 or S2lo.a <= 0 or S1hi.a <= 0 or S2hi.a <= 0:
            stitch_violations.append((i, "S1S2_not_positive_on_hull_edges_REEVAL"))
            continue
        Glo_hull = Gu_iv(Wpt, Ylo_pt, mp.iv.sqrt(S1lo), mp.iv.sqrt(S2lo))
        Ghi_hull = Gu_iv(Wpt, Yhi_pt, mp.iv.sqrt(S1hi), mp.iv.sqrt(S2hi))
        if not (Glo_hull.a > 0 and Ghi_hull.b < 0):
            stitch_violations.append((i, f"sign_bracketing_failed_on_hull_REEVAL: "
                                          f"Glo={Glo_hull} Ghi={Ghi_hull}"))

    if stitch_violations:
        for i, reason in stitch_violations:
            print(f"  STITCHING {i}: {reason}")
        fail(f"{len(stitch_violations)} invalid stitching(s) by independent re-evaluation")

    print(f"Tile-to-tile stitchings (independently re-evaluated on the hull): "
          f"PASS for all {len(stitching)} boundaries")

    print(f"\n=== LEVEL 2 AUDIT (independent numerical re-evaluation): PASS ===")
    print(f"Time: {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else '/tmp/pilot_mesh_archive.json'
    sys.exit(main(path))
