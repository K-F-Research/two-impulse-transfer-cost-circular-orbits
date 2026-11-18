# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
"""
full_mesh_bridge.py

Adaptive mesh certification over [0.00015, 0.15], using the same
construction and tolerances as the preliminary pilot run on the
sub-interval [0.00015, 0.005]:

  - initial relative half-width 0.75% (within the 0.5-1% target range);
  - three explicit subdivision criteria, with safety thresholds that are
    CONFIGURABLE and never relaxed silently:
      * G_Y < 0 strictly on the whole box;
      * signs of G strictly opposite at the two edges of Y_B;
      * sup(Q) < -MARGIN_TOL (a sign margin beyond an absolute tolerance);
    plus two additional safeguards:
      * width(Y_B) < YB_WIDTH_MAX;
      * |denominator of G| separated from 0 by at least DENG_SEP_MIN;
  - if a criterion fails, the interval is SUBDIVIDED (tolerances are
    never relaxed silently);
  - rigorous tile-to-tile stitching: at each shared boundary w_s, the
    hull of the two tiles' Y-boxes is re-checked for G_Y<0 and for the
    sign change of G at w_s -- geometric overlap alone is not sufficient;
  - a complete JSON archive of every certified quantity, suitable for
    independent audit by a second, separate program that recomputes
    nothing and only re-reads the archive (pilot_mesh_audit.py and
    pilot_mesh_audit_level2.py).
"""
import sympy as sp
import mpmath as mp
import pickle, json, time, sys

mp.mp.dps = 50
mp.iv.dps = 50

# ---------------------------------------------------------------------
# Safety thresholds (CONFIGURABLE; never relaxed silently)
# ---------------------------------------------------------------------
MARGIN_TOL = mp.mpf('1e-12')     # minimal required sign margin: -sup(Q) > MARGIN_TOL
YB_WIDTH_MAX = mp.mpf('0.05')    # maximal acceptable width of Y_B (safeguard)
DENG_SEP_MIN = mp.mpf('0.5')     # |den_G| must stay > DENG_SEP_MIN on the box
REL_HALF_WIDTH_INIT = mp.mpf('0.0075')  # 0.75% of the current point
MAX_SUBDIV_DEPTH = 40
SCRIPT_VERSION = "full_mesh_bridge.py v1 (2026-11-09), same method and tolerances as the pilot run (sha256 of pilot script: 29f1a122d5770a75322dddb5137e743d6d688af334e0f389db2b6faf7afee6b4)"
IV_DPS = 50

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
G_mp = sp.lambdify((w, Y), G, modules='mpmath')
Gw_mp = sp.lambdify((w, Y), G_w, modules='mpmath')
GY_mp = sp.lambdify((w, Y), G_Y, modules='mpmath')


def I(a, b=None):
    if b is None:
        return mp.iv.mpf(str(a))
    return mp.iv.mpf([str(a), str(b)])


def s(x):
    """Serialize an mpf/iv value as an exact string (for the JSON archive).

    An interval endpoint (x.a, x.b) under mp.iv is itself a zero-width
    ivmpf (not a plain mpf) -- nstr() applied to it directly prints
    '[v, v]' instead of 'v'. The internal _mpi_ attribute (a (lo, hi)
    tuple in raw mpf representation) is used instead, to extract the
    true flat bounds, which works both for a wide interval and for a
    zero-width endpoint."""
    if hasattr(x, '_mpi_'):
        lo = mp.mp.mpf(x._mpi_[0])
        hi = mp.mp.mpf(x._mpi_[1])
        if lo == hi:
            return mp.nstr(lo, 40)
        return [mp.nstr(lo, 40), mp.nstr(hi, 40)]
    return mp.nstr(x, 40)


def certify_branch_box(Wbox, Ylo, Yhi):
    Ybox = mp.iv.mpf([Ylo, Yhi])
    GY_box = GY_iv(Wbox, Ybox)
    if not (GY_box.b < 0):
        return False, "GY_not_negative", None
    Ylo_pt = I(Ylo)
    S1lo = S1_iv(Wbox, Ylo_pt); S2lo = S2_iv(Wbox, Ylo_pt)
    if S1lo.a <= 0 or S2lo.a <= 0:
        return False, "S1S2_not_positive_at_Ylo", None
    Glo = Gu_iv(Wbox, Ylo_pt, mp.iv.sqrt(S1lo), mp.iv.sqrt(S2lo))
    if not (Glo.a > 0):
        return False, "G_not_positive_at_Ylo", None
    Yhi_pt = I(Yhi)
    S1hi = S1_iv(Wbox, Yhi_pt); S2hi = S2_iv(Wbox, Yhi_pt)
    if S1hi.a <= 0 or S2hi.a <= 0:
        return False, "S1S2_not_positive_at_Yhi", None
    Ghi = Gu_iv(Wbox, Yhi_pt, mp.iv.sqrt(S1hi), mp.iv.sqrt(S2hi))
    if not (Ghi.b < 0):
        return False, "G_not_negative_at_Yhi", None
    return True, None, (GY_box, Glo, Ghi)


def root_box(w0_val, hint, target_digits=25):
    seed = mp.findroot(lambda z: G_mp(w0_val, z), (hint - mp.mpf('1e-4'), hint + mp.mpf('1e-4')))
    half = mp.mpf(10) ** (-target_digits)
    ok, reason, info = certify_branch_box(I(w0_val), seed - half, seed + half)
    doublings = 0
    while not ok and doublings < 300:
        half *= 2
        ok, reason, info = certify_branch_box(I(w0_val), seed - half, seed + half)
        doublings += 1
    if not ok:
        return None
    lo_h, hi_h = (half / 2 if doublings > 0 else half * mp.mpf('1e-3')), half
    for _ in range(40):
        mid_h = (lo_h + hi_h) / 2
        ok_mid, _, _ = certify_branch_box(I(w0_val), seed - mid_h, seed + mid_h)
        if ok_mid:
            hi_h = mid_h
        else:
            lo_h = mid_h
    return seed - hi_h, seed + hi_h, seed


def branch_box(Wbox, Y0lo, Y0hi, seed, w_lo, w_hi, w0_val):
    yp = -Gw_mp(w0_val, seed) / GY_mp(w0_val, seed)
    halfw = (w_hi - w_lo) / 2
    radius = max((Y0hi - Y0lo)/2, abs(yp)*halfw * mp.mpf('1.05'))
    ok, reason, info = certify_branch_box(Wbox, seed - radius, seed + radius)
    tries = 0
    while not ok and tries < 60:
        radius *= 1.3
        ok, reason, info = certify_branch_box(Wbox, seed - radius, seed + radius)
        tries += 1
    if not ok:
        return None
    lo_r, hi_r = (radius / mp.mpf('1.3') if tries > 0 else radius * mp.mpf('0.5')), radius
    for _ in range(30):
        mid_r = (lo_r + hi_r) / 2
        ok_mid, _, _ = certify_branch_box(Wbox, seed - mid_r, seed + mid_r)
        if ok_mid:
            hi_r = mid_r
        else:
            lo_r = mid_r
    return seed - hi_r, seed + hi_r


def certify_tile(a, b, hint):
    """Attempts to certify the tile [a,b]. Returns (record_dict, reason_if_failed)."""
    w0 = (a + b) / 2
    half_width = (b - a) / 2
    Wbox = I(a, b)

    rb = root_box(w0, hint)
    if rb is None:
        return None, "root_box_failed"
    Y0lo, Y0hi, seed = rb

    bb = branch_box(Wbox, Y0lo, Y0hi, seed, a, b, w0)
    if bb is None:
        return None, "branch_box_failed"
    Ylo, Yhi = bb
    YB = mp.iv.mpf([Ylo, Yhi])

    GY_box = GY_iv(Wbox, YB)
    S1B = S1_iv(Wbox, YB); S2B = S2_iv(Wbox, YB)
    denG_box = denG_iv(Wbox, YB)

    if not (GY_box.b < 0):
        return None, "GY_not_negative_on_tile"
    if not (S1B.a > 0 and S2B.a > 0):
        return None, "S1S2_not_positive_on_tile"
    denG_min_abs = min(abs(denG_box.a), abs(denG_box.b)) if (denG_box.a > 0 or denG_box.b < 0) else mp.mpf('0')
    if not (denG_box.a > 0 or denG_box.b < 0) or denG_min_abs < DENG_SEP_MIN:
        return None, "denG_not_separated_enough"
    width_YB = Yhi - Ylo
    if width_YB > YB_WIDTH_MAX:
        return None, "YB_too_wide"

    # re-extract the sign change of G at the edges explicitly for the
    # archive (already checked inside certify_branch_box above)
    Ylo_pt = I(Ylo); Yhi_pt = I(Yhi)
    S1lo = S1_iv(Wbox, Ylo_pt); S2lo = S2_iv(Wbox, Ylo_pt)
    Glo = Gu_iv(Wbox, Ylo_pt, mp.iv.sqrt(S1lo), mp.iv.sqrt(S2lo))
    S1hi = S1_iv(Wbox, Yhi_pt); S2hi = S2_iv(Wbox, Yhi_pt)
    Ghi = Gu_iv(Wbox, Yhi_pt, mp.iv.sqrt(S1hi), mp.iv.sqrt(S2hi))

    U1B = mp.iv.sqrt(S1B); U2B = mp.iv.sqrt(S2B)
    Y0box = mp.iv.mpf([Y0lo, Y0hi])
    S1_0 = S1_iv(I(w0), Y0box); S2_0 = S2_iv(I(w0), Y0box)
    U1_0 = mp.iv.sqrt(S1_0); U2_0 = mp.iv.sqrt(S2_0)
    Q0 = P_iv(I(w0), Y0box, U1_0, U2_0)

    Gw_box = Gw_iv(Wbox, YB)
    Yp_box = -Gw_box / GY_box
    Qw_box = Qw_iv(Wbox, YB, U1B, U2B)
    QY_box = QY_iv(Wbox, YB, U1B, U2B)
    T_box = Qw_box + QY_box * Yp_box

    Qfinal = Q0 + T_box * mp.iv.mpf([-half_width, half_width])
    margin = -Qfinal.b

    if not (margin > MARGIN_TOL):
        return None, "margin_insufficient"

    record = {
        "w_lo": s(a), "w_hi": s(b), "w0": s(w0),
        "Y0": [s(Y0lo), s(Y0hi)],
        "YB": [s(Ylo), s(Yhi)],
        "GY": [s(GY_box.a), s(GY_box.b)],
        "G_lo": [s(Glo.a), s(Glo.b)],
        "G_hi": [s(Ghi.a), s(Ghi.b)],
        "S1": [s(S1B.a), s(S1B.b)],
        "S2": [s(S2B.a), s(S2B.b)],
        "denG": [s(denG_box.a), s(denG_box.b)],
        "T": [s(T_box.a), s(T_box.b)],
        "Q0": [s(Q0.a), s(Q0.b)],
        "Q": [s(Qfinal.a), s(Qfinal.b)],
        "margin": s(margin),
        "seed_Y": s(seed),
    }
    return record, None


def process_interval(a, b, hint, depth, accepted, failures):
    if depth > MAX_SUBDIV_DEPTH:
        failures.append({"w_lo": s(a), "w_hi": s(b), "reason": "max_depth_exceeded"})
        raise RuntimeError(f"max depth exceeded at [{a},{b}]")
    record, reason = certify_tile(a, b, hint)
    if record is not None:
        accepted.append(record)
        return mp.mpf(record["seed_Y"])
    # subdivide
    mid = (a + b) / 2
    seed_mid = process_interval(a, mid, hint, depth + 1, accepted, failures)
    seed_right = process_interval(mid, b, seed_mid, depth + 1, accepted, failures)
    return seed_right


if __name__ == "__main__":
    t0 = time.time()
    left = mp.mpf('0.00015')
    right = mp.mpf('0.15')    # full domain of this region (method and tolerances unchanged from the pilot run)

    b_slope = mp.mpf('-1') / 3 - mp.sqrt(2) / 6
    hint0 = mp.mpf('0.5') + b_slope * left

    accepted = []
    failures = []

    # build the initial mesh at constant relative half-width
    # REL_HALF_WIDTH_INIT, subdividing automatically whenever a criterion fails
    w_cur = left
    hint = hint0
    tile_boundaries = [left]
    while w_cur < right:
        step = REL_HALF_WIDTH_INIT * 2 * w_cur
        w_next = min(w_cur + step, right)
        tile_boundaries.append(w_next)
        w_cur = w_next

    for i in range(len(tile_boundaries) - 1):
        a, b = tile_boundaries[i], tile_boundaries[i + 1]
        try:
            hint = process_interval(a, b, hint, 0, accepted, failures)
        except RuntimeError as e:
            print(f"Unrecoverable failure on [{a},{b}]: {e}")
            failures.append({"w_lo": s(a), "w_hi": s(b), "reason": str(e)})

    accepted.sort(key=lambda r: mp.mpf(r["w_lo"]))

    print(f"Accepted tiles: {len(accepted)}")
    print(f"Unrecoverable failures: {len(failures)}")
    def _parse_margin(m):
        return mp.mpf(m[0]) if isinstance(m, list) else mp.mpf(m)
    print(f"Minimum margin: {min(_parse_margin(r['margin']) for r in accepted) if accepted else 'N/A'}")
    print(f"Time: {time.time()-t0:.2f}s")

    # ---------------------------------------------------------------
    # Tile-to-tile stitching: hull of the Y_B boxes at each boundary
    # ---------------------------------------------------------------
    stitching_results = []
    gap_or_overlap = False
    for i in range(1, len(accepted)):
        prev_r, cur_r = accepted[i-1], accepted[i]
        if mp.mpf(prev_r["w_hi"]) != mp.mpf(cur_r["w_lo"]):
            gap_or_overlap = True
            stitching_results.append({"between": [i-1, i], "status": "GAP_OR_OVERLAP_IN_W"})
            continue
        w_s = mp.mpf(prev_r["w_hi"])
        Y_prev = (mp.mpf(prev_r["YB"][0]), mp.mpf(prev_r["YB"][1]))
        Y_cur = (mp.mpf(cur_r["YB"][0]), mp.mpf(cur_r["YB"][1]))
        Ylo_hull = min(Y_prev[0], Y_cur[0])
        Yhi_hull = max(Y_prev[1], Y_cur[1])
        ok, reason, info = certify_branch_box(I(w_s), Ylo_hull, Yhi_hull)
        stitching_results.append({
            "between": [i-1, i], "w_s": s(w_s),
            "Y_hull": [s(Ylo_hull), s(Yhi_hull)],
            "status": "OK" if ok else f"FAIL:{reason}",
        })
        if not ok:
            gap_or_overlap = True

    n_stitch_fail = sum(1 for sres in stitching_results if sres["status"] != "OK")
    print(f"Tile-to-tile stitchings tested: {len(stitching_results)}, failures: {n_stitch_fail}")

    # ---------------------------------------------------------------
    # Complete archive
    # ---------------------------------------------------------------
    archive = {
        "script_version": SCRIPT_VERSION,
        "mp_dps": mp.mp.dps,
        "iv_dps": IV_DPS,
        "params": {
            "MARGIN_TOL": s(MARGIN_TOL),
            "YB_WIDTH_MAX": s(YB_WIDTH_MAX),
            "DENG_SEP_MIN": s(DENG_SEP_MIN),
            "REL_HALF_WIDTH_INIT": s(REL_HALF_WIDTH_INIT),
            "MAX_SUBDIV_DEPTH": MAX_SUBDIV_DEPTH,
        },
        "domain": [s(left), s(right)],
        "n_tiles": len(accepted),
        "n_failures": len(failures),
        "tiles": accepted,
        "failures": failures,
        "stitching": stitching_results,
        "elapsed_seconds": time.time() - t0,
    }
    with open('/tmp/full_mesh_archive.json', 'w') as f:
        json.dump(archive, f, indent=1)
    print("Archive written: /tmp/full_mesh_archive.json")
