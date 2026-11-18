# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
# Junction certificate between the Hohmann bridge (variables t, a),
# certified by hohmann_full_mesh.py and audited by
# audit_hohmann_archive.py, and the central mesh (variables p, phi) of
# mesh_certificate_R2_v6.py, at the shared boundary w = 1.144701.
#
# The central mesh does not persist individual phi-boxes (only
# aggregate counts), so a fresh tile ending exactly at w = 1.144701 is
# re-certified here, independently, using the same certify_tile()
# construction as the central mesh.
#
# A phi-only overlap between the two sides is not sufficient to
# identify the two stationary roots as the same physical point: two
# phi-boxes can overlap while their corresponding p-intervals are
# disjoint. The junction is therefore established in the full (p, phi)
# plane:
#   1. the p-boxes on both sides have a non-empty intersection
#      (P_overlap);
#   2. on P_overlap x hull(phi_Hohmann, phi_central), the partial
#      derivative C_pphi has strictly constant sign, giving uniqueness
#      of phi at fixed p on the overlap;
#   3. each side independently establishes existence of its own
#      stationary root on its own box, valid a fortiori on the
#      narrower overlap sub-box.
# Uniqueness at fixed p, combined with existence on both sides, forces
# the two roots to be the same physical (p, phi) point.

import pickle, json, sys, hashlib
import sympy as sp
import mpmath as mp

mp.iv.dps = 50
mp.mp.dps = 60

ARCHIVE_PATH = '/tmp/hohmann_full_mesh_archive.json'
FROZEN_PATH = '/tmp/hohmann_FROZEN_v1.pkl'
SECURITY_EPS_EXPONENT = 40

W_JOIN = mp.mpf('1.144701')  # exact boundary of the central mesh (UP side)
P_H = sp.Rational(4, 3)      # p_H = 4/3 (R = 2 Hohmann point)

# ---------------------------------------------------------------------
# 1. Symbolic (p, phi) model and certify_tile(), as in
#    mesh_certificate_R2_v6.py (same definitions, same precisions).
# ---------------------------------------------------------------------
p, phi = sp.symbols('p phi', positive=True)
R = sp.Integer(2)
r1 = sp.Integer(1)
r2 = R
c = sp.cos(phi)
s = sp.sin(phi)
x = p / r1 - 1
y = p / r2 - 1
e2 = x**2 + (x * c + y)**2 / s**2
vtheta1 = sp.sqrt(p) / r1
vtheta2 = sp.sqrt(p) / r2
vc1 = 1 / sp.sqrt(r1)
vc2 = 1 / sp.sqrt(r2)
v1_2 = 2 / r1 + (e2 - 1) / p
v2_2 = 2 / r2 + (e2 - 1) / p
S1 = v1_2 + vc1**2 - 2 * vc1 * vtheta1
S2 = v2_2 + vc2**2 - 2 * vc2 * vtheta2

C = sp.sqrt(S1) + sp.sqrt(S2)
Cp = sp.diff(C, p)
Cphi = sp.diff(C, phi)
Cpp = sp.diff(Cp, p)
Cpphi = sp.diff(Cp, phi)
Cphiphi = sp.diff(Cphi, phi)
H = Cphiphi * Cpp - Cpphi**2
Hp = sp.diff(H, p)
Hphi = sp.diff(H, phi)

f_Cp = sp.lambdify((p, phi), Cp, 'mpmath')
f_H = sp.lambdify((p, phi), H, 'mpmath')
f_Cpphi = sp.lambdify((p, phi), Cpphi, 'mpmath')
f_Cpp = sp.lambdify((p, phi), Cpp, 'mpmath')
f_Hp = sp.lambdify((p, phi), Hp, 'mpmath')
f_Hphi = sp.lambdify((p, phi), Hphi, 'mpmath')

iv_Cpp = sp.lambdify((p, phi), Cpp, mp.iv)
iv_Cpphi = sp.lambdify((p, phi), Cpphi, mp.iv)
iv_Hp = sp.lambdify((p, phi), Hp, mp.iv)
iv_Hphi = sp.lambdify((p, phi), Hphi, mp.iv)
iv_Cp = sp.lambdify((p, phi), Cp, mp.iv)
iv_H = sp.lambdify((p, phi), H, mp.iv)
iv_S1 = sp.lambdify((p, phi), S1, mp.iv)
iv_S2 = sp.lambdify((p, phi), S2, mp.iv)

_s1s2_check_counts = {"passed": 0, "failed": 0}


def check_S1_S2_positive(Bp, B_phi):
    try:
        S1_box = iv_S1(Bp, B_phi)
        S2_box = iv_S2(Bp, B_phi)
    except (ValueError, ZeroDivisionError):
        _s1s2_check_counts["failed"] += 1
        return False
    ok = (mp.mpf(S1_box.a) > 0) and (mp.mpf(S2_box.a) > 0)
    if ok:
        _s1s2_check_counts["passed"] += 1
    else:
        _s1s2_check_counts["failed"] += 1
    return ok


def iv_point(val):
    return mp.iv.mpf(val)


def newton_phi(p0, phi_guess):
    ph = phi_guess
    for _ in range(25):
        val = f_Cp(p0, ph)
        der = f_Cpphi(p0, ph)
        step = val / der
        ph = ph - step
        if abs(step) < mp.mpf('1e-40'):
            break
    return ph


def _newton_contract(Bp, dp, p0, phi_center_guess, phi_bracket_halfwidth):
    B_phi = mp.iv.mpf([phi_center_guess - phi_bracket_halfwidth,
                        phi_center_guess + phi_bracket_halfwidth])
    B_phi_start = B_phi
    prev_width = None
    try:
        for _ in range(40):
            phi_mid = (mp.mpf(B_phi.a) + mp.mpf(B_phi.b)) / 2
            Cp0 = iv_Cp(iv_point(p0), iv_point(phi_mid))
            Cpp_box = iv_Cpp(Bp, mp.iv.mpf(phi_mid))
            num = Cp0 + Cpp_box * dp
            den = iv_Cpphi(Bp, B_phi)
            den_a, den_b = mp.mpf(den.a), mp.mpf(den.b)
            if den_a <= 0 <= den_b:
                return None
            cand = mp.iv.mpf(phi_mid) - num / den
            lo = max(mp.mpf(cand.a), mp.mpf(B_phi.a))
            hi = min(mp.mpf(cand.b), mp.mpf(B_phi.b))
            if lo > hi:
                return None
            new_B = mp.iv.mpf([lo, hi])
            width_new = hi - lo
            B_phi = new_B
            if prev_width is not None and width_new >= prev_width:
                break
            prev_width = width_new
    except (ValueError, ZeroDivisionError):
        return None

    Bs_lo, Bs_hi = mp.mpf(B_phi_start.a), mp.mpf(B_phi_start.b)
    Bf_lo, Bf_hi = mp.mpf(B_phi.a), mp.mpf(B_phi.b)
    if not (Bf_lo > Bs_lo and Bf_hi < Bs_hi):
        return None
    return B_phi


def certify_tile(w_lo, w_hi, phi_center_guess, phi_bracket_halfwidth):
    Bp = mp.iv.mpf([w_lo, w_hi]) ** 2
    p0 = (mp.mpf(Bp.a) + mp.mpf(Bp.b)) / 2
    dp = Bp - mp.iv.mpf(p0)

    hw_cap = mp.mpf('0.2')
    B_phi = None
    hw = phi_bracket_halfwidth
    for _ in range(10):
        if hw > hw_cap:
            break
        B_phi = _newton_contract(Bp, dp, p0, phi_center_guess, hw)
        if B_phi is not None:
            break
        hw *= 3
    if B_phi is None:
        return False, None, None, None, None

    if not check_S1_S2_positive(Bp, B_phi):
        return False, B_phi, None, None, None

    Bf_lo, Bf_hi = mp.mpf(B_phi.a), mp.mpf(B_phi.b)
    phi0 = (Bf_lo + Bf_hi) / 2
    H0_iv = iv_H(iv_point(p0), iv_point(phi0))

    dp_half_iv = abs(dp)
    dphi_half_iv = abs(B_phi - mp.iv.mpf(phi0))
    Hp_box = iv_Hp(Bp, B_phi)
    Hphi_box = iv_Hphi(Bp, B_phi)
    Hp_mag_iv = abs(Hp_box)
    Hphi_mag_iv = abs(Hphi_box)
    remainder_iv = Hp_mag_iv * dp_half_iv + Hphi_mag_iv * dphi_half_iv
    remainder = mp.mpf(remainder_iv.b)
    H0_a, H0_b = mp.mpf(H0_iv.a), mp.mpf(H0_iv.b)
    if H0_a <= 0 <= H0_b:
        H0_mag_worst = mp.mpf('0')
    else:
        H0_mag_worst = min(abs(H0_a), abs(H0_b))
    if not (remainder < H0_mag_worst):
        return False, B_phi, H0_iv, remainder, None
    sign = 1 if H0_a > 0 else -1
    return True, B_phi, H0_iv, remainder, sign


# ---------------------------------------------------------------------
# 2. Reference root (same starting point as mesh_certificate_R2_v6.py),
#    used only to seed the Newton guess, never as a certified value.
# ---------------------------------------------------------------------
w1_reported = mp.mpf('0.30698081390187630871')
p0_root = w1_reported**2
phi0_root = mp.mpf('164.3710517') * mp.pi / 180


def F(v):
    pp, ph = v[0], v[1]
    return mp.matrix([f_Cp(pp, ph), f_H(pp, ph)])


def Jmat(v):
    pp, ph = v[0], v[1]
    return mp.matrix([[f_Cpp(pp, ph), f_Cpphi(pp, ph)],
                       [f_Hp(pp, ph), f_Hphi(pp, ph)]])


z = mp.matrix([p0_root, phi0_root])
for _ in range(30):
    delta = mp.lu_solve(Jmat(z), -F(z))
    z = z + delta
    if mp.sqrt(delta[0]**2 + delta[1]**2) < mp.mpf('1e-30'):
        break
p_c, phi_c = z[0], z[1]
w_c = mp.sqrt(p_c)
print(f"Reference root (seed only): p_c={p_c}  phi_c={phi_c}  w_c={w_c}")

# ---------------------------------------------------------------------
# 3. Fresh, independent certification of a tile ending exactly at
#    W_JOIN = 1.144701 (UP boundary of the central mesh).
# ---------------------------------------------------------------------
# Deliberate overlap zone: the central tile is widened to cover a
# range of w well beyond the serialization security margin (1e-40),
# so that its phi-box extends comfortably inside the last Hohmann
# tile rather than merely touching its boundary.
dw = mp.mpf('2e-4')
w_lo_join = W_JOIN - dw
w_hi_join = W_JOIN
phi_guess_join = newton_phi(((w_lo_join + w_hi_join) / 2)**2, phi_c)
halfwidth0 = mp.mpf('5e-4')

print(f"\nCertifying a fresh central tile on [w={w_lo_join}, {w_hi_join}] "
      f"(central-mesh boundary), phi guess={phi_guess_join} ...")

ok = False
hw = halfwidth0
for attempt in range(12):
    ok, B_phi_central, H0_iv, remainder, sign = certify_tile(
        w_lo_join, w_hi_join, phi_guess_join, hw)
    if ok:
        break
    hw *= 3
    if hw > mp.mpf('0.2'):
        break

if not ok:
    print("\n*** FAILED: could not certify a fresh tile at w = 1.144701 ***")
    sys.exit(1)

print(f"Central tile certified: B_phi = {B_phi_central}")
print(f"  H0 = {H0_iv}  (sign={sign})  remainder={remainder}")
print(f"  S1, S2 checks: passed={_s1s2_check_counts['passed']} "
      f"failed={_s1s2_check_counts['failed']}")

# p-box on the central side (same formula as in certify_tile):
# Bp = [w_lo, w_hi]^2, kept here for the (p, phi) cross-check.
Bp_central = mp.iv.mpf([w_lo_join, w_hi_join]) ** 2
print(f"p-box (central side, = [w_lo,w_hi]^2): {Bp_central}")

# ---------------------------------------------------------------------
# 4. Last tile of the Hohmann bridge, converted from (t, a) to (p, phi).
# ---------------------------------------------------------------------
with open(FROZEN_PATH, 'rb') as f:
    FROZEN = pickle.load(f)
_h = hashlib.sha256()
for _k in sorted(FROZEN.keys()):
    _h.update(_k.encode())
    _h.update(str(FROZEN[_k]).encode())
frozen_sha256 = _h.hexdigest()

with open(ARCHIVE_PATH, 'r') as f:
    archive = json.load(f)
assert archive["frozen_sha256"] == frozen_sha256, "frozen expression fingerprint mismatch"

tiles = archive["tiles"]
tiles_sorted = sorted(tiles, key=lambda r: mp.mpf(r["t_lo"]))
last = tiles_sorted[-1]

eps = mp.mpf(10) ** (-SECURITY_EPS_EXPONENT)


def read_lo(x): return mp.mpf(x) - eps
def read_hi(x): return mp.mpf(x) + eps


t_hi_last = read_hi(last["t_hi"])
t_lo_last = read_lo(last["t_lo"])
a_m_last = mp.mpf(last["a_m"])
r_last = mp.mpf(last["r"])

print(f"\nLast Hohmann tile: t_lo={t_lo_last}  t_hi={t_hi_last}  "
      f"a_m={a_m_last}  r={r_last}")

a_lo = a_m_last - r_last - eps
a_hi = a_m_last + r_last + eps

p_H_mp = mp.mpf(sp.Rational(4, 3))
# p = p_H + a*t, evaluated at t = t_hi_last (certified upper boundary of
# the last tile), a ranging over [a_lo, a_hi].
p_box_lo = p_H_mp + a_lo * t_hi_last
p_box_hi = p_H_mp + a_hi * t_hi_last
if p_box_lo > p_box_hi:
    p_box_lo, p_box_hi = p_box_hi, p_box_lo

# phi = 2*arcsin(sqrt(t)), evaluated over the full tile
# [t_lo_last, t_hi_last] for consistency with the whole tile, not only
# its boundary.
phi_from_t_lo = 2 * mp.asin(mp.sqrt(t_lo_last))
phi_from_t_hi = 2 * mp.asin(mp.sqrt(t_hi_last))
phi_box_lo = min(phi_from_t_lo, phi_from_t_hi) - eps
phi_box_hi = max(phi_from_t_lo, phi_from_t_hi) + eps

p_full_candidates = [p_H_mp + a_lo * t_lo_last, p_H_mp + a_lo * t_hi_last,
                      p_H_mp + a_hi * t_lo_last, p_H_mp + a_hi * t_hi_last]
p_full_lo = min(p_full_candidates)
p_full_hi = max(p_full_candidates)

print(f"\n(t, a) -> (p, phi) conversion of the last Hohmann tile:")
print(f"  p (at boundary t=t_hi_last) in [{p_box_lo}, {p_box_hi}]")
print(f"  p (over the full tile)      in [{p_full_lo}, {p_full_hi}]")
print(f"  phi (= 2*arcsin(sqrt(t)))    in [{phi_box_lo}, {phi_box_hi}]")

B_phi_hohmann = mp.iv.mpf([phi_box_lo, phi_box_hi])

# ---------------------------------------------------------------------
# 5. Overlap check in phi.
# ---------------------------------------------------------------------
c_lo, c_hi = mp.mpf(B_phi_central.a), mp.mpf(B_phi_central.b)
h_lo, h_hi = mp.mpf(B_phi_hohmann.a), mp.mpf(B_phi_hohmann.b)
overlap_lo = max(c_lo, h_lo)
overlap_hi = min(c_hi, h_hi)
overlap_ok = overlap_lo <= overlap_hi

print(f"\nphi-box (central mesh, freshly certified): [{c_lo}, {c_hi}]")
print(f"phi-box (Hohmann bridge, last tile, converted): [{h_lo}, {h_hi}]")
print(f"Overlap: [{overlap_lo}, {overlap_hi}]" if overlap_ok else "Overlap: EMPTY")

if not overlap_ok:
    print("\n*** JUNCTION FAILED: the two phi-boxes do not overlap ***")
    sys.exit(1)

overlap_width = overlap_hi - overlap_lo
COMFORT_MARGIN = mp.mpf('1e-6')  # well above eps = 1e-40
comfortable = overlap_width > COMFORT_MARGIN
print(f"Overlap width: {overlap_width}  (margin > {COMFORT_MARGIN}: {comfortable})")
if not comfortable:
    print("\n*** WARNING: non-empty but insufficiently wide overlap ***")
    sys.exit(1)

print("\nNon-empty, comfortable phi-overlap confirmed. This alone is not")
print("sufficient to identify the two roots: two phi-boxes can overlap while")
print("corresponding to disjoint p-intervals. The additional requirement is a")
print("non-empty p-overlap together with uniqueness of phi at fixed p on it.")

# ---------------------------------------------------------------------
# 6. Rigorous check in the (p, phi) plane: non-empty p-overlap and
#    uniqueness of phi at fixed p on that overlap.
# ---------------------------------------------------------------------
print("\n=== 6. Junction check in the (p, phi) plane ===")

p_lo_h, p_hi_h = mp.mpf(min(p_box_lo, p_box_hi)), mp.mpf(max(p_box_lo, p_box_hi))
p_lo_c, p_hi_c = mp.mpf(Bp_central.a), mp.mpf(Bp_central.b)

P_overlap_lo = max(p_lo_h, p_lo_c)
P_overlap_hi = min(p_hi_h, p_hi_c)
p_overlap_ok = P_overlap_lo <= P_overlap_hi

print(f"p-box (Hohmann bridge, last tile, at boundary t=t_hi_last): [{p_lo_h}, {p_hi_h}]")
print(f"p-box (central mesh, fresh tile): [{p_lo_c}, {p_hi_c}]")
if not p_overlap_ok:
    print("\n*** JUNCTION FAILED: the p-boxes do not overlap ***")
    sys.exit(1)
print(f"p-overlap: [{P_overlap_lo}, {P_overlap_hi}]  (width {P_overlap_hi - P_overlap_lo})")

# Hull of the two phi-boxes: the common root, if it exists, must lie
# inside it.
hull_phi_lo = min(c_lo, h_lo)
hull_phi_hi = max(c_hi, h_hi)
hull_phi = mp.iv.mpf([hull_phi_lo, hull_phi_hi])
P_overlap_box = mp.iv.mpf([P_overlap_lo, P_overlap_hi])

# Uniqueness of phi at fixed p on the overlap: C_pphi must have
# strictly constant sign (bounded away from 0) on
# P_overlap x hull_phi. Reuses iv_Cpphi, already validated on the full
# central tile, hence valid a fortiori on this sub-box.
Cpphi_box = iv_Cpphi(P_overlap_box, hull_phi)
cpphi_lo, cpphi_hi = mp.mpf(Cpphi_box.a), mp.mpf(Cpphi_box.b)
cpphi_separated = (cpphi_lo > 0) or (cpphi_hi < 0)

print(f"C_pphi on [p-overlap] x [phi-hull]: [{cpphi_lo}, {cpphi_hi}]")
print(f"  bounded away from 0 (uniqueness of phi at fixed p): {cpphi_separated}")

if not cpphi_separated:
    print("\n*** JUNCTION FAILED: C_pphi does not have constant sign on the")
    print("    overlap box; uniqueness of phi at fixed p is not established.")
    sys.exit(1)

print("\n=== JUNCTION ESTABLISHED (full (p, phi) check) ===")
print("1. The p-boxes overlap (not only the phi-boxes):")
print(f"   p-overlap = [{P_overlap_lo}, {P_overlap_hi}].")
print("2. On this p-overlap, combined with the hull of the two phi-boxes,")
print("   C_pphi has strictly constant sign: for every p in the overlap,")
print("   Cp(p, .) = 0 has at most one root phi in the hull. Each side's own")
print("   existence certificate (Hohmann: Ghat_a = t*Cpp > 0 on the last")
print("   tile; central: interval-Newton enclosure on the fresh tile) keeps")
print("   its own root inside its own phi-box, itself contained in the")
print("   hull. Since the hull admits at most one root for each p in the")
print("   overlap, the two roots coincide: the Hohmann bridge and the")
print("   central mesh describe the same physical (p, phi) stationary point.")
