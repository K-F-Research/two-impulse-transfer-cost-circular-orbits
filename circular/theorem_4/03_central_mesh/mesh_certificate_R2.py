# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
# Certifies, by interval arithmetic, that the curvature determinant H(p,phi)
# changes sign exactly once on the stationary branch over [w=0.15, w=1.144701]
# (w := sqrt(p)), with H<0 immediately below the sign change and H>0
# immediately above it, consistent with the proof document.
#
# The construction walks outward in both directions from a Krawczyk-certified
# interior zero K_1 of the branch's curvature, tile by tile, certifying on
# each tile: (1) existence and uniqueness of the stationary-branch root
# phi(p) via an interval-Newton contraction, (2) strict positivity of S1, S2
# (well-definedness of the cost function's square roots) on the full tile
# box, (3) the sign of H via a first-order (mean-value) Taylor enclosure
# with an interval-arithmetic remainder bound, and (4) continuity of the
# branch across tile boundaries via a nonempty-intersection test between
# consecutive tiles' certified enclosures.
#
# Four refinements over the construction's first working version, all
# confirmed without changing the method or the tolerances:
#
#   (A) mp.iv.dps (interval working precision) is set explicitly to 50
#       before any interval evaluation; only mp.mp.dps was set previously.
#
#   (B) Point values (Cp0 in the interval-Newton contraction, H0 in the
#       tile certifier, and the point values used in the Krawczyk
#       re-check) are evaluated end-to-end by dedicated interval functions
#       (sp.lambdify(..., mp.iv)) on a zero-width interval, rather than
#       computed by an ordinary-float function and wrapped in mp.iv.mpf(...)
#       afterward -- the latter does not recover rounding error already
#       lost during the floating-point computation.
#
#   (C) Tile-to-tile branch continuity is certified, not merely followed
#       numerically: after each tile is certified, its phi-enclosure B_phi
#       is checked for a nonempty intersection with the previous tile's
#       certified B_phi. Since each tile's certification already proves
#       existence and uniqueness of the root of Cp(p,.)=0 over its entire
#       p-block -- including the single boundary point shared with the
#       neighboring tile -- a nonempty intersection of the two enclosures
#       at that point forces, by uniqueness, that both tiles certify the
#       SAME root, hence the same branch. This is the rigorous tile-to-tile
#       stitching argument, obtained by a simple interval-intersection
#       test rather than an expensive recomputation.
#
#   (D) S1>0 and S2>0 (well-definedness of the square roots inside the
#       cost function C=sqrt(S1)+sqrt(S2)) are checked EXPLICITLY on every
#       tile's full box, not only implicitly via the exception that
#       mp.iv.sqrt raises on a negative radicand. A global pass/fail
#       counter is kept for the audit trail.
#
#   (E) Non-singularity of the Krawczyk preconditioner J0 is asserted and
#       printed explicitly (|det(J0)| above a conservative safety
#       threshold) before any tile computation.
#
#   (F) Beyond the floating-point diagnostic of (E), the preconditioner
#       ACTUALLY USED, Y = J0^-1, is certified non-singular by interval
#       arithmetic: det(Y) is enclosed by an interval that does not
#       contain 0. det(J0) (ordinary float) remains a diagnostic only;
#       det(Y) (interval) is what certifies non-singularity.
#
#   (G) Every extraction of an interval bound to mp.mpf (Hp_mag, Hphi_mag,
#       H0_mag_worst, and the S1/S2 positivity tests) is now performed at
#       mp.mp.dps=60, strictly above the interval working precision
#       mp.iv.dps=50. Converting mp.mpf(x) at a target precision BELOW the
#       source precision rounds to nearest and can, in the wrong
#       direction, turn a certified upper bound into a slightly smaller
#       value -- breaking the majorization/minorization guarantee. With
#       the target precision strictly above the source precision, every
#       such conversion is exact (no rounding), with no loss of rigor.
#       Sign tests are performed directly on the interval bounds themselves
#       (already the case previously), never on a floating-point midpoint.
#
#   (H) The Taylor-remainder bound (Hp_mag, Hphi_mag, remainder) is computed
#       entirely in mp.iv arithmetic (directed rounding guaranteed by
#       construction) up to the final comparison; only the final extraction
#       of the upper bound to mp.mpf is performed, and it is exact by (G).
#
#   (I) Disk checkpointing per half-mesh (DOWN walk, then UP walk): each
#       half is saved to disk as soon as it completes; if the script is
#       restarted after an interruption, it re-reads the checkpoint(s)
#       already present instead of recomputing from scratch. This affects
#       resilience to an interruption only, not the method or the
#       tolerances.

import pickle
import time
import sympy as sp
import mpmath as mp

mp.iv.dps = 50  # (A): interval precision set explicitly, never fixed before
mp.mp.dps = 60  # (G): STRICTLY above mp.iv.dps=50, so that every conversion
                # mp.mpf(interval_endpoint) is EXACT (target precision >
                # source precision => no rounding, hence no loss of rigor
                # on the majorizations/minorizations feeding the sign tests
                # and the Taylor remainder).

# ---------------------------------------------------------------------
# 1. Symbolic model (R = 2), built directly from Section 1's definitions
#    -- identical construction to krawczyk_certificate_R2.py.
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
# (B): interval evaluators for the VALUES themselves (not only their
# derivatives), used for every certified point value.
iv_Cp = sp.lambdify((p, phi), Cp, mp.iv)
iv_H = sp.lambdify((p, phi), H, mp.iv)

# ---------------------------------------------------------------------
# Two further checks, added without changing the architecture or the
# already-validated tolerances:
#
#   (D) EXPLICIT positivity of S1,S2 on each tile (Bp x B_phi), not only
#       implicit via the exception raised by mp.iv.sqrt on a negative
#       radicand inside C=sqrt(S1)+sqrt(S2). A global pass/fail counter is
#       kept for the audit trail.
#   (E) non-singularity of the Krawczyk preconditioner: explicit printout
#       of det(J0) and its magnitude, with a safety threshold.
# ---------------------------------------------------------------------
iv_S1 = sp.lambdify((p, phi), S1, mp.iv)
iv_S2 = sp.lambdify((p, phi), S2, mp.iv)

_s1s2_check_counts = {"passed": 0, "failed": 0}


def check_S1_S2_positive(Bp, B_phi):
    """(D) EXPLICIT check (not only implicit) of S1,S2 > 0 on the full box
    Bp x B_phi. Returns True/False and increments the global audit
    counter."""
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
    """Wraps a scalar mpf as a zero-width interval, ready to be passed to
    an iv_* function -- distinct from 'compute in float then wrap the
    RESULT'; here it is the INPUT that is an exact (point) interval, so
    the subsequent evaluation by iv_Cp/iv_H is an interval evaluation
    end to end."""
    return mp.iv.mpf(val)


# ---------------------------------------------------------------------
# 2. Root refinement (independent Newton solve, same method as
#    krawczyk_certificate_R2.py, starting point only, not assumed)
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
print("Refined root: p_c =", p_c, " phi_c =", phi_c, " w_c =", w_c)


def newton_phi(p0, phi_guess):
    """Ordinary high-precision Newton for Cp(p0,phi)=0, to get a good
    center before certifying the tile with interval Newton. Used ONLY to
    locate a starting guess -- never to supply a certified value (per
    fix B above)."""
    ph = phi_guess
    for _ in range(25):
        val = f_Cp(p0, ph)
        der = f_Cpphi(p0, ph)
        step = val / der
        ph = ph - step
        if abs(step) < mp.mpf('1e-40'):
            break
    return ph


# ---------------------------------------------------------------------
# 3-4. Per-tile certification and the adaptive walk
# ---------------------------------------------------------------------
def _newton_contract(Bp, dp, p0, phi_center_guess, phi_bracket_halfwidth):
    """One interval-Newton contraction attempt for the branch phi(w) at a
    fixed starting bracket half-width, using a mean-value enclosure of the
    numerator to avoid the interval-dependency blow-up of the raw,
    deeply-nested expression for Cp. Returns the final B_phi if it
    contracted strictly inside the starting bracket, else None."""
    B_phi = mp.iv.mpf([phi_center_guess - phi_bracket_halfwidth,
                        phi_center_guess + phi_bracket_halfwidth])
    B_phi_start = B_phi
    prev_width = None
    try:
        for _ in range(40):
            phi_mid = (mp.mpf(B_phi.a) + mp.mpf(B_phi.b)) / 2
            # (B): Cp0 evaluated by INTERVAL function at a point, not
            # by the ordinary-float f_Cp wrapped afterward.
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
        # the bracket strayed outside the physical domain (e.g. a sqrt of
        # a negative interval) -- treat exactly like a failed contraction,
        # the caller will try a narrower starting bracket.
        return None

    Bs_lo, Bs_hi = mp.mpf(B_phi_start.a), mp.mpf(B_phi_start.b)
    Bf_lo, Bf_hi = mp.mpf(B_phi.a), mp.mpf(B_phi.b)
    if not (Bf_lo > Bs_lo and Bf_hi < Bs_hi):
        return None
    return B_phi


def certify_tile(w_lo, w_hi, phi_center_guess, phi_bracket_halfwidth):
    """Certify the branch phi(w) and the sign of H on p in [w_lo,w_hi]^2,
    using mean-value (first-order Taylor) enclosures throughout. The right
    starting bracket half-width for the interval-Newton branch solve varies
    by orders of magnitude depending on the tile's position and width, so a
    small geometric scan picks a working one directly (each candidate only
    costs two cheap calls, Cpp/Cpphi, not the expensive Hp/Hphi, evaluated
    once at the end). Returns (ok, B_phi, H0_iv, remainder, sign); H0_iv/rem
    are still returned on remainder failure so the caller can rescale the
    tile width analytically (remainder scales ~linearly with tile width).

    H0_iv is now a genuine interval (point-width, but interval-EVALUATED --
    fix B), not a float wrapped after the fact; 'sign' is only returned
    when H0_iv does not straddle zero, i.e. the sign itself is certified."""
    Bp = mp.iv.mpf([w_lo, w_hi]) ** 2
    p0 = (mp.mpf(Bp.a) + mp.mpf(Bp.b)) / 2
    dp = Bp - mp.iv.mpf(p0)

    hw_cap = mp.mpf('0.2')  # well inside one physical branch, far from any domain edge
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

    Bf_lo, Bf_hi = mp.mpf(B_phi.a), mp.mpf(B_phi.b)

    # (D): EXPLICIT check of S1,S2 > 0 on the full box, before proceeding
    # -- in addition to the existing implicit mechanism (exception on a
    # negative sqrt, already present in _newton_contract).
    if not check_S1_S2_positive(Bp, B_phi):
        return False, B_phi, None, None, None

    phi0 = (Bf_lo + Bf_hi) / 2
    # (B): H0 evaluated by INTERVAL function at the point (p0,phi0),
    # not by ordinary-float f_H wrapped afterward.
    H0_iv = iv_H(iv_point(p0), iv_point(phi0))

    # (H): the Taylor-remainder enclosure (bounding |Hp|, |Hphi|, and the
    # remainder itself) stays ENTIRELY in mp.iv interval arithmetic
    # (directed rounding guaranteed by construction) up to the final
    # comparison. An earlier version extracted Hp_box.a/.b to mp.mpf and
    # THEN performed max()/multiplication/addition in ORDINARY floating
    # point (round-to-nearest) -- correct at the working precision used
    # here (negligible effect), but not an end-to-end directed-rounding
    # guarantee. Here abs/*/+ are all carried out by mp.iv (correctly
    # outward-rounded output); only the FINAL extraction of the upper
    # bound remainder_iv.b to mp.mpf is performed, and it is exact since
    # mp.mp.dps=60 > mp.iv.dps=50 (fix G).
    dp_half_iv = abs(dp)                       # dp = Bp - p0, already an interval
    dphi_half_iv = abs(B_phi - mp.iv.mpf(phi0))
    Hp_box = iv_Hp(Bp, B_phi)
    Hphi_box = iv_Hphi(Bp, B_phi)
    Hp_mag_iv = abs(Hp_box)
    Hphi_mag_iv = abs(Hphi_box)
    remainder_iv = Hp_mag_iv * dp_half_iv + Hphi_mag_iv * dphi_half_iv
    remainder = mp.mpf(remainder_iv.b)  # certified upper bound of the remainder
    # H0 magnitude used for the remainder test must itself be certified:
    # take the endpoint of H0_iv closest to zero (worst case), not a
    # float midpoint. H0_a/H0_b are EXACT extractions (fix G), and
    # abs()/min() on two exact values introduces no rounding.
    H0_a, H0_b = mp.mpf(H0_iv.a), mp.mpf(H0_iv.b)
    if H0_a <= 0 <= H0_b:
        # H0 itself not certified away from zero -- cannot proceed regardless
        # of the remainder; treat as a failed tile (forces rescaling).
        H0_mag_worst = mp.mpf('0')
    else:
        H0_mag_worst = min(abs(H0_a), abs(H0_b))
    if not (remainder < H0_mag_worst):
        return False, B_phi, H0_iv, remainder, None
    # sign is certified iff H0_iv itself does not straddle 0 (checked above)
    sign = 1 if H0_a > 0 else -1
    return True, B_phi, H0_iv, remainder, sign


def walk(w_start, w_end, direction, phi_guess0, w_c_local):
    """direction = +1 (increasing w) or -1 (decreasing w). Tile width is a
    persistent-adaptive fraction k_eff of the current distance to w_c:
    k_eff shrinks after a tile needed analytic rescaling (remembering the
    smaller safe scale instead of resetting to the global cap every tile,
    which would force an expensive from-scratch rescale on every single
    tile in regions where the cap is unsafe) and eases back up after an
    unscaled success.

    (C): after each tile is certified, its B_phi enclosure is checked
    for nonempty intersection with the PREVIOUS tile's certified B_phi.
    Each tile's certify_tile already proves existence+uniqueness of the
    root of Cp(p,.)=0 over its ENTIRE Bp block, including the single
    shared boundary point p with the neighboring tile -- so if the two
    tiles' B_phi enclosures overlap at all, uniqueness forces that overlap
    to contain the (unique) shared-boundary root for BOTH tiles, i.e. they
    certify the SAME branch. An empty intersection would mean the walk
    jumped to a different branch and must abort rather than silently
    continue on numerics alone."""
    w = w_start
    phi_guess = phi_guess0
    h_signs = set()
    n_tiles = 0
    h_min, h_max = None, None
    dw_max = mp.mpf('3e-3')
    dw_min = mp.mpf('1e-11')
    k_eff = mp.mpf('0.3')
    k_cap = mp.mpf('0.3')
    t0 = time.time()
    tol = mp.mpf('1e-20')
    prev_B_phi = None
    n_boundary_checks = 0
    while True:
        if direction == +1 and w >= w_end - tol:
            break
        if direction == -1 and w <= w_end + tol:
            break
        dist = abs(w - w_c_local)
        dw = min(dw_max, max(dw_min, k_eff * dist))
        if direction == +1:
            w_lo, w_hi = w, min(w + dw, w_end)
        else:
            w_lo, w_hi = max(w - dw, w_end), w
        p_mid = ((w_lo + w_hi) / 2) ** 2
        phi_center = newton_phi(p_mid, phi_guess)
        # the branch itself varies by roughly O(dw) in phi across the tile;
        # the starting bracket for the interval-Newton phi-solve must be
        # comfortably wider than that true variation.
        halfwidth = max(mp.mpf('1e-6'), 25 * dw)
        cur_scale = mp.mpf(1)
        ok, B_phi, H0, rem, sign = certify_tile(w_lo, w_hi, phi_center, halfwidth)
        attempt = 0
        while not ok and attempt < 20:
            attempt += 1
            if rem is not None and rem > 0 and H0 is not None and not (
                    mp.mpf(H0.a) <= 0 <= mp.mpf(H0.b)):
                H0_mag = min(abs(mp.mpf(H0.a)), abs(mp.mpf(H0.b)))
                cur_scale = cur_scale * (mp.mpf('0.4') * H0_mag / rem)
            else:
                cur_scale = cur_scale * mp.mpf('0.3')
            if cur_scale <= 0:
                cur_scale = mp.mpf('1e-8')
            if direction == +1:
                cur_w_hi = w_lo + (w_hi - w_lo) * cur_scale
                cur_w_lo = w_lo
            else:
                cur_w_lo = w_hi - (w_hi - w_lo) * cur_scale
                cur_w_hi = w_hi
            cur_halfwidth = max(mp.mpf('1e-6'), 25 * (cur_w_hi - cur_w_lo))
            ok, B_phi, H0, rem, sign = certify_tile(
                cur_w_lo, cur_w_hi, phi_center, cur_halfwidth)
        if ok:
            if cur_scale != 1:
                w_lo, w_hi = cur_w_lo, cur_w_hi
                k_eff = max(dw_min, k_eff * cur_scale * mp.mpf('0.8'))
            else:
                k_eff = min(k_cap, k_eff * mp.mpf('1.3'))
        if not ok:
            raise RuntimeError(f"failed to certify tile near w={w_lo},{w_hi}")

        # (C): tile-to-tile branch continuity/uniqueness check.
        if prev_B_phi is not None:
            lo = max(mp.mpf(prev_B_phi.a), mp.mpf(B_phi.a))
            hi = min(mp.mpf(prev_B_phi.b), mp.mpf(B_phi.b))
            if lo > hi:
                raise RuntimeError(
                    f"branch continuity FAILED between consecutive tiles "
                    f"near w={w_lo},{w_hi}: enclosures "
                    f"{prev_B_phi} and {B_phi} do not overlap -- "
                    f"cannot certify same branch")
        if prev_B_phi is not None:
            n_boundary_checks += 1
        prev_B_phi = B_phi

        h_signs.add(sign)
        H0_mid = (mp.mpf(H0.a) + mp.mpf(H0.b)) / 2
        if h_min is None or H0_mid < h_min:
            h_min = H0_mid
        if h_max is None or H0_mid > h_max:
            h_max = H0_mid
        n_tiles += 1
        phi_guess = (mp.mpf(B_phi.a) + mp.mpf(B_phi.b)) / 2
        if direction == +1:
            w = w_hi
        else:
            w = w_lo
    elapsed = time.time() - t0
    return n_tiles, h_signs, elapsed, h_min, h_max, n_boundary_checks


# ---------------------------------------------------------------------
# Re-verify that eps = 1e-7 still gives a valid Krawczyk box (a wider,
# still-valid choice than the eps = 1e-9 of krawczyk_certificate_R2.py),
# then walk outward from its two edges.
# ---------------------------------------------------------------------
iv_S1 = sp.lambdify((p, phi), S1, mp.iv)
iv_S2 = sp.lambdify((p, phi), S2, mp.iv)
eps = mp.mpf('1e-7')
Bp_edge = mp.iv.mpf([p_c - eps, p_c + eps])
Bphi_edge = mp.iv.mpf([phi_c - eps, phi_c + eps])
J0 = Jmat(z)
# (E): EXPLICIT non-singularity of the Krawczyk preconditioner.
det_J0 = J0[0, 0]*J0[1, 1] - J0[0, 1]*J0[1, 0]
print(f"\nKrawczyk preconditioner J0 = {J0}")
print(f"det(J0) = {det_J0}  |det(J0)| = {abs(det_J0)}")
DET_J0_MIN = mp.mpf('1e-6')  # safety threshold, very conservative relative
                              # to the observed magnitude (~2e4)
assert abs(det_J0) > DET_J0_MIN, (
    f"Krawczyk preconditioner close to singular: "
    f"|det(J0)|={abs(det_J0)} <= {DET_J0_MIN}")
print(f"Krawczyk preconditioner non-singular (|det(J0)| > {DET_J0_MIN}): PASS "
      f"(floating-point diagnostic only -- see the interval certificate below)")
Yinv = J0**-1
J_box = mp.iv.matrix([[iv_Cpp(Bp_edge, Bphi_edge), iv_Cpphi(Bp_edge, Bphi_edge)],
                       [iv_Hp(Bp_edge, Bphi_edge), iv_Hphi(Bp_edge, Bphi_edge)]])
z0 = mp.iv.matrix([p_c, phi_c])
box = mp.iv.matrix([Bp_edge, Bphi_edge])
Yiv = mp.iv.matrix(Yinv)

# (F): interval certificate of non-singularity for the preconditioner
# ACTUALLY USED, Y = J0^-1 (not only J0 itself in ordinary floating
# point). det(J0) above remains a diagnostic; it is det(Yiv), enclosed by
# an interval and strictly separated from 0, that constitutes the
# non-singularity certificate.
det_Yiv = Yiv[0, 0] * Yiv[1, 1] - Yiv[0, 1] * Yiv[1, 0]
print(f"det(Y) interval = {det_Yiv}  (Y = J0^-1, the preconditioner used)")
Yiv_nonsingular = (det_Yiv.a > 0) or (det_Yiv.b < 0)
assert Yiv_nonsingular, (
    "Non-singularity of the preconditioner Y=J0^-1 NOT certified by interval "
    f"(det(Y) interval = {det_Yiv} contains or touches 0)")
print(f"Preconditioner Y=J0^-1 non-singular (interval certificate, "
      f"det(Y) strictly separated from 0): PASS")

I2 = mp.iv.eye(2)
# (B): F0iv evaluated by INTERVAL functions at the point z, not by
# wrapping the ordinary-float F(z) result afterward.
F0iv = mp.iv.matrix([iv_Cp(iv_point(z[0]), iv_point(z[1])),
                      iv_H(iv_point(z[0]), iv_point(z[1]))])
K = (z0 - Yiv * F0iv) + (I2 - Yiv * J_box) * (box - z0)
krawczyk_ok = (K[0].a > box[0].a) and (K[0].b < box[0].b) and \
              (K[1].a > box[1].a) and (K[1].b < box[1].b)
print(f"\nKrawczyk re-check at eps=1e-7: {krawczyk_ok}")
assert krawczyk_ok, "eps=1e-7 box does not validate -- cannot use as mesh starting edge"

w_edge_lo = mp.sqrt(p_c - eps)
w_edge_hi = mp.sqrt(p_c + eps)
phi_at_lo = newton_phi(w_edge_lo ** 2, phi_c)
phi_at_hi = newton_phi(w_edge_hi ** 2, phi_c)

print("K_1 (rebuilt, eps=1e-7) edges: w_lo=%s  w_hi=%s" % (w_edge_lo, w_edge_hi))

# ---------------------------------------------------------------------
# (I): checkpointing. Each half-mesh (DOWN, then UP) is saved to disk as
# soon as it completes; if the script restarts after an interruption, it
# re-reads
# the checkpoint(s) already present instead of recomputing everything.
# The global audit counters (_s1s2_check_counts) are also restored to
# their exact state at the end of the DOWN walk, so the final count
# (DOWN+UP) stays correct even after a resumption.
# ---------------------------------------------------------------------
import os
CKPT_DOWN = '/tmp/mesh_certificate_R2_checkpoint_down.pkl'
CKPT_UP = '/tmp/mesh_certificate_R2_checkpoint_up.pkl'

if os.path.exists(CKPT_DOWN):
    with open(CKPT_DOWN, 'rb') as fh:
        ck = pickle.load(fh)
    n_down, t_down, hmin_down, hmax_down, nb_down = (
        ck['n_down'], ck['t_down'], ck['hmin_down'], ck['hmax_down'], ck['nb_down'])
    signs_down = {ck['sign_down']}
    _s1s2_check_counts['passed'] = ck['s1s2_passed']
    _s1s2_check_counts['failed'] = ck['s1s2_failed']
    print(f"\n[RESUMED] DOWN walk re-read from checkpoint {CKPT_DOWN} "
          f"(avoids a full recomputation after an interruption): tiles={n_down}, "
          f"signs={signs_down}, boundary_checks={nb_down}")
else:
    print("\nWalking DOWN from K_1 edge to w=0.15 ...")
    n_down, signs_down, t_down, hmin_down, hmax_down, nb_down = walk(
        w_edge_lo, mp.mpf('0.15'), -1, phi_at_lo, w_c)
    print(f"  tiles={n_down}, signs={signs_down}, time={t_down:.1f}s, "
          f"H range=[{mp.nstr(hmin_down, 6)},{mp.nstr(hmax_down, 6)}], "
          f"boundary_checks={nb_down}")
    with open(CKPT_DOWN, 'wb') as fh:
        pickle.dump({
            'n_down': n_down, 'sign_down': list(signs_down)[0],
            't_down': t_down, 'hmin_down': hmin_down, 'hmax_down': hmax_down,
            'nb_down': nb_down,
            's1s2_passed': _s1s2_check_counts['passed'],
            's1s2_failed': _s1s2_check_counts['failed'],
        }, fh)
    print(f"[CHECKPOINT] DOWN walk saved to {CKPT_DOWN}")

if os.path.exists(CKPT_UP):
    with open(CKPT_UP, 'rb') as fh:
        ck = pickle.load(fh)
    n_up, t_up, hmin_up, hmax_up, nb_up = (
        ck['n_up'], ck['t_up'], ck['hmin_up'], ck['hmax_up'], ck['nb_up'])
    signs_up = {ck['sign_up']}
    _s1s2_check_counts['passed'] = ck['s1s2_passed']
    _s1s2_check_counts['failed'] = ck['s1s2_failed']
    print(f"\n[RESUMED] UP walk re-read from checkpoint {CKPT_UP} "
          f"(avoids a full recomputation after an interruption): tiles={n_up}, "
          f"signs={signs_up}, boundary_checks={nb_up}")
else:
    print("\nWalking UP from K_1 edge to w=1.144701 ...")
    n_up, signs_up, t_up, hmin_up, hmax_up, nb_up = walk(
        w_edge_hi, mp.mpf('1.144701'), +1, phi_at_hi, w_c)
    print(f"  tiles={n_up}, signs={signs_up}, time={t_up:.1f}s, "
          f"H range=[{mp.nstr(hmin_up, 6)},{mp.nstr(hmax_up, 6)}], "
          f"boundary_checks={nb_up}")
    with open(CKPT_UP, 'wb') as fh:
        pickle.dump({
            'n_up': n_up, 'sign_up': list(signs_up)[0],
            't_up': t_up, 'hmin_up': hmin_up, 'hmax_up': hmax_up,
            'nb_up': nb_up,
            's1s2_passed': _s1s2_check_counts['passed'],
            's1s2_failed': _s1s2_check_counts['failed'],
        }, fh)
    print(f"[CHECKPOINT] UP walk saved to {CKPT_UP}")

print("\n=== SUMMARY ===")
assert len(signs_down) == 1 and len(signs_up) == 1
assert signs_down != signs_up
print("Exactly one sign change confirmed (inside K_1), consistent with the proof document.")
print(f"Tile-to-tile branch continuity certified throughout both walks (fix C): "
      f"{nb_down + nb_up} boundary checks, all passed (else RuntimeError would have aborted).")
print("All point-values (Cp0, H0, F0) interval-evaluated end-to-end (fix B).")
print(f"mp.iv.dps = {mp.iv.dps} (fix A).")
print(f"S1,S2 > 0 explicitly checked on EVERY tile (fix D): "
      f"{_s1s2_check_counts['passed']} passed, {_s1s2_check_counts['failed']} failed "
      f"(failed tiles were rescaled/retried, not silently accepted).")
print(f"Krawczyk preconditioner det(J0) = {det_J0}, |det(J0)| = {abs(det_J0)} "
      f"(fix E v3, floating-point diagnostic only).")
print(f"Krawczyk preconditioner det(Y) interval = {det_Yiv} "
      f"(fix F v4, interval certificate of non-singularity for the "
      f"preconditioner actually used Y=J0^-1): PASS.")
print(f"mp.mp.dps = {mp.mp.dps} > mp.iv.dps = {mp.iv.dps} "
      f"(fix G v4: mp.mpf(endpoint) conversions without loss of precision).")
print("Taylor remainder (Hp_mag, Hphi_mag, remainder) computed entirely in "
      "mp.iv directed-rounding arithmetic, exact final extraction "
      "(fix H v5).")
print("Disk checkpointing per half-mesh (DOWN then UP), automatically "
      "re-read on restart (fix I v6).")

with open('mesh_certificate_R2_v6_result.pkl', 'wb') as fh:
    pickle.dump({
        'n_down': n_down, 'sign_down': list(signs_down)[0],
        'n_up': n_up, 'sign_up': list(signs_up)[0],
        'hmin_down': mp.nstr(hmin_down, 10), 'hmax_down': mp.nstr(hmax_down, 10),
        'hmin_up': mp.nstr(hmin_up, 10), 'hmax_up': mp.nstr(hmax_up, 10),
        't_down': t_down, 't_up': t_up,
        'n_boundary_checks_down': nb_down, 'n_boundary_checks_up': nb_up,
        's1s2_checks_passed': _s1s2_check_counts['passed'],
        's1s2_checks_failed': _s1s2_check_counts['failed'],
        'det_J0': mp.nstr(det_J0, 20),
        'det_Yiv': str(det_Yiv),
        'mp_mp_dps': mp.mp.dps, 'mp_iv_dps': mp.iv.dps,
    }, fh)
print("\nResult saved to mesh_certificate_R2_v6_result.pkl")
