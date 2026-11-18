"""
2i-OTCA-C v1.0 -- Independent numerical cross-check (Lambert p-iteration method)
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See the project LICENSE file / ATBD Annex A for the full license text.

Independent cross-check of the exact two-impulse transfer cost between two
circular coplanar orbits, via the classical p-iteration / semi-latus-rectum
parametrization of the Lambert geometry problem -- mathematically equivalent
to the velocity-triangle formulation of Shen & Tsiotras, "Optimal Two-Impulse
Rendezvous Between Two Circular Orbits Using Multiple-Revolution Lambert's
Solutions," AAS/AIAA Astrodynamics Specialist Conference, August 1999,
Eqs. (3)-(4) and Appendix A -- but parametrized by p instead of the
semimajor axis a, which avoids their alpha/beta branch-sign ambiguity.

Independence note: this script is written from scratch from the published
equations; it does not import or call otca_c.core, and uses a different
free parameter (p) and a different root-finding/optimization method
(golden-section on an analytically bracketed admissible interval) than
whatever otca_c.core implements internally.

Derivation (short-path / single-revolution branch, theta in (0, 180) deg,
theta measured as the geometric transfer angle swept from r1 to r2):
  e*cos(f1) = p/r1 - 1
  e*sin(f1) = [(p/r1 - 1)*cos(theta) - (p/r2 - 1)] / sin(theta)
  e = hypot(e*cos f1, e*sin f1)
  Only p for which e < 1 correspond to a bound ellipse through both points
  on this branch; this admissible interval (p_lo, p_hi) is bracketed
  numerically (e(p)=1 has exactly two positive roots on the short-path
  branch) before minimizing the cost within it.
  h = sqrt(mu*p); v_r,i = (mu/h)*e*sin(f_i); v_t,i = h/r_i
  dV_i = sqrt(v_i^2 + vc_i^2 - 2*v_i*vc_i*cos(gamma_i)), gamma_i = atan2(v_r,i, v_t,i)
  dV = dV_1 + dV_2, minimized over p in (p_lo, p_hi).

mpmath (50 digits) is used throughout because e*sin(f1) involves a
subtraction that cancels catastrophically in float64 as theta -> 180 deg
(pure Hohmann limit); the quantity itself stays finite and well-defined,
only naive double-precision evaluation loses accuracy there. This mirrors
the mpmath precision-guard already used elsewhere in this project
(ATBD Annex B.2).
"""
import mpmath as mp
mp.mp.dps = 50


def _e_of_p(p, r1, r2, ct, st):
    ecf1 = p / r1 - 1
    esf1 = ((p / r1 - 1) * ct - (p / r2 - 1)) / st
    return mp.sqrt(ecf1**2 + esf1**2)


def transfer_cost_mp(p, r1, r2, theta, mu):
    ct, st = mp.cos(theta), mp.sin(theta)
    ecf1 = p / r1 - 1
    esf1 = ((p / r1 - 1) * ct - (p / r2 - 1)) / st
    e = mp.sqrt(ecf1**2 + esf1**2)
    if e >= 1:
        return mp.nan
    f1 = mp.atan2(esf1, ecf1)
    f2 = f1 + theta
    h = mp.sqrt(mu * p)

    def speed_and_fpa(r, f):
        vr = (mu / h) * e * mp.sin(f)
        vt = h / r
        v = mp.sqrt(vr**2 + vt**2)
        gamma = mp.atan2(vr, vt)
        return v, gamma

    v1, g1 = speed_and_fpa(r1, f1)
    v2, g2 = speed_and_fpa(r2, f2)
    vc1 = mp.sqrt(mu / r1)
    vc2 = mp.sqrt(mu / r2)
    dv1 = mp.sqrt(max(v1**2 + vc1**2 - 2 * v1 * vc1 * mp.cos(g1), mp.mpf(0)))
    dv2 = mp.sqrt(max(v2**2 + vc2**2 - 2 * v2 * vc2 * mp.cos(g2), mp.mpf(0)))
    return dv1 + dv2


def _bracket_admissible_interval(r1, r2, ct, st, s):
    """Find (p_lo, p_hi) with e(p)=1 by scanning then bisecting, robust for
    any r1, r2, theta on the short-path branch -- no case-specific tuning."""
    n_scan = 3000
    ps = [s * mp.mpf(k + 1) / n_scan * 2 for k in range(n_scan)]
    vals = [_e_of_p(p, r1, r2, ct, st) - 1 for p in ps]
    roots = []
    for i in range(len(ps) - 1):
        if vals[i] * vals[i + 1] < 0:
            root = mp.findroot(lambda p: _e_of_p(p, r1, r2, ct, st) - 1,
                                (ps[i], ps[i + 1]), solver='bisect')
            roots.append(root)
    if len(roots) < 2:
        raise RuntimeError(f"Could not bracket admissible p-interval "
                            f"(found {len(roots)} root(s) of e(p)=1)")
    roots.sort()
    return roots[0], roots[-1]


def _golden_section_min(f, lo, hi, tol):
    gr = (mp.sqrt(5) - 1) / 2
    c = hi - gr * (hi - lo)
    d_ = lo + gr * (hi - lo)
    fc, fd = f(c), f(d_)
    while abs(hi - lo) > tol:
        if fc < fd:
            hi, d_, fd = d_, c, fc
            c = hi - gr * (hi - lo)
            fc = f(c)
        else:
            lo, c, fc = c, d_, fd
            d_ = lo + gr * (hi - lo)
            fd = f(d_)
    xm = (lo + hi) / 2
    return xm, f(xm)


def min_cost_transfer_mp(r1, r2, theta_deg, mu):
    r1m, r2m, mum = mp.mpf(r1), mp.mpf(r2), mp.mpf(mu)
    theta = mp.radians(mp.mpf(theta_deg))
    ct, st = mp.cos(theta), mp.sin(theta)
    d = mp.sqrt(r1m**2 + r2m**2 - 2 * r1m * r2m * ct)
    s = (r1m + r2m + d) / 2
    p_lo, p_hi = _bracket_admissible_interval(r1m, r2m, ct, st, s)
    eps = (p_hi - p_lo) * mp.mpf('1e-9')
    f = lambda p: transfer_cost_mp(p, r1m, r2m, theta, mum)
    p_opt, dv_opt = _golden_section_min(f, p_lo + eps, p_hi - eps, (p_hi - p_lo) * mp.mpf('1e-30'))
    return float(p_opt), float(dv_opt), (float(p_lo), float(p_hi))


if __name__ == "__main__":
    mu = 398600.4418
    r1, r2 = 6678.0, 10017.0
    print("Sanity check vs closed-form Hohmann (theta -> 180 deg):")
    for theta in [170.0, 179.0, 179.9, 179.999]:
        p_opt, dv, (plo, phi) = min_cost_transfer_mp(r1, r2, theta, mu)
        print(f"  theta={theta:<10.6f} p_opt={p_opt:12.4f}  bracket=({plo:.2f},{phi:.2f})  dV_min={dv:.12f} km/s")
    vc1 = (mu / r1) ** 0.5
    vc2 = (mu / r2) ** 0.5
    a_h = (r1 + r2) / 2
    v1_t = (mu * (2 / r1 - 1 / a_h)) ** 0.5
    v2_t = (mu * (2 / r2 - 1 / a_h)) ** 0.5
    dv_hohmann = abs(v1_t - vc1) + abs(vc2 - v2_t)
    print(f"  Hohmann closed-form reference: dV = {dv_hohmann:.12f} km/s")
