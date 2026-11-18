/*
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Pure-computation JavaScript port of src/otca_c/ (core.py, hohmann.py,
geometry.py, derived.py, lambert.py), kept in its own file (not inlined in
index.html) so the exact same code can run both in the browser (loaded via
<script src="otca.js">) and under Node for the automated parity check
against the Python library (../tests/test_js_python_parity.py).

This is a faithful, function-by-function port, not an independent
reimplementation -- every function below names its Python counterpart in
a comment. Two places deliberately deviate from a literal line-by-line
port, both documented at the point of deviation:

  1. core.py's C_p (dC/dp) is derived symbolically with SymPy; there is no
     symbolic engine here, so costDerivative() below uses a tight central
     finite difference instead. The root-finding STRATEGY (adaptive
     bracket expansion, sign-change verification, a bracketing
     root-finder) is otherwise identical to find_p_star.
  2. JavaScript has no scipy.optimize.brentq, so brentq() below is a
     from-scratch implementation of Brent's method (bisection + secant +
     inverse quadratic interpolation, the same algorithm scipy's own
     brentq implements), used both for p* and for the Lambert z root.

Everything else -- the adaptive bracket-expansion patterns, the explicit
ConvergenceError-style failures (never a silent wrong answer, never a
value clamped onto a boundary), the corrected Lambert A-formula, the
explicit ellipticity/ratio checks on the Lambert solution -- mirrors the
Python source exactly, including the same reasoning in the comments,
because that reasoning is what the independent Lambert cross-check and
RK4 validation in the Python library caught as real, fixed bugs.
*/
(function (root, factory) {
  if (typeof module !== "undefined" && module.exports) {
    module.exports = factory();
  } else {
    root.Otca = factory();
  }
})(typeof self !== "undefined" ? self : this, function () {
"use strict";

// --- Errors -----------------------------------------------------------
// Mirrors otca_c.errors.ConvergenceError: raised whenever a numerical
// search fails to locate what the underlying theorems guarantee exists,
// rather than returning a silently wrong or non-finite value.
class ConvergenceError extends Error {
  constructor(message) {
    super(message);
    this.name = "ConvergenceError";
  }
}

// --- Generic bracketing root-finder (Brent's method) -------------------
// Used both by findPStar (root of costDerivative) and solveLambert (root
// of tof(z)-tofS). Mirrors scipy.optimize.brentq's algorithm (bisection
// safeguarding secant/inverse-quadratic steps), since no such routine is
// built into JavaScript.
function brentq(f, a, b, xtol, maxiter) {
  xtol = xtol === undefined ? 1e-12 : xtol;
  maxiter = maxiter === undefined ? 100 : maxiter;
  let fa = f(a), fb = f(b);
  if (fa === 0) return { root: a, iterations: 0, converged: true };
  if (fb === 0) return { root: b, iterations: 0, converged: true };
  if (fa * fb > 0) {
    throw new ConvergenceError(
      "brentq: f(a) and f(b) must have opposite signs (a=" + a + ", b=" + b +
      ", f(a)=" + fa + ", f(b)=" + fb + ")"
    );
  }
  let c = a, fc = fa, d = b - a, e = d;
  for (let iter = 0; iter < maxiter; iter++) {
    if (Math.abs(fc) < Math.abs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    const tol1 = 2 * Number.EPSILON * Math.abs(b) + xtol / 2;
    const xm = (c - b) / 2;
    if (Math.abs(xm) <= tol1 || fb === 0) {
      return { root: b, iterations: iter, converged: true };
    }
    if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
      const s = fb / fa;
      let p, q;
      if (a === c) {
        p = 2 * xm * s;
        q = 1 - s;
      } else {
        const q1 = fa / fc, r = fb / fc;
        p = s * (2 * xm * q1 * (q1 - r) - (b - a) * (r - 1));
        q = (q1 - 1) * (r - 1) * (s - 1);
      }
      if (p > 0) q = -q; else p = -p;
      if (2 * p < Math.min(3 * xm * q - Math.abs(tol1 * q), Math.abs(e * q))) {
        e = d; d = p / q;
      } else {
        d = xm; e = d;
      }
    } else {
      d = xm; e = d;
    }
    a = b; fa = fb;
    if (Math.abs(d) > tol1) b += d; else b += (xm > 0 ? tol1 : -tol1);
    fb = f(b);
    if ((fb > 0) === (fc > 0)) { c = a; fc = fa; e = d = b - a; }
  }
  return { root: b, iterations: maxiter, converged: false };
}

// --- core.py: S1/S2/C and the p* root-finder ---------------------------

function costS1S2(p, phi, R) {
  const c = Math.cos(phi), s = Math.sin(phi);
  const x = p - 1, y = p / R - 1;
  const e2 = x * x + Math.pow(x * c + y, 2) / (s * s);
  const v1_sq = 2 + (e2 - 1) / p;
  const v2_sq = 2 / R + (e2 - 1) / p;
  const v_theta1 = Math.sqrt(p);
  const v_theta2 = Math.sqrt(p) / R;
  const v_c1 = 1, v_c2 = 1 / Math.sqrt(R);
  const S1 = v1_sq + v_c1 * v_c1 - 2 * v_c1 * v_theta1;
  const S2 = v2_sq + v_c2 * v_c2 - 2 * v_c2 * v_theta2;
  return { S1: S1, S2: S2, e2: e2 };
}

function cost(p, phi, R) {
  const r = costS1S2(p, phi, R);
  return Math.sqrt(r.S1) + Math.sqrt(r.S2);
}

// costDerivative(p,phi,R) ~ core.py's cost_derivative (_f_Cp), by central
// finite difference -- see module docstring, deviation (1).
function costDerivative(p, phi, R) {
  const h = Math.max(Math.abs(p) * 1e-6, 1e-12);
  const plus = cost(p + h, phi, R);
  const minus = Math.abs(p - h) > 0 ? cost(p - h, phi, R) : NaN;
  return (plus - minus) / (2 * h);
}

// findPStar(phi,R): adaptive-bracket root search on costDerivative,
// mirroring core.py's find_p_star exactly (contract p_lo toward 0,
// expand p_hi toward +infinity, verify the existence/uniqueness-result-
// guaranteed sign
// pattern at each candidate endpoint, never accept a non-finite
// endpoint) -- not the earlier golden-section search this file used
// before, which (like core.py's own earlier fixed-bracket find_p_star)
// could silently fail for delta_theta extremely close to 0/360 deg.
function findPStar(phi, R, pLo, pHi, xtol, maxExpansions) {
  pLo = pLo === undefined ? 1e-9 : pLo;
  pHi = pHi === undefined ? 1e9 : pHi;
  xtol = xtol === undefined ? 1e-13 : xtol;
  maxExpansions = maxExpansions === undefined ? 300 : maxExpansions;
  const hardFloor = 1e-300, hardCeiling = 1e300;

  function finiteSign(p) {
    const v = costDerivative(p, phi, R);
    return Number.isFinite(v) ? v : null;
  }

  let gLo = finiteSign(pLo), expansionsLo = 0;
  while ((gLo === null || gLo >= 0) && pLo > hardFloor && expansionsLo < maxExpansions) {
    pLo /= 10.0;
    gLo = finiteSign(pLo);
    expansionsLo++;
  }
  if (gLo === null || gLo >= 0) {
    throw new ConvergenceError(
      "could not find p_lo with a finite, negative C_p (the " +
      "existence/uniqueness result's p->0+ " +
      "limit) for phi=" + phi + ", R=" + R + " after " + expansionsLo +
      " contractions down to p_lo=" + pLo + " (hard floor " + hardFloor +
      "); this input may be on or past a genuine domain boundary."
    );
  }

  let gHi = finiteSign(pHi), expansionsHi = 0;
  while ((gHi === null || gHi <= 0) && pHi < hardCeiling && expansionsHi < maxExpansions) {
    pHi *= 10.0;
    gHi = finiteSign(pHi);
    expansionsHi++;
  }
  if (gHi === null || gHi <= 0) {
    throw new ConvergenceError(
      "could not find p_hi with a finite, positive C_p (the " +
      "existence/uniqueness result's " +
      "p->+infinity limit) for phi=" + phi + ", R=" + R + " after " +
      expansionsHi + " expansions up to p_hi=" + pHi + " (hard ceiling " +
      hardCeiling + "); this input may be on or past a genuine domain " +
      "boundary."
    );
  }

  const result = brentq(function (p) { return costDerivative(p, phi, R); }, pLo, pHi, xtol, 200);
  return result.root;
}

// solveHohmann(R): exact closed form, hohmann.py.
function solveHohmann(R) {
  const a = (1 + R) / 2;
  const e = (R - 1) / (R + 1);
  const p = a * (1 - e * e);
  const vc1 = 1, vc2 = 1 / Math.sqrt(R);
  const vPeri = Math.sqrt(2 / 1 - 1 / a);
  const vApo = Math.sqrt(2 / R - 1 / a);
  const dv1 = Math.abs(vPeri - vc1);
  const dv2 = Math.abs(vc2 - vApo);
  return { p: p, e: e, a: a, dv1: dv1, dv2: dv2, total: dv1 + dv2 };
}

// geometryCheck(p,phi,R,e,S1,S2): geometry.py's independent_check.
// IMPORTANT: `phi` must be the SIGNED model parameter (can be negative
// for delta_theta > 180 deg) -- see api.py's own comment on this, a real
// bug found and fixed during this project's development.
function geometryCheck(p, phi, R, e, S1, S2) {
  const c = Math.cos(phi), s = Math.sin(phi);
  const deltaTheta = Math.PI - phi;
  const x = p - 1, y = p / R - 1;
  const eSinNu1 = -(x * c + y) / s;
  const nu1 = Math.atan2(eSinNu1, x);
  const nu2 = nu1 + deltaTheta;
  const r1c = p / (1 + e * Math.cos(nu1));
  const r2c = p / (1 + e * Math.cos(nu2));
  const a = p / (1 - e * e);
  const h = Math.sqrt(p);
  function speedComponents(r, nu) {
    const vTheta = h / r;
    const v2total = 2 / r - 1 / a;
    const vRadialSq = Math.max(v2total - vTheta * vTheta, 0);
    // copysign(sqrt(vRadialSq), sin(nu)) -- at nu=0/pi exactly, vRadialSq
    // is itself exactly 0, so the sign choice there is immaterial; the
    // `|| 1` only avoids Math.sign(0)*0 ambiguity, it never changes a
    // nonzero result.
    const vRadial = Math.sign(Math.sin(nu) || 1) * Math.sqrt(vRadialSq);
    return [vTheta, vRadial];
  }
  const c1 = speedComponents(r1c, nu1), c2 = speedComponents(r2c, nu2);
  const vc1 = 1, vc2 = 1 / Math.sqrt(R);
  const dv1sq = Math.pow(c1[0] - vc1, 2) + c1[1] * c1[1];
  const dv2sq = Math.pow(c2[0] - vc2, 2) + c2[1] * c2[1];
  return {
    r1res: r1c - 1, r2res: r2c - R,
    dv1res: dv1sq - S1, dv2res: dv2sq - S2,
    nu1: nu1, nu2: nu2,
  };
}

function circularSpeedKmS(rKm, muKm3s2) { return Math.sqrt(muKm3s2 / rKm); }

// --- derived.py ---------------------------------------------------------

const ECCENTRICITY_CONDITIONING_THRESHOLD = 1e-6;

function eccentricityConditioningWarning(e) {
  if (!(e >= 0.0 && e < 1.0)) return null;
  if (1.0 - e < ECCENTRICITY_CONDITIONING_THRESHOLD) {
    return (
      "eccentricity e=" + e + " is extremely close to 1 (1-e=" +
      (1.0 - e).toExponential(3) + ", below this API's numerical " +
      "conditioning threshold of " + ECCENTRICITY_CONDITIONING_THRESHOLD +
      "): r_apoapsis_km, semi_major_axis_km, orbital_period_s and " +
      "time_of_flight_s divide by (1-e) or (1-e^2) and can carry a much " +
      "larger relative error here than elsewhere in the domain, even " +
      "though the underlying orbit is a genuine, physically valid " +
      "ellipse. delta_v1_km_s, delta_v2_km_s, delta_v_total_km_s and " +
      "r_periapsis_km do not divide by (1-e) and are not affected by " +
      "this warning."
    );
  }
  return null;
}

function periapsisRadiusKm(pKm, e) { return pKm / (1 + e); }

function apoapsisRadiusKm(pKm, e) {
  if (e >= 1.0) {
    throw new Error("apoapsis radius undefined for e=" + e + " >= 1 (not an ellipse)");
  }
  return pKm / (1 - e);
}

function orbitalPeriodS(aKm, muKm3s2) {
  if (!(aKm > 0)) throw new Error("semi-major axis must be positive, got " + aKm + " km");
  return 2 * Math.PI * Math.sqrt(Math.pow(aKm, 3) / muKm3s2);
}

function trueToEccentric(nu, e) {
  return 2 * Math.atan2(Math.sqrt(1 - e) * Math.sin(nu / 2), Math.sqrt(1 + e) * Math.cos(nu / 2));
}

function eccentricToTrue(E, e) {
  return 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(E / 2), Math.sqrt(1 - e) * Math.cos(E / 2));
}

// meanToEccentric: solves Kepler's equation M = E - e*sin(E) for E, given
// M (mean anomaly, the one quantity that advances at a CONSTANT rate in
// time -- n*t -- unlike true or eccentric anomaly, which is exactly why
// animating a transfer at its true physical speed needs this inverse
// step). Not a port of any otca_c module (the Python library only ever
// goes anomaly -> time, via time_of_flight_s/_true_to_eccentric, never
// the reverse) -- added here for the demo's own real-time animation.
//
// A first version of this function used plain Newton iteration from the
// initial guess E0=M. That DIVERGES catastrophically at high
// eccentricity near periapsis: the derivative dM/dE = 1-e*cos(E) goes to
// 1-e (near zero for e close to 1) right where M itself is small, so a
// Newton step there is divided by a near-zero number and overshoots
// wildly (confirmed: e=0.999, a case near periapsis sent E from ~0.04 to
// >1e10 instead of converging). Kepler's equation IS globally monotonic
// in E for every 0<=e<1 (dM/dE = 1-e*cos(E) >= 1-e > 0 always), and M
// stays within (-pi,pi] whenever E does -- so this is solved the same
// way as every other root-find in this file: brentq on a bracket that
// is guaranteed valid by that monotonicity, never a bare Newton step.
function meanToEccentric(M, e) {
  const twoPi = 2 * Math.PI;
  let Mf = M % twoPi;
  if (Mf > Math.PI) Mf -= twoPi;
  if (Mf < -Math.PI) Mf += twoPi;
  const f = function (E) { return E - e * Math.sin(E) - Mf; };
  const margin = 1e-9;
  const result = brentq(f, -Math.PI - margin, Math.PI + margin, 1e-13, 100);
  return result.root + (M - Mf);
}

function timeOfFlightS(aKm, e, nu1, nu2, muKm3s2) {
  if (!(e >= 0.0 && e < 1.0)) {
    throw new Error("timeOfFlightS requires an ellipse, 0<=e<1; got e=" + e);
  }
  const n = Math.sqrt(muKm3s2 / Math.pow(aKm, 3));
  function meanAnomalyUnwrapped(nu) {
    const k = Math.floor((nu + Math.PI) / (2 * Math.PI));
    const nuFolded = nu - k * 2 * Math.PI;
    const E = trueToEccentric(nuFolded, e);
    const M = E - e * Math.sin(E);
    return M + k * 2 * Math.PI;
  }
  const M1 = meanAnomalyUnwrapped(nu1), M2 = meanAnomalyUnwrapped(nu2);
  return (M2 - M1) / n;
}

function propellantMassFraction(deltaVKmS, ispS) {
  if (!(ispS > 0)) throw new Error("ispS must be positive, got " + ispS);
  const g0 = 9.80665;
  return 1 - Math.exp(-(deltaVKmS * 1000) / (ispS * g0));
}

// --- lambert.py -----------------------------------------------------------

function stumpffC(z) {
  if (Math.abs(z) < 1e-6) return 1 / 2 - z / 24 + z * z / 720;
  if (z > 0) { const sq = Math.sqrt(z); return (1 - Math.cos(sq)) / z; }
  const sq = Math.sqrt(-z); return (1 - Math.cosh(sq)) / z;
}

function stumpffS(z) {
  if (Math.abs(z) < 1e-6) return 1 / 6 - z / 120 + z * z / 5040;
  if (z > 0) { const sq = Math.sqrt(z); return (sq - Math.sin(sq)) / Math.pow(sq, 3); }
  const sq = Math.sqrt(-z); return (Math.sinh(sq) - sq) / Math.pow(sq, 3);
}

// solveLambert: full port of lambert.py's solve_lambert, including the
// corrected A formula, the adaptive clipped-expansion bracket search
// (replacing the earlier Newton+clamp approach this file used before),
// and the explicit ellipticity/ratio checks on the recovered orbit (no
// silent e=0.0 fallback).
function solveLambert(r1VecKm, r2VecKm, tofS, muKm3s2, prograde, maxIter, tolS) {
  prograde = prograde === undefined ? true : prograde;
  maxIter = maxIter === undefined ? 100 : maxIter;
  tolS = tolS === undefined ? 1e-6 : tolS;

  const r1 = Math.hypot(r1VecKm[0], r1VecKm[1]);
  const r2 = Math.hypot(r2VecKm[0], r2VecKm[1]);
  const crossZ = r1VecKm[0] * r2VecKm[1] - r1VecKm[1] * r2VecKm[0];
  let cosDnu = (r1VecKm[0] * r2VecKm[0] + r1VecKm[1] * r2VecKm[1]) / (r1 * r2);
  cosDnu = Math.max(-1, Math.min(1, cosDnu));
  let dnu = Math.acos(cosDnu);
  if (prograde) { if (crossZ < 0) dnu = 2 * Math.PI - dnu; }
  else { if (crossZ > 0) dnu = 2 * Math.PI - dnu; }

  if (Math.abs(crossZ) < 1e-12 * r1 * r2 && cosDnu > 0) {
    throw new ConvergenceError(
      "Lambert solver: r1_vec and r2_vec are colinear and same-pointing " +
      "(dnu=" + (dnu * 180 / Math.PI).toFixed(6) + " deg ~ 0 or 360 deg) -- " +
      "the transfer plane is geometrically undefined; this degenerate " +
      "zero-angle case is already excluded upstream by compute_transfer's " +
      "own delta_theta normalization."
    );
  }

  // A = sqrt(2*r1*r2)*cos(dnu/2) -- the corrected, numerically stable
  // form (see lambert.py's module docstring for the full derivation and
  // the sign-error this corrects for).
  const A = Math.sqrt(2.0 * r1 * r2) * Math.cos(dnu / 2.0);
  if (Math.abs(A) < 1e-9 * Math.sqrt(r1 * r2)) {
    throw new ConvergenceError(
      "Lambert solver: A~0 (dnu=" + (dnu * 180 / Math.PI).toFixed(6) +
      " deg is 0, 180, or 360 deg -- a degenerate or exactly-Hohmann " +
      "transfer plane); this case is excluded upstream by " +
      "compute_transfer's own delta_theta normalization and its exact " +
      "Hohmann closed form."
    );
  }

  const sqrtMu = Math.sqrt(muKm3s2);
  const zLoHard = -4.0 * Math.PI * Math.PI;
  const zHi = 4.0 * Math.PI * Math.PI;

  function tofOf(z) {
    if (!(zLoHard < z && z < zHi)) return null;
    const C = stumpffC(z);
    if (C <= 0.0) return null;
    const S = stumpffS(z);
    const y = r1 + r2 + A * (z * S - 1.0) / Math.sqrt(C);
    if (y < 0) return null;
    const chi = Math.sqrt(y / C);
    return (Math.pow(chi, 3) * S + A * Math.sqrt(y)) / sqrtMu;
  }

  const f0 = tofOf(0.0);
  if (f0 === null) {
    throw new ConvergenceError(
      "Lambert solver: z=0 (parabolic) is itself invalid for this " +
      "geometry (r1=" + r1.toFixed(3) + " km, r2=" + r2.toFixed(3) +
      " km, dnu=" + (dnu * 180 / Math.PI).toFixed(3) + " deg) -- cannot " +
      "even start the bracket search; please report this as a bug."
    );
  }
  const target = tofS;

  function signed(z) {
    const t = tofOf(z);
    return t === null ? null : t - target;
  }

  const fAt0 = f0 - target;
  let zLo, zHiBracket;
  if (Math.abs(fAt0) <= tolS) {
    zLo = 0.0; zHiBracket = 0.0;
  } else {
    const directionSign = fAt0 < 0 ? 1.0 : -1.0;
    const zBound = directionSign > 0 ? zHi : zLoHard;
    const margin = zBound !== 0 ? 1e-9 * Math.abs(zBound) : 1e-9;
    let step = 1.0;
    let zPrev = 0.0, fPrev = fAt0;
    zLo = null; zHiBracket = null;
    for (let i = 0; i < 200; i++) {
      const zCandidate = directionSign * step;
      const zNext = Math.max(zLoHard + margin, Math.min(zBound - margin, zCandidate));
      if (zNext === zPrev) break;
      const fNext = signed(zNext);
      if (fNext !== null && (fPrev === null || fNext * fPrev <= 0)) {
        zLo = Math.min(zPrev, zNext);
        zHiBracket = Math.max(zPrev, zNext);
        break;
      }
      if (fNext !== null) { zPrev = zNext; fPrev = fNext; }
      step *= 2.0;
    }
    if (zLo === null) {
      throw new ConvergenceError(
        "Lambert solver: could not bracket a root for tof_s=" + tofS +
        " s (r1=" + r1.toFixed(3) + " km, r2=" + r2.toFixed(3) + " km, " +
        "dnu=" + (dnu * 180 / Math.PI).toFixed(3) + " deg) within the " +
        "valid domain of z -- this time of flight may be outside what a " +
        "single-revolution transfer can achieve for this geometry."
      );
    }
  }

  function f(z) {
    const t = tofOf(z);
    if (t === null) {
      throw new ConvergenceError(
        "Lambert solver: internal bracket [" + zLo + ", " + zHiBracket +
        "] contains an invalid interior point at z=" + z + " even though " +
        "both of its own endpoints were valid -- tof(z) is not " +
        "well-behaved on this bracket; please report this as a bug rather " +
        "than trust the result."
      );
    }
    return t - target;
  }

  let zRoot, iterations;
  if (zLo === zHiBracket) {
    zRoot = zLo; iterations = 0;
  } else {
    const result = brentq(f, zLo, zHiBracket, 1e-12, maxIter);
    if (!result.converged) {
      throw new ConvergenceError(
        "Lambert solver's bracketing root-finder did not converge within " +
        maxIter + " iterations (r1=" + r1.toFixed(3) + " km, r2=" +
        r2.toFixed(3) + " km, dnu=" + (dnu * 180 / Math.PI).toFixed(3) +
        " deg, tof_s=" + tofS + ")."
      );
    }
    zRoot = result.root; iterations = result.iterations;
  }

  const C = stumpffC(zRoot), S = stumpffS(zRoot);
  const y = r1 + r2 + A * (zRoot * S - 1.0) / Math.sqrt(C);

  const fCoef = 1.0 - y / r1;
  const g = A * Math.sqrt(y / muKm3s2);
  const gdot = 1.0 - y / r2;

  const v1x = (r2VecKm[0] - fCoef * r1VecKm[0]) / g;
  const v1y = (r2VecKm[1] - fCoef * r1VecKm[1]) / g;
  const v2x = (gdot * r2VecKm[0] - r1VecKm[0]) / g;
  const v2y = (gdot * r2VecKm[1] - r1VecKm[1]) / g;

  const h = r1VecKm[0] * v1y - r1VecKm[1] * v1x;
  const pKm = h * h / muKm3s2;
  const v1Sq = v1x * v1x + v1y * v1y;
  const denom = 2.0 / r1 - v1Sq / muKm3s2;
  if (!(denom > 0)) {
    throw new ConvergenceError(
      "Lambert solver: vis-viva denominator (2/r1 - v1^2/mu) = " + denom +
      " is not positive -- the recovered velocity implies a " +
      "non-elliptical (parabolic or hyperbolic) orbit, outside this " +
      "package's always-elliptical domain; refusing to silently report " +
      "a>0 or e<1 in that case."
    );
  }
  const aKm = 1.0 / denom;
  let ratio = pKm / aKm;
  if (ratio > 1.0) {
    if (ratio - 1.0 > 1e-10) {
      throw new ConvergenceError(
        "Lambert solver: p/a = " + ratio + " exceeds 1.0 by more than the " +
        "tolerated floating-point roundoff (1e-10) -- this is not a " +
        "circular orbit artifact, it indicates an actual inconsistency in " +
        "the recovered solution; refusing to silently report e=0."
      );
    }
    ratio = 1.0;
  }
  const e = Math.sqrt(1.0 - ratio);

  return {
    pKm: pKm, aKm: aKm, e: e,
    v1: [v1x, v1y], v2: [v2x, v2y],
    z: zRoot, iterations: iterations,
  };
}

function lambertCrossCheck(r1Km, r2Km, deltaThetaRad, tofS, muKm3s2, pStarKm, e, dv1Km, dv2Km) {
  const r1vec = [r1Km, 0];
  const r2vec = [r2Km * Math.cos(deltaThetaRad), r2Km * Math.sin(deltaThetaRad)];
  const sol = solveLambert(r1vec, r2vec, tofS, muKm3s2, true);
  const vc1 = Math.sqrt(muKm3s2 / r1Km), vc2 = Math.sqrt(muKm3s2 / r2Km);
  const vc1vec = [0, vc1];
  const vc2vec = [-Math.sin(deltaThetaRad) * vc2, Math.cos(deltaThetaRad) * vc2];
  const lambertDv1 = Math.hypot(sol.v1[0] - vc1vec[0], sol.v1[1] - vc1vec[1]);
  const lambertDv2 = Math.hypot(vc2vec[0] - sol.v2[0], vc2vec[1] - sol.v2[1]);
  return {
    pResidual: sol.pKm - pStarKm,
    eResidual: sol.e - e,
    dv1Residual: lambertDv1 - dv1Km,
    dv2Residual: lambertDv2 - dv2Km,
    iterations: sol.iterations,
  };
}

// --- RK4 numerical propagation (mirrors otca_c/rk4.py) -----------------
// The THIRD independent cross-check: takes the transfer orbit's own
// departure state vector (position/velocity, from the closed-form
// conic-to-Cartesian formulas below -- no root-finding), integrates
// Newton's law forward by direct fixed-step RK4 for the already-computed
// time of flight, and compares the result to the already-known analytic
// arrival state. A step-doubling convergence study (Richardson-style)
// accompanies every call, exactly as in the Python source: a single RK4
// run at one step count proves nothing about accuracy by itself, so this
// always runs at both n and 2n steps and verifies the observed order is
// consistent with genuine 4th-order behavior before trusting the result.
// This is a line-by-line port of rk4.py, not an independent
// reimplementation -- see that file's own module docstring for the full
// rationale (why this differs from both geometryCheck and
// lambertCrossCheck above).
const RK4_NOISE_FLOOR_KM = 1e-6;
const RK4_NOISE_FLOOR_KM_S = 1e-9;

function rk4StateVectorFromElements(pKm, e, nu1Rad, thetaRad, muKm3s2) {
  const nu = nu1Rad + thetaRad;
  const rKm = pKm / (1 + e * Math.cos(nu));
  const hKm2S = Math.sqrt(muKm3s2 * pKm);
  const vTangentialKmS = hKm2S / rKm;
  const vRadialKmS = Math.sqrt(muKm3s2 / pKm) * e * Math.sin(nu);
  const c = Math.cos(thetaRad), s = Math.sin(thetaRad);
  const xKm = rKm * c, yKm = rKm * s;
  const vxKmS = vRadialKmS * c - vTangentialKmS * s;
  const vyKmS = vRadialKmS * s + vTangentialKmS * c;
  return { pos: [xKm, yKm], vel: [vxKmS, vyKmS] };
}

function rk4Accel(pos, muKm3s2) {
  const x = pos[0], y = pos[1];
  const r3 = Math.pow(x * x + y * y, 1.5);
  return [-muKm3s2 * x / r3, -muKm3s2 * y / r3];
}

function rk4Step(pos, vel, muKm3s2, hS) {
  function deriv(p, v) { return [v, rk4Accel(p, muKm3s2)]; }
  const d1 = deriv(pos, vel);
  const d2 = deriv(
    [pos[0] + 0.5 * hS * d1[0][0], pos[1] + 0.5 * hS * d1[0][1]],
    [vel[0] + 0.5 * hS * d1[1][0], vel[1] + 0.5 * hS * d1[1][1]]
  );
  const d3 = deriv(
    [pos[0] + 0.5 * hS * d2[0][0], pos[1] + 0.5 * hS * d2[0][1]],
    [vel[0] + 0.5 * hS * d2[1][0], vel[1] + 0.5 * hS * d2[1][1]]
  );
  const d4 = deriv(
    [pos[0] + hS * d3[0][0], pos[1] + hS * d3[0][1]],
    [vel[0] + hS * d3[1][0], vel[1] + hS * d3[1][1]]
  );
  const newPos = [
    pos[0] + (hS / 6) * (d1[0][0] + 2 * d2[0][0] + 2 * d3[0][0] + d4[0][0]),
    pos[1] + (hS / 6) * (d1[0][1] + 2 * d2[0][1] + 2 * d3[0][1] + d4[0][1]),
  ];
  const newVel = [
    vel[0] + (hS / 6) * (d1[1][0] + 2 * d2[1][0] + 2 * d3[1][0] + d4[1][0]),
    vel[1] + (hS / 6) * (d1[1][1] + 2 * d2[1][1] + 2 * d3[1][1] + d4[1][1]),
  ];
  return { pos: newPos, vel: newVel };
}

function rk4Propagate(pos0, vel0, muKm3s2, tofS, nSteps) {
  const hS = tofS / nSteps;
  let pos = pos0, vel = vel0;
  const r0 = Math.hypot(pos0[0], pos0[1]);
  const energy0 = 0.5 * (vel0[0] * vel0[0] + vel0[1] * vel0[1]) - muKm3s2 / r0;
  const angMom0 = pos0[0] * vel0[1] - pos0[1] * vel0[0];
  let maxEnergyDrift = 0, maxAngMomDrift = 0;
  for (let i = 0; i < nSteps; i++) {
    const st = rk4Step(pos, vel, muKm3s2, hS);
    pos = st.pos; vel = st.vel;
    const r = Math.hypot(pos[0], pos[1]);
    if (!(r > 0) || !isFinite(r)) {
      throw new ConvergenceError(
        "RK4 propagation diverged (non-finite or non-positive radius) " +
        "after a step with h_s=" + hS + " s, n_steps=" + nSteps
      );
    }
    const energy = 0.5 * (vel[0] * vel[0] + vel[1] * vel[1]) - muKm3s2 / r;
    const angMom = pos[0] * vel[1] - pos[1] * vel[0];
    if (energy0 !== 0) maxEnergyDrift = Math.max(maxEnergyDrift, Math.abs((energy - energy0) / energy0));
    if (angMom0 !== 0) maxAngMomDrift = Math.max(maxAngMomDrift, Math.abs((angMom - angMom0) / angMom0));
  }
  return { pos: pos, vel: vel, maxEnergyDrift: maxEnergyDrift, maxAngMomDrift: maxAngMomDrift };
}

function rk4Validate(pStarKm, e, nu1Rad, deltaThetaRad, tofS, muKm3s2, r2Km, opts) {
  opts = opts || {};
  let nSteps = opts.nSteps || 4000;
  const positionTolKm = opts.positionTolKm === undefined ? 1e-3 : opts.positionTolKm;
  const velocityTolKmS = opts.velocityTolKm === undefined ? 1e-6 : opts.velocityTolKm;
  const maxDoublings = opts.maxDoublings === undefined ? 8 : opts.maxDoublings;

  if (!(e >= 0 && e < 1)) {
    throw new ConvergenceError("rk4Validate requires an ellipse, 0<=e<1; got e=" + e);
  }
  if (!(tofS > 0)) {
    throw new ConvergenceError("time of flight must be positive, got tofS=" + tofS);
  }

  const dep = rk4StateVectorFromElements(pStarKm, e, nu1Rad, 0, muKm3s2);
  const arr = rk4StateVectorFromElements(pStarKm, e, nu1Rad, deltaThetaRad, muKm3s2);
  function dist(u, v) { return Math.hypot(u[0] - v[0], u[1] - v[1]); }

  let runA = rk4Propagate(dep.pos, dep.vel, muKm3s2, tofS, nSteps);
  let posResidualA = dist(runA.pos, arr.pos);

  for (let attempt = 0; attempt < maxDoublings; attempt++) {
    const finerSteps = 2 * nSteps;
    const runB = rk4Propagate(dep.pos, dep.vel, muKm3s2, tofS, finerSteps);
    const posResidualB = dist(runB.pos, arr.pos);
    const velResidualB = dist(runB.vel, arr.vel);

    let observedOrder = null;
    if (posResidualA > RK4_NOISE_FLOOR_KM && posResidualB > RK4_NOISE_FLOOR_KM) {
      const ratio = posResidualA / posResidualB;
      if (ratio <= 1.0) {
        throw new ConvergenceError(
          "doubling the RK4 step count did not reduce the residual against the analytic arrival state (" +
          posResidualA + " km at n_steps=" + nSteps + " vs " + posResidualB + " km at n_steps=" + finerSteps +
          " -- ratio=" + ratio + "), inconsistent with a correctly implemented convergent integrator."
        );
      }
      observedOrder = Math.log2(ratio);
      if (observedOrder < 2.5) {
        throw new ConvergenceError(
          "RK4 step-doubling did not show the expected ~4th-order convergence (observed order=" +
          observedOrder.toFixed(2) + ")."
        );
      }
    }

    if ((posResidualB < positionTolKm && velResidualB < velocityTolKmS) ||
        (posResidualB <= RK4_NOISE_FLOOR_KM && velResidualB <= RK4_NOISE_FLOOR_KM_S)) {
      return {
        positionResidualKm: posResidualB,
        velocityResidualKmS: velResidualB,
        nSteps: finerSteps,
        observedConvergenceOrder: observedOrder,
        maxEnergyDrift: Math.max(runA.maxEnergyDrift, runB.maxEnergyDrift),
        maxAngMomDrift: Math.max(runA.maxAngMomDrift, runB.maxAngMomDrift),
      };
    }

    nSteps = finerSteps;
    runA = runB;
    posResidualA = posResidualB;
  }

  throw new ConvergenceError(
    "RK4 validation did not reach tolerance within " + maxDoublings +
    " step-count doublings (reached n_steps=" + nSteps + ", position_residual_km=" + posResidualA + ")."
  );
}

return {
  ConvergenceError: ConvergenceError,
  brentq: brentq,
  costS1S2: costS1S2,
  cost: cost,
  costDerivative: costDerivative,
  findPStar: findPStar,
  solveHohmann: solveHohmann,
  geometryCheck: geometryCheck,
  circularSpeedKmS: circularSpeedKmS,
  ECCENTRICITY_CONDITIONING_THRESHOLD: ECCENTRICITY_CONDITIONING_THRESHOLD,
  eccentricityConditioningWarning: eccentricityConditioningWarning,
  periapsisRadiusKm: periapsisRadiusKm,
  apoapsisRadiusKm: apoapsisRadiusKm,
  orbitalPeriodS: orbitalPeriodS,
  trueToEccentric: trueToEccentric,
  eccentricToTrue: eccentricToTrue,
  meanToEccentric: meanToEccentric,
  timeOfFlightS: timeOfFlightS,
  propellantMassFraction: propellantMassFraction,
  stumpffC: stumpffC,
  stumpffS: stumpffS,
  solveLambert: solveLambert,
  lambertCrossCheck: lambertCrossCheck,
  rk4Validate: rk4Validate,
};

});
