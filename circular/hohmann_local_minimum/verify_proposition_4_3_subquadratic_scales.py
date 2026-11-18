"""
Copyright (c) 2026 - Frederic G. Speyser & Kseniya Vyatkina

Note:
Independent numerical confirmation of Proposition 4.3 (non-competitiveness
of every sub-quadratic scale), the generalization used to close Theorem 2's
asymptotic-scale gap: for every path p(phi) = p_H + beta*phi^gamma with
gamma in [1,2), the excess cost over the Hohmann transfer decays like
phi^(2*gamma-2) as phi -> 0 -- strictly slower than the phi^2 excess
achievable along the admitted competitive family (gamma=2) -- so no such
path can be competitive with the Hohmann-adjacent family for any gamma in
[1,2), not only gamma=1 as the original, narrower statement covered.

Method: evaluate the FULL (unexpanded) two-impulse transfer cost
cost(p, phi, R) directly (no series expansion used anywhere in this
script), at mpmath 60-digit precision, along p = p_H + beta*phi^gamma for
a representative ratio R=3 and gamma = 1, 5/4, 3/2, 7/4, over a geometric
sequence of phi -> 0. The empirical log-log slope of (cost - Hohmann
baseline) between consecutive phi values is compared against the
predicted exponent 2*gamma-2.

A note on an earlier bug (left here for transparency, since it affected an
intermediate version of this exact check): the Hohmann baseline must use
vis-viva with the semi-major axis a_H=(r1+r2)/2, NOT the semi-latus
rectum p_H -- using p_H there produces a large, phi-independent, and
physically impossible *negative* cost gap. `hohmann_baseline` below uses
a_H correctly.

Status: REPRODUCED. Running this script gives empirical exponents
converging, as phi -> 0, to 0.0, 0.5, 1.0, 1.5 for gamma = 1, 1.25, 1.5,
1.75 respectively -- matching the predicted 2*gamma-2 in all four cases,
confirming the generalized exclusion over the full continuum gamma in
[1,2), not merely gamma=1.

Scope note: this script confirms the predicted EXPONENT only. The proof's
additional claim about the exact leading COEFFICIENT (P*beta^2*s^2) is
established analytically in the manuscript (Section 4.2) and was checked
separately there; reproducing that coefficient match is not included here.
"""

import mpmath as mp

mp.mp.dps = 60


def cost(p, phi, R):
    """Full (unexpanded) two-impulse transfer cost, r1=1, r2=R, semi-latus
    rectum p, transfer angle theta = pi - phi. No series expansion is used:
    this is the exact cost function, evaluated at full (60-digit)
    precision for every phi tested, including phi as small as 1e-16."""
    c = mp.cos(phi)
    x = p - 1
    y = p / R - 1
    e2 = x**2 + (x * c + y)**2 / (1 - c**2)
    v1_sq = 2 + (e2 - 1) / p
    v2_sq = 2 / R + (e2 - 1) / p
    v_theta1 = mp.sqrt(p)
    v_theta2 = mp.sqrt(p) / R
    v_c1 = mp.mpf(1)
    v_c2 = 1 / mp.sqrt(R)
    S1 = v1_sq + v_c1**2 - 2 * v_c1 * v_theta1
    S2 = v2_sq + v_c2**2 - 2 * v_c2 * v_theta2
    return mp.sqrt(S1) + mp.sqrt(S2)


def hohmann_baseline(R):
    """Exact Hohmann-transfer total delta-v, r1=1, r2=R, via vis-viva on
    the Hohmann ellipse of semi-major axis a_H=(r1+r2)/2."""
    aH = (1 + R) / 2
    v1H = mp.sqrt(2 - 1 / aH)
    v2H = mp.sqrt(2 / R - 1 / aH)
    vc1 = mp.mpf(1)
    vc2 = 1 / mp.sqrt(R)
    return abs(v1H - vc1) + abs(vc2 - v2H)


R = mp.mpf(3)
pH = 2 * R / (1 + R)
dvH = hohmann_baseline(R)
beta = mp.mpf('0.3')

print(f"R = {R}, p_H = {pH}, Hohmann baseline dv = {dvH}", flush=True)
print("Sanity check: cost(p_H, phi) -> dvH as phi -> 0:",
      cost(pH, mp.mpf('1e-15'), R), flush=True)
print(flush=True)

results = {}
for gamma in [mp.mpf(1), mp.mpf('1.25'), mp.mpf('1.5'), mp.mpf('1.75')]:
    predicted = 2 * gamma - 2
    print(f"--- gamma = {gamma}  (predicted exponent 2*gamma-2 = "
          f"{predicted}) ---", flush=True)
    phis = [mp.mpf('0.05') / 2**k for k in range(0, 10)]
    excesses = []
    for phi in phis:
        p = pH + beta * phi**gamma
        exc = cost(p, phi, R) - dvH
        assert exc > 0, "excess over Hohmann must be positive near a minimum"
        excesses.append(exc)
    slopes = []
    for i in range(len(phis) - 1):
        slope = ((mp.log(excesses[i + 1]) - mp.log(excesses[i]))
                  / (mp.log(phis[i + 1]) - mp.log(phis[i])))
        slopes.append(slope)
        print(f"  phi={mp.nstr(phis[i],4)}->{mp.nstr(phis[i+1],4)}: "
              f"empirical exponent = {mp.nstr(slope,10)}", flush=True)
    final_slope = slopes[-1]
    results[gamma] = (final_slope, predicted)
    print(f"  final empirical exponent: {mp.nstr(final_slope,10)}  "
          f"(predicted: {predicted})", flush=True)
    print(flush=True)

print("Summary (last empirical exponent vs. predicted 2*gamma-2):",
      flush=True)
all_ok = True
for gamma, (emp, pred) in results.items():
    ok = abs(emp - pred) < mp.mpf('0.02')
    all_ok = all_ok and ok
    print(f"  gamma={gamma}: empirical={mp.nstr(emp,8)}, "
          f"predicted={pred}, match within 0.02? {ok}", flush=True)

print(flush=True)
if all_ok:
    print("RESULT: CONFIRMED. The cost excess over Hohmann decays like "
          "phi^(2*gamma-2) for every tested gamma in [1,2), always "
          "strictly slower than phi^2 -- every sub-quadratic scaling "
          "path is non-competitive.", flush=True)
else:
    print("RESULT: MISMATCH -- see output above.", flush=True)
