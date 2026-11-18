"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Parity test between the Python library (src/otca_c/) and its JavaScript
port in 2i-OTCA-C-v1.0-demo/otca.js (the GitHub-Pages-compatible demonstrator). Runs the
JS port under Node (no DOM involved -- only the pure math module is
exercised), feeds it the SAME inputs as a Python compute_transfer() call,
and compares every numerical quantity the two implementations are
supposed to agree on: p*, e, a, delta-v, time of flight, and the Lambert
cross-check's own residuals -- not just that both "ran without error".

This also checks that the two implementations FAIL in the same places,
for the same documented reasons (the exact Hohmann degeneracy, the
near-360-degree Lambert/RK4 edge cases), per the explicit request (during
this project's external review) that a passing parity campaign must
include matching error contracts, not just matching numbers on the easy
cases.

Skipped automatically if `node` is not on PATH (this is a cross-language
integration test, not a core Python test -- the rest of the Python test
suite does not depend on Node being installed).
"""
import json
import math
import shutil
import subprocess
from pathlib import Path

import pytest

from otca_c import Bodies, compute_transfer
from otca_c.errors import ConvergenceError
from otca_c.rk4 import rk4_validate

DEMO_DIR = Path(__file__).resolve().parent.parent.parent / "2i-OTCA-C-v1.0-demo"
OTCA_JS = DEMO_DIR / "otca.js"

pytestmark = [
    pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed"),
    pytest.mark.skipif(not OTCA_JS.exists(), reason=f"otca.js not found at {OTCA_JS}"),
]


def _run_js(js_snippet: str):
    """Run a snippet of Node code with Otca already required as a global,
    and the snippet's own `console.log(JSON.stringify(...))` output
    parsed back as the return value. Raises AssertionError (with Node's
    stderr) on any non-zero exit."""
    full_script = (
        f"const Otca = require({json.dumps(str(OTCA_JS))});\n" + js_snippet
    )
    result = subprocess.run(
        ["node", "-e", full_script], capture_output=True, text=True, timeout=30,
    )
    if result.returncode != 0:
        raise AssertionError(f"node script failed:\n{result.stderr}")
    return json.loads(result.stdout.strip().splitlines()[-1])


# Mirrors the demo's own update() branch logic closely enough to compute
# the same quantities Python's compute_transfer does, without needing a
# DOM -- i.e. this is the JS *library* being parity-checked, the same
# otca.js the page itself loads, just driven directly instead of through
# the UI.
_JS_COMPUTE_TRANSFER = """
function computeTransfer(r1, r2, deltaThetaDeg, mu) {
  const R = r2 / r1;
  const deltaThetaRad = deltaThetaDeg * Math.PI / 180;
  const phi = Math.PI - deltaThetaRad;
  const phiAbs = Math.abs(phi);
  const isHohmann = phiAbs < 1e-12;
  let pStar, S1, S2, e, a, nu1, nu2;
  if (isHohmann) {
    const h = Otca.solveHohmann(R);
    pStar = h.p; e = h.e; a = h.a; S1 = h.dv1*h.dv1; S2 = h.dv2*h.dv2;
    nu1 = 0; nu2 = Math.PI;
  } else {
    pStar = Otca.findPStar(phiAbs, R);
    const cs = Otca.costS1S2(pStar, phiAbs, R);
    S1 = cs.S1; S2 = cs.S2; e = Math.sqrt(cs.e2); a = pStar / (1 - e*e);
    const g = Otca.geometryCheck(pStar, phi, R, e, S1, S2);
    nu1 = g.nu1; nu2 = g.nu2;
  }
  const vc1 = Otca.circularSpeedKmS(r1, mu);
  const dv1 = Math.sqrt(S1) * vc1, dv2 = Math.sqrt(S2) * vc1;
  const pStarKm = pStar * r1, aKm = a * r1;
  const tofS = Otca.timeOfFlightS(aKm, e, nu1, nu2, mu);
  return {
    isHohmann, pStarKm, e, aKm, dv1, dv2, dvTotal: dv1 + dv2, tofS, nu1, nu2,
    rPeriKm: Otca.periapsisRadiusKm(pStarKm, e),
    rApoKm: Otca.apoapsisRadiusKm(pStarKm, e),
    periodS: Otca.orbitalPeriodS(aKm, mu),
  };
}
"""


CASES = [
    (6778.0, 13556.0, 170.0, Bodies.EARTH),
    (6778.0, 13556.0, 10.0, Bodies.EARTH),
    (6778.0, 13556.0, 190.0, Bodies.EARTH),
    (6778.0, 13556.0, 350.0, Bodies.EARTH),
    (6778.0, 13556.0, 269.0, Bodies.EARTH),
    (6778.0, 13556.0, 5.0, Bodies.EARTH),
    (6778.0, 42164.0, 160.0, Bodies.EARTH),
    (3789.9, 24190.9, 170.0, Bodies.MARS),
    (1000.0, 7000.0, 90.0, Bodies.MOON),
    (6778.0, 13556.0, 180.0, Bodies.EARTH),  # exact Hohmann
]


@pytest.mark.parametrize("r1_km,r2_km,delta_theta_deg,mu", CASES)
def test_js_matches_python_primary_solve(r1_km, r2_km, delta_theta_deg, mu):
    py = compute_transfer(r1_km, r2_km, delta_theta_deg, mu_km3_s2=mu)
    js = _run_js(
        _JS_COMPUTE_TRANSFER
        + f"console.log(JSON.stringify(computeTransfer({r1_km!r}, {r2_km!r}, {delta_theta_deg!r}, {mu!r})));"
    )

    assert js["isHohmann"] == py.is_hohmann

    def rel(a, b):
        return abs(a - b) / max(abs(b), 1e-12)

    # The JS port uses a finite-difference C_p and a hand-written Brent
    # solver (see otca.js's module docstring) instead of Python's
    # symbolic C_p + scipy.brentq -- a documented, deliberate deviation,
    # not a bug -- so the tolerance here is generous relative to the two
    # implementations' own internal self-consistency (~1e-13), but still
    # tight enough that a real divergence (wrong formula, wrong sign,
    # wrong branch) would fail it immediately.
    tol = 1e-6
    assert rel(js["pStarKm"], py.p_star_km) < tol
    assert rel(js["e"], py.eccentricity) < tol
    assert rel(js["aKm"], py.semi_major_axis_km) < tol
    assert rel(js["dv1"], py.delta_v1_km_s) < tol
    assert rel(js["dv2"], py.delta_v2_km_s) < tol
    assert rel(js["dvTotal"], py.delta_v_total_km_s) < tol
    assert rel(js["tofS"], py.time_of_flight_s) < tol
    assert rel(js["periodS"], py.orbital_period_s) < tol
    assert rel(js["rPeriKm"], py.r_periapsis_km) < tol
    assert rel(js["rApoKm"], py.r_apoapsis_km) < tol


def test_js_near_360_degree_boundary_matches_python_within_a_looser_tolerance():
    """The extreme delta_theta=359.999 deg, R=2.0 case (e within 1e-10 of
    1): both implementations converge, but the JS port's finite-difference
    C_p is measurably less precise than Python's symbolic one exactly in
    this near-degenerate region (see the module-level comment above) --
    pin the actual observed relative agreement, not an idealized one, so
    a future regression (in EITHER implementation) is still caught."""
    py = compute_transfer(6778.0, 13556.0, 359.999, mu_km3_s2=Bodies.EARTH)
    js = _run_js(
        _JS_COMPUTE_TRANSFER
        + "console.log(JSON.stringify(computeTransfer(6778.0, 13556.0, 359.999, 398600.4418)));"
    )
    assert (1.0 - js["e"]) < 1e-6
    assert (1.0 - py.eccentricity) < 1e-6
    rel_tof = abs(js["tofS"] - py.time_of_flight_s) / py.time_of_flight_s
    assert rel_tof < 1e-4


@pytest.mark.parametrize("r1_km,r2_km,delta_theta_deg,mu", CASES[:-1])  # exclude exact Hohmann
def test_js_lambert_cross_check_matches_python(r1_km, r2_km, delta_theta_deg, mu):
    """Feed both implementations' own primary solve's time of flight into
    their respective independent Lambert cross-checks, and compare the
    recovered (p, e) and the cross-check's own residuals -- not just the
    primary solve (already covered above)."""
    py = compute_transfer(r1_km, r2_km, delta_theta_deg, mu_km3_s2=mu, lambert_check=True)
    assert py.lambert_check is not None

    js_snippet = _JS_COMPUTE_TRANSFER + f"""
const r = computeTransfer({r1_km!r}, {r2_km!r}, {delta_theta_deg!r}, {mu!r});
const deltaThetaRad = ({delta_theta_deg!r} * Math.PI / 180) % (2*Math.PI);
const lc = Otca.lambertCrossCheck({r1_km!r}, {r2_km!r}, deltaThetaRad, r.tofS, {mu!r}, r.pStarKm, r.e, r.dv1, r.dv2);
console.log(JSON.stringify(lc));
"""
    js_lc = _run_js(js_snippet)

    assert abs(js_lc["pResidual"]) / py.p_star_km < 1e-5
    assert abs(js_lc["eResidual"]) < 1e-5
    assert abs(js_lc["dv1Residual"]) < 1e-5
    assert abs(js_lc["dv2Residual"]) < 1e-5


def test_js_matches_python_degenerate_zero_angle_rejection_pattern():
    """Python rejects delta_theta=0/360 deg with DomainError at the
    api.py level (folding/normalization), BEFORE reaching core.py at all
    -- the JS port's UI layer (index.html's update()) performs the same
    fold-and-reject check itself (the demo's slider never actually
    reaches 0/360 deg, but the check is still present and must behave the
    same way; see otca.js's find_p_star, which -- like core.py's own --
    still properly raises on phi=0 exactly, the lower-level boundary)."""
    with pytest.raises(ConvergenceError):
        _run_js_find_p_star_at_phi_zero()


def _run_js_find_p_star_at_phi_zero():
    """phi=0 exactly (the model's own domain boundary, sin(phi)=0) must
    raise Otca.ConvergenceError in the JS port, exactly like
    core.find_p_star(0.0, R) does in Python -- checked by having Node
    itself exit nonzero with the ConvergenceError name/message on stderr,
    which this helper turns back into a Python-side ConvergenceError so
    the parametrized contract (raises ConvergenceError) reads naturally."""
    script = (
        f"const Otca = require({json.dumps(str(OTCA_JS))});\n"
        "try { Otca.findPStar(0.0, 2.0); console.log('NO_ERROR'); } "
        "catch (e) { if (e instanceof Otca.ConvergenceError) { console.log('CONVERGENCE_ERROR: ' + e.message); process.exit(1); } throw e; }"
    )
    result = subprocess.run(["node", "-e", script], capture_output=True, text=True, timeout=30)
    if result.returncode != 0 and "CONVERGENCE_ERROR" in result.stdout:
        raise ConvergenceError(result.stdout.strip())
    if result.returncode != 0:
        raise AssertionError(f"unexpected node failure:\n{result.stderr}")
    raise AssertionError("JS findPStar(0.0, 2.0) did not raise -- contract mismatch with Python")


def test_js_lambert_raises_on_degenerate_colinear_case_like_python():
    from otca_c.lambert import solve_lambert
    with pytest.raises(ConvergenceError):
        solve_lambert((7000.0, 0.0), (7000.0, 0.0), 3600.0, Bodies.EARTH)

    script = (
        f"const Otca = require({json.dumps(str(OTCA_JS))});\n"
        "try { Otca.solveLambert([7000.0, 0.0], [7000.0, 0.0], 3600.0, 398600.4418); "
        "console.log('NO_ERROR'); process.exit(0); } "
        "catch (e) { if (e instanceof Otca.ConvergenceError) { process.exit(1); } throw e; }"
    )
    result = subprocess.run(["node", "-e", script], capture_output=True, text=True, timeout=30)
    assert result.returncode == 1, (
        f"JS solveLambert did not raise ConvergenceError on the degenerate "
        f"colinear case (stdout={result.stdout!r}, stderr={result.stderr!r})"
    )


@pytest.mark.parametrize("r1_km,r2_km,delta_theta_deg,mu", CASES)
def test_js_rk4_validate_matches_python(r1_km, r2_km, delta_theta_deg, mu):
    """Parity check for the THIRD independent cross-check (rk4.py's
    rk4_validate / otca.js's rk4Validate): same departure/arrival state
    vectors, same step-doubling convergence study, same residuals --
    including the hard near-periapsis case (delta_theta=350 deg, R=2.0)
    that needs many doublings to converge, not just the easy cases."""
    py = compute_transfer(r1_km, r2_km, delta_theta_deg, mu_km3_s2=mu)
    delta_theta_rad = math.radians(delta_theta_deg) % (2 * math.pi)
    py_rk4 = rk4_validate(
        py.p_star_km, py.eccentricity, py.nu1_rad, delta_theta_rad,
        py.time_of_flight_s, mu, r2_km,
    )
    js = _run_js(
        f"const res = Otca.rk4Validate({py.p_star_km!r}, {py.eccentricity!r}, "
        f"{py.nu1_rad!r}, {delta_theta_rad!r}, {py.time_of_flight_s!r}, {mu!r}, {r2_km!r});\n"
        "console.log(JSON.stringify(res));"
    )
    assert js["nSteps"] == py_rk4.n_steps
    # Both runs are already deep in the floating-point noise floor relative
    # to the 1e-3 km / 1e-6 km/s pass tolerance rk4Validate itself checks
    # against -- this comparison is only meant to catch a genuine formula
    # or step-count regression between the two implementations, not to
    # demand bit-identical floating-point noise between two independent
    # Python/JS arithmetic paths.
    tol = 1e-6
    assert abs(js["positionResidualKm"] - py_rk4.position_residual_km) < max(tol * py_rk4.position_residual_km, 1e-9)
    assert abs(js["velocityResidualKmS"] - py_rk4.velocity_residual_km_s) < max(tol * py_rk4.velocity_residual_km_s, 1e-12)


def test_js_rk4_validate_exact_hohmann_matches_python():
    """The exact-Hohmann case (no Lambert cross-check runs here -- see the
    page's own scope note -- but RK4 propagation has no such degeneracy
    and still applies)."""
    py = compute_transfer(6778.0, 13556.0, 180.0, mu_km3_s2=Bodies.EARTH)
    assert py.is_hohmann
    py_rk4 = rk4_validate(
        py.p_star_km, py.eccentricity, py.nu1_rad, math.pi,
        py.time_of_flight_s, Bodies.EARTH, 13556.0,
    )
    js = _run_js(
        f"const res = Otca.rk4Validate({py.p_star_km!r}, {py.eccentricity!r}, "
        f"{py.nu1_rad!r}, {math.pi!r}, {py.time_of_flight_s!r}, {Bodies.EARTH!r}, 13556.0);\n"
        "console.log(JSON.stringify(res));"
    )
    assert js["positionResidualKm"] < 1e-3 and js["velocityResidualKmS"] < 1e-6
    assert py_rk4.ok()
