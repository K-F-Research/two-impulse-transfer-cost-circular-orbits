"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Regression test against the archived independent-comparison grid
(`Independent_comparison/grid_compare_full_results.json`, Run B, 230
cases: R in {1.02 .. 20}, angular deviation from Hohmann in {1 .. 30} deg,
r1=6678 km, mu=Earth). That file is the frozen record of a cross-check
already performed once (an independently written, differently formulated
p-iteration Lambert solver run at 50-digit precision vs. this package's
`api.compute_transfer`); see `Independent_comparison/FINAL_REPORT.md` for
the full methodology.

This test does NOT re-run the independent solver (that would duplicate
`Independent_comparison/grid_compare_full.py`, a non-pytest script with
its own mpmath dependency chain). Instead it re-runs ONLY this package's
own `api.compute_transfer` at the same 230 (R, dev) pairs and checks the
result against the `dv_core` value recorded in that file at the time the
cross-check was performed. This is a regression guard, distinct in
purpose from the cross-check itself: it catches any future unintended
change to `core.py`/`api.py` that would silently alter results, by
comparing against a frozen snapshot rather than against an independent
method. A failure here means "the exact solver's output has changed
since this snapshot was taken", not "the independent method disagrees".

The reference data file is copied verbatim into `tests/` (not read
from `Independent_comparison/` by relative path) so this test suite has
no dependency on the repository layout outside the installable package's
own `tests/` directory.
"""
import json
import math
from pathlib import Path

import pytest

from otca_c import api

_DATA_PATH = Path(__file__).parent / "grid_compare_full_results.json"
_MU_EARTH = 398600.4418
_R1_KM = 6678.0

with open(_DATA_PATH, "r", encoding="utf-8") as _f:
    _GRID_CASES = json.load(_f)

assert len(_GRID_CASES) == 230, (
    f"expected the archived Run B grid to contain exactly 230 cases "
    f"(23 R values x 10 deviation values, per Independent_comparison/"
    f"FINAL_REPORT.md), found {len(_GRID_CASES)}; the reference file may "
    f"have been replaced with a different run -- update this constant "
    f"deliberately if so, do not silently relax it."
)


@pytest.mark.parametrize(
    "case",
    _GRID_CASES,
    ids=[f"R={c['R']}_dev={c['dev']}deg" for c in _GRID_CASES],
)
def test_exact_solver_matches_archived_snapshot(case: dict) -> None:
    r2_km = case["R"] * _R1_KM
    delta_theta_deg = 180.0 - case["dev"]
    result = api.compute_transfer(
        _R1_KM, r2_km, delta_theta_deg, mu_km3_s2=_MU_EARTH, validate=False
    )
    assert math.isfinite(result.delta_v_total_km_s)
    # Tight tolerance: this is a same-machine, same-algorithm re-run, not
    # a cross-method comparison -- any drift beyond floating-point noise
    # indicates an unintended change in core.py/api.py since the snapshot
    # was archived.
    assert result.delta_v_total_km_s == pytest.approx(case["dv_core"], rel=1e-9)


def test_archived_snapshot_itself_matches_its_own_claimed_tolerance() -> None:
    """Sanity check on the fixture file's own internal consistency: every
    row's claimed `rel_diff` (independent method vs. core, at archival
    time) must itself be below the 5.7e-12 worst case stated in
    `Independent_comparison/FINAL_REPORT.md`. This does not re-verify the
    independent method (out of scope here), only that the archived numbers
    are not corrupted or hand-edited since FINAL_REPORT.md was written."""
    worst = max(case["rel_diff"] for case in _GRID_CASES)
    assert worst < 6e-12, (
        f"archived worst-case rel_diff={worst:.3e} exceeds the 5.7e-12 "
        "figure documented in Independent_comparison/FINAL_REPORT.md; "
        "the fixture file may have been edited or regenerated from a "
        "different run without updating that report."
    )
