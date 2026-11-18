"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Licensed under the Research, Education and Non-Commercial Use License
(NCSL-1.0). See LICENSE-CODE.md at the repository root for the full terms.

Tests for `curvature_regime.py`, including a pin of its certified
constants (`W_C_CERTIFIED`, `PHI_C_CERTIFIED_DEG`) against the
repository's own certificate/script files, so that a future
re-certification with a different value is caught here rather than
silently diverging from this module's copy.
"""
import math
from pathlib import Path

import pytest

from otca_c.curvature_regime import (
    PHI_C_CERTIFIED_DEG,
    W_C_CERTIFIED,
    classify_curvature_regime,
)
from otca_c.errors import DomainError

# Repository root, if this package is checked out inside the full
# repository (sibling to `circular/`); None if not (e.g. a standalone
# package install) -- the pin tests below skip gracefully in that case
# rather than failing on an environment that never had the source
# certificate available to begin with.
_REPO_ROOT_CANDIDATES = [
    Path(__file__).resolve().parents[1],  # .../2i-OTCA-C-v1.0-algorithm
    Path(__file__).resolve().parents[2] if len(Path(__file__).resolve().parents) > 2 else None,
]


def _find_certificate_script() -> "Path | None":
    for root in _REPO_ROOT_CANDIDATES:
        if root is None:
            continue
        candidate = (
            root / "circular" / "theorem_4" / "03_central_mesh"
            / "mesh_certificate_R2.py"
        )
        if candidate.is_file():
            return candidate
        candidate2 = root.parent / "circular" / "theorem_4" / "03_central_mesh" / "mesh_certificate_R2.py"
        if candidate2.is_file():
            return candidate2
    return None


def test_wc_and_phic_pinned_against_certificate_script():
    script = _find_certificate_script()
    if script is None:
        pytest.skip("certificate script not found in this checkout; skipping pin test")
    text = script.read_text(encoding="utf-8")
    assert "0.30698081390187630871" in text, (
        "W_C_CERTIFIED no longer matches the literal in "
        "mesh_certificate_R2.py -- update curvature_regime.py's "
        "W_C_CERTIFIED deliberately if this is an intentional "
        "re-certification, do not silently let this drift."
    )
    assert "164.3710517" in text, (
        "PHI_C_CERTIFIED_DEG no longer matches the literal in "
        "mesh_certificate_R2.py -- same caution as W_C_CERTIFIED above."
    )


def test_r_not_2_raises_domain_error():
    with pytest.raises(DomainError):
        classify_curvature_regime(6678.0, 6678.0 * 2.5, 170.0)
    with pytest.raises(DomainError):
        classify_curvature_regime(6678.0, 6678.0 * 1.9999, 170.0)


def test_at_inflection_point():
    r1_km = 6678.0
    r2_km = r1_km * 2.0
    delta_theta_at = 180.0 - PHI_C_CERTIFIED_DEG
    result = classify_curvature_regime(r1_km, r2_km, delta_theta_at)
    assert result.regime == "at_inflection"
    assert math.isclose(result.w, W_C_CERTIFIED, rel_tol=1e-6)


def test_pre_and_post_inflection_sides():
    r1_km = 6678.0
    r2_km = r1_km * 2.0

    near_hohmann = classify_curvature_regime(r1_km, r2_km, 179.0)
    assert near_hohmann.regime == "post_inflection"
    assert near_hohmann.w > W_C_CERTIFIED

    near_far_boundary = classify_curvature_regime(r1_km, r2_km, 2.0)
    assert near_far_boundary.regime == "pre_inflection"
    assert near_far_boundary.w < W_C_CERTIFIED


def test_exact_hohmann_does_not_raise():
    r1_km = 6678.0
    r2_km = r1_km * 2.0
    result = classify_curvature_regime(r1_km, r2_km, 180.0)
    assert result.regime == "post_inflection"  # Hohmann is deep in H>0 territory


def test_invalid_radii_raise_domain_error():
    with pytest.raises(DomainError):
        classify_curvature_regime(-1.0, 6678.0 * 2.0, 170.0)
