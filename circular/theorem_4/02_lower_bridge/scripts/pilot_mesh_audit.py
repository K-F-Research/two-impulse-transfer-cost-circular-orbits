# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
"""
pilot_mesh_audit.py

A second, INDEPENDENT program: it recomputes no symbolic or interval
proof at all -- it only re-reads the JSON archive produced by
full_mesh_bridge.py and checks, from the bounds already recorded there,
that every criterion is satisfied and that coverage and stitching are
correct. It fails EXPLICITLY (message and non-zero exit code) on the
first defect found.

Checks performed:
  1. Exact coverage of the stated domain, with no gap and no overlap.
  2. For each tile: G_Y<0, S1>0, S2>0, the denominator of G separated
     from 0 by at least DENG_SEP_MIN, G_lo>0, G_hi<0, width(Y_B) <=
     YB_WIDTH_MAX, margin -sup(Q) > MARGIN_TOL (recomputed from the
     archived bounds, not from a new symbolic evaluation).
  3. Every archived tile-to-tile stitching record is "OK" and the
     records cover every internal boundary (n-1 stitchings for n tiles).
"""
import json
import mpmath as mp
import sys

mp.mp.dps = 50


def mpf_of(x):
    """x is either a string (exact value) or a [lo, hi] list."""
    if isinstance(x, list):
        return mp.mpf(x[0]), mp.mpf(x[1])
    v = mp.mpf(x)
    return v, v


def fail(msg):
    print(f"AUDIT FAIL: {msg}")
    sys.exit(1)


def main(path):
    with open(path) as f:
        archive = json.load(f)

    print(f"Script version archive: {archive['script_version']}")
    print(f"mp.dps={archive['mp_dps']}  iv.dps={archive['iv_dps']}")
    print(f"Parameters: {archive['params']}")

    MARGIN_TOL = mp.mpf(archive['params']['MARGIN_TOL'])
    YB_WIDTH_MAX = mp.mpf(archive['params']['YB_WIDTH_MAX'])
    DENG_SEP_MIN = mp.mpf(archive['params']['DENG_SEP_MIN'])

    domain_lo, domain_hi = (mp.mpf(v) for v in archive['domain'])
    tiles = archive['tiles']
    n = len(tiles)
    print(f"\nStated domain: [{domain_lo}, {domain_hi}]")
    print(f"Number of tiles: {n}")

    if archive['n_failures'] != 0:
        fail(f"{archive['n_failures']} unrecoverable failure(s) recorded in the archive")

    if n == 0:
        fail("no tile accepted -- nothing to audit")

    # --- 1. Exact coverage, no gap, no overlap ---
    tiles_sorted = sorted(tiles, key=lambda r: mp.mpf(r['w_lo']))
    if mp.mpf(tiles_sorted[0]['w_lo']) != domain_lo:
        fail(f"the first tile does not start at the domain's left edge: "
             f"{tiles_sorted[0]['w_lo']} != {domain_lo}")
    if mp.mpf(tiles_sorted[-1]['w_hi']) != domain_hi:
        fail(f"the last tile does not end at the domain's right edge: "
             f"{tiles_sorted[-1]['w_hi']} != {domain_hi}")
    for i in range(1, n):
        prev_hi = mp.mpf(tiles_sorted[i-1]['w_hi'])
        cur_lo = mp.mpf(tiles_sorted[i]['w_lo'])
        if prev_hi != cur_lo:
            fail(f"gap or overlap between tile {i-1} and {i}: "
                 f"w_hi={prev_hi} != w_lo={cur_lo}")
    print("Coverage: PASS (contiguous, no gap, no overlap)")

    # --- 2. Per-tile criteria, recomputed from the archived bounds ---
    violations = []
    margins = []
    for idx, r in enumerate(tiles_sorted):
        GY_lo, GY_hi = mpf_of(r['GY'])
        S1_lo, S1_hi = mpf_of(r['S1'])
        S2_lo, S2_hi = mpf_of(r['S2'])
        denG_lo, denG_hi = mpf_of(r['denG'])
        Glo_lo, Glo_hi = mpf_of(r['G_lo'])
        Ghi_lo, Ghi_hi = mpf_of(r['G_hi'])
        YB_lo, YB_hi = mpf_of(r['YB'])
        Q_lo, Q_hi = mpf_of(r['Q'])

        if not (GY_hi < 0):
            violations.append((idx, "GY_not_negative"))
        if not (S1_lo > 0):
            violations.append((idx, "S1_not_positive"))
        if not (S2_lo > 0):
            violations.append((idx, "S2_not_positive"))
        den_sep_ok = (denG_lo > 0 or denG_hi < 0) and min(abs(denG_lo), abs(denG_hi)) >= DENG_SEP_MIN
        if not den_sep_ok:
            violations.append((idx, "denG_not_separated"))
        if not (Glo_lo > 0):
            violations.append((idx, "G_lo_not_positive"))
        if not (Ghi_hi < 0):
            violations.append((idx, "G_hi_not_negative"))
        if (YB_hi - YB_lo) > YB_WIDTH_MAX:
            violations.append((idx, "YB_too_wide"))
        margin_recomputed = -Q_hi
        if not (margin_recomputed > MARGIN_TOL):
            violations.append((idx, f"margin_insufficient ({margin_recomputed})"))
        margins.append(margin_recomputed)

    if violations:
        for idx, reason in violations:
            print(f"  TILE {idx}: {reason}")
        fail(f"{len(violations)} criterion violation(s) detected")
    print(f"Per-tile criteria (G_Y, S1, S2, denG, G_lo, G_hi, width of YB, margin): "
          f"PASS for all {n} tiles")
    print(f"Minimum margin (recomputed from the archive): {min(margins)}")
    print(f"Maximum margin: {max(margins)}")

    # --- 3. Tile-to-tile stitching ---
    stitching = archive['stitching']
    if len(stitching) != n - 1:
        fail(f"number of archived stitchings ({len(stitching)}) != n_tiles-1 ({n-1})")
    for i, st in enumerate(stitching):
        if st['status'] != 'OK':
            fail(f"stitching {i} (between tiles {st['between']}): status = {st['status']}")
        if st['between'] != [i, i+1]:
            fail(f"stitching {i} does not link consecutive tiles: {st['between']}")
    print(f"Tile-to-tile stitching: PASS for all {len(stitching)} internal boundaries")

    print(f"\n=== GLOBAL AUDIT: PASS ===")
    print(f"Domain [{domain_lo},{domain_hi}] certified by {n} tiles, "
          f"{len(stitching)} stitchings verified, 0 violation.")
    return 0


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else '/tmp/pilot_mesh_archive.json'
    sys.exit(main(path))
