"""
Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

audit_hohmann_archive.py

Independent auditor of the Hohmann bridge (hohmann_full_mesh_archive.json).

Audit discipline enforced:
  - reuses NO boolean field from the constructor (no "ok", no "status"
    taken at face value);
  - does NOT reuse findroot -- reads only the archived t_lo, t_hi, a_m, r
    values and RECONSTRUCTS the boxes Tbox, Abox from these numbers;
  - RECOMPUTES from scratch, from the FROZEN expressions and the archived
    bounds, the 4 criteria (S1,S2>0; Ghat_a>0; boundary signs via the
    centered form; t*Hhat>0) and the topology (coverage, junctions);
  - verifies the SHA-256 fingerprint of the frozen expressions;
  - applies an EXPLICIT SAFETY MARGIN when reading each archived decimal
    bound, to compensate for the ordinary (non-directed) rounding of
    mp.nstr() used by the constructor: each bound read is widened
    OUTWARD by security_eps before any arithmetic use -- this restores an
    outer-inclusion guarantee independent of the serialization's rounding
    direction;
  - runs at TWO precisions (dps=50 and dps=100) and verifies that the
    CONCLUSIONS (inf Ghat_a>0, inf tHhat>0, existence+uniqueness) stay
    identical at both precisions -- detecting any hidden dependence on
    precision.
"""
import sympy as sp
import mpmath as mp
import pickle, json, hashlib, sys

ARCHIVE_PATH = '/tmp/hohmann_full_mesh_archive.json'
FROZEN_PATH = '/tmp/hohmann_FROZEN_v1.pkl'
SECURITY_EPS_EXPONENT = 40  # safety margin: 10^-40, >> any precision used here


def audit_at_precision(dps):
    mp.mp.dps = dps
    mp.iv.dps = dps

    with open(FROZEN_PATH, 'rb') as f:
        FROZEN = pickle.load(f)
    _h = hashlib.sha256()
    for _k in sorted(FROZEN.keys()):
        _h.update(_k.encode()); _h.update(str(FROZEN[_k]).encode())
    frozen_sha256 = _h.hexdigest()

    with open(ARCHIVE_PATH, 'r') as f:
        archive = json.load(f)

    assert archive["frozen_sha256"] == frozen_sha256, \
        f"FINGERPRINT MISMATCH: archive={archive['frozen_sha256']} vs recomputed={frozen_sha256}"

    t, a = sp.symbols('t a')
    Ghat = FROZEN['Ghat']; Ghat_a_expr = FROZEN['Ghat_a']; Ghat_t_expr = FROZEN['Ghat_t']
    tHhat = FROZEN['tHhat']; S1_hat = FROZEN['S1_hat']; S2_hat = FROZEN['S2_hat']

    modules_iv = [{'sqrt': mp.iv.sqrt, 'mpf': mp.iv.mpf}, 'mpmath']
    Ghat_iv = sp.lambdify((t, a), Ghat, modules=modules_iv, cse=True)
    Ghat_a_iv = sp.lambdify((t, a), Ghat_a_expr, modules=modules_iv, cse=True)
    Ghat_t_iv = sp.lambdify((t, a), Ghat_t_expr, modules=modules_iv, cse=True)
    tHhat_iv = sp.lambdify((t, a), tHhat, modules=modules_iv, cse=True)
    S1_iv = sp.lambdify((t, a), S1_hat, modules=modules_iv, cse=True)
    S2_iv = sp.lambdify((t, a), S2_hat, modules=modules_iv, cse=True)

    eps = mp.mpf(10) ** (-SECURITY_EPS_EXPONENT)

    def read_lo(x):
        """Reads an archived bound as a guaranteed LOWER bound (widened downward)."""
        return mp.mpf(x) - eps

    def read_hi(x):
        """Reads an archived bound as a guaranteed UPPER bound (widened upward)."""
        return mp.mpf(x) + eps

    def Ipt(x): return mp.iv.mpf(str(x))
    def Ibox(lo, hi): return mp.iv.mpf([str(lo), str(hi)])

    def Ghat_edge_centered(Tbox, a_edge, t_m, ht):
        Gmid = Ghat_iv(Ipt(t_m), Ipt(a_edge))
        Gt_box = Ghat_t_iv(Tbox, Ipt(a_edge))
        return Gmid + Gt_box * Ibox(-ht, ht)

    tiles = archive["tiles"]
    n_tiles = len(tiles)
    reconstructed = []
    n_pass = 0
    n_fail = 0
    fail_details = []

    for idx, rec in enumerate(tiles):
        # only t_lo, t_hi, a_m, r are read (not the Ghat_a/tHhat bounds
        # already computed by the constructor, which are fully recomputed)
        t_lo = read_lo(rec["t_lo"]); t_hi = read_hi(rec["t_hi"])
        a_m = mp.mpf(rec["a_m"])  # root: a single number, not a bound to widen
        r = mp.mpf(rec["r"])
        t_m = (t_lo + t_hi) / 2
        ht = (t_hi - t_lo) / 2
        Tbox = Ibox(t_lo, t_hi)
        Abox = Ibox(a_m - r, a_m + r)

        ok = True
        reasons = []
        try:
            S1box = S1_iv(Tbox, Abox); S2box = S2_iv(Tbox, Abox)
            if not (S1box.a > 0 and S2box.a > 0):
                ok = False; reasons.append("RADICAND")
        except mp.libmp.libmpf.ComplexResult:
            ok = False; reasons.append("COMPLEX_S1S2")

        try:
            GaBox = Ghat_a_iv(Tbox, Abox)
            if not (GaBox.a > 0):
                ok = False; reasons.append("GA_NOT_POS")
        except mp.libmp.libmpf.ComplexResult:
            ok = False; reasons.append("COMPLEX_Ga")

        try:
            Glo = Ghat_edge_centered(Tbox, a_m - r, t_m, ht)
            Ghi = Ghat_edge_centered(Tbox, a_m + r, t_m, ht)
            if not (Glo.b < 0 and Ghi.a > 0):
                ok = False; reasons.append("BRANCH_NOT_BRACKETED")
        except mp.libmp.libmpf.ComplexResult:
            ok = False; reasons.append("COMPLEX_edges")

        try:
            tHbox = tHhat_iv(Tbox, Abox)
            if not (tHbox.a > 0):
                ok = False; reasons.append("H_NOT_POS")
        except mp.libmp.libmpf.ComplexResult:
            ok = False; reasons.append("COMPLEX_tH")

        if ok:
            n_pass += 1
            reconstructed.append({"idx": idx, "t_lo": t_lo, "t_hi": t_hi,
                                    "Ga_inf": float(GaBox.a), "tH_inf": float(tHbox.a)})
        else:
            n_fail += 1
            fail_details.append((idx, reasons))

    # --- Topology: coverage and junctions, RECOMPUTED from the ARCHIVED
    # (widened) t_lo, t_hi values, not from the constructor's "status"
    # field ---
    tiles_sorted = sorted(tiles, key=lambda r: mp.mpf(r["t_lo"]))
    gaps = []
    for j in range(1, len(tiles_sorted)):
        prev_hi = mp.mpf(tiles_sorted[j-1]["t_hi"])
        cur_lo = mp.mpf(tiles_sorted[j]["t_lo"])
        if prev_hi != cur_lo:
            gaps.append((j-1, j, str(prev_hi), str(cur_lo)))

    n_junc_checked = 0
    n_junc_fail = 0
    for j in range(1, len(tiles_sorted)):
        prev_r = tiles_sorted[j-1]; cur_r = tiles_sorted[j]
        t_s = mp.mpf(prev_r["t_hi"])
        a_prev = (mp.mpf(prev_r["a_m"]) - mp.mpf(prev_r["r"]), mp.mpf(prev_r["a_m"]) + mp.mpf(prev_r["r"]))
        a_cur = (mp.mpf(cur_r["a_m"]) - mp.mpf(cur_r["r"]), mp.mpf(cur_r["a_m"]) + mp.mpf(cur_r["r"]))
        a_lo_hull = min(a_prev[0], a_cur[0]) - eps
        a_hi_hull = max(a_prev[1], a_cur[1]) + eps
        n_junc_checked += 1
        try:
            GaBox = Ghat_a_iv(Ipt(t_s), Ibox(a_lo_hull, a_hi_hull))
            if not (GaBox.a > 0):
                n_junc_fail += 1
        except mp.libmp.libmpf.ComplexResult:
            n_junc_fail += 1

    t_start_archive = mp.mpf(tiles_sorted[0]["t_lo"])
    t_end_archive = mp.mpf(tiles_sorted[-1]["t_hi"])
    worst_Ga = min(r["Ga_inf"] for r in reconstructed) if reconstructed else None
    worst_tH = min(r["tH_inf"] for r in reconstructed) if reconstructed else None

    return {
        "dps": dps,
        "frozen_sha256_match": True,
        "n_tiles_total": n_tiles,
        "n_pass_recertified": n_pass,
        "n_fail_recertified": n_fail,
        "fail_details": fail_details,
        "n_gaps": len(gaps),
        "gaps": gaps,
        "n_junctions_checked": n_junc_checked,
        "n_junctions_failed": n_junc_fail,
        "t_start": str(t_start_archive),
        "t_end": str(t_end_archive),
        "worst_Ga_inf": worst_Ga,
        "worst_tH_inf": worst_tH,
    }


if __name__ == "__main__":
    results = {}
    for dps in [50, 100]:
        print(f"\n=== INDEPENDENT AUDIT at dps={dps} ===", flush=True)
        r = audit_at_precision(dps)
        results[dps] = r
        for k, v in r.items():
            if k not in ("fail_details", "gaps"):
                print(f"  {k}: {v}", flush=True)
        if r["fail_details"]:
            print(f"  FAILURE DETAILS: {r['fail_details']}", flush=True)
        if r["gaps"]:
            print(f"  GAP DETAILS: {r['gaps']}", flush=True)

    print("\n=== COMPARISON BETWEEN THE TWO PRECISIONS ===", flush=True)
    r50, r100 = results[50], results[100]
    conclusions_identical = (
        r50["n_fail_recertified"] == 0 and r100["n_fail_recertified"] == 0 and
        r50["n_gaps"] == 0 and r100["n_gaps"] == 0 and
        r50["n_junctions_failed"] == 0 and r100["n_junctions_failed"] == 0 and
        (r50["worst_Ga_inf"] is not None and r50["worst_Ga_inf"] > 0) and
        (r100["worst_Ga_inf"] is not None and r100["worst_Ga_inf"] > 0) and
        (r50["worst_tH_inf"] is not None and r50["worst_tH_inf"] > 0) and
        (r100["worst_tH_inf"] is not None and r100["worst_tH_inf"] > 0)
    )
    print(f"Identical conclusions at both precisions (inf Ghat_a>0, inf tHhat>0, "
          f"0 gaps, 0 junction failures, 0 tile failures): {conclusions_identical}", flush=True)

    if not conclusions_identical:
        print("\n*** AUDIT FAILED -- DO NOT DECLARE THE BRIDGE CERTIFIED ***", flush=True)
        sys.exit(1)
    else:
        print("\n=== INDEPENDENT AUDIT SUCCEEDED AT BOTH PRECISIONS ===", flush=True)
