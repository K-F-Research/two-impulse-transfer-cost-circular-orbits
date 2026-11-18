# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
"""
hohmann_full_mesh.py

Full adaptive mesh of the Hohmann bridge, from T_START to T_TARGET (= the
exact t*, the boundary of the already-certified central mesh). Method
validated on 8 test tiles:

    Ghat_a = tCpphat                         : NATURAL evaluation (already
                                                 well conditioned)
    S1_hat, S2_hat, tHhat                    : NATURAL evaluation
    Ghat at the edges a=a_lo,a_hi (boundaries): form CENTERED IN T ONLY
        Ghat(T,a_edge) subseteq Ghat(t_m,a_edge) + Ghat_t(T,a_edge)*[-ht,ht]
        with the central value ITSELF evaluated as a degenerate interval
        (never as an ordinary mp.mp float).

FROZEN expressions (hohmann_FROZEN_v1.pkl, SHA-256 recorded in the log):
do not modify them further -- any future improvement must create a
separate v2.

Convention (re-documented explicitly):
    Ghat(t,a)  := Cp(p_H + a*t, t)
    Ghat_a     := d/da Ghat (t fixed)        = t*Cpp(p(t),t)   = tCpphat
    Ghat_t     := d/dt Ghat (a fixed)        = a*Cpp(p(t),t) + Cpt(p(t),t)
                                              = (a*tCpphat + tCpthat)/t

Failure classification (never interpreted as a false property of the
function -- only as "not established with this box"):
    BRANCH_NOT_BRACKETED   : Ghat signs not opposite at the edges of A
    GA_CONTAINS_ZERO       : Ghat_a not strictly positive on the box
    H_LOWER_NONPOSITIVE    : lower bound of t*Hhat not strictly positive
    RADICAND_UNCERTAIN     : S1 or S2 not strictly positive on the box
    COMPLEX_RESULT         : interval overestimation -> negative radicand
                             (NOT a proof that the function is irregular)

Atomic checkpointing after every top-level segment (resilience to VM
interruptions, already experienced on this project).
"""
import sympy as sp
import mpmath as mp
import pickle, json, time, hashlib, os

mp.mp.dps = 50
mp.iv.dps = 50
IV_DPS = 50
SCRIPT_VERSION = "hohmann_full_mesh.py v1 (2026-11-09)"

with open('/tmp/hohmann_FROZEN_v1.pkl', 'rb') as f:
    FROZEN = pickle.load(f)

_h = hashlib.sha256()
for _k in sorted(FROZEN.keys()):
    _h.update(_k.encode()); _h.update(str(FROZEN[_k]).encode())
FROZEN_SHA256 = _h.hexdigest()

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
Ghat_mp = sp.lambdify((t, a), Ghat, modules='mpmath', cse=True)

# ---------------------------------------------------------------------
# Parameters (CONFIGURABLE, never relaxed silently)
# ---------------------------------------------------------------------
T_START = mp.mpf('1e-4')
T_TARGET = mp.mpf('0.19055605997847746151659843273598244615251833188744')  # exact t*
HT_INIT = mp.mpf('1e-3')         # initial half-width in t per top-level tile
R_CANDIDATES = [mp.mpf(x) for x in ['0.03','0.02','0.01','0.005','0.003','0.002',
                                     '0.001','0.0005','0.0002','0.0001']]
MAX_SUBDIV_DEPTH = 30
CHECKPOINT_PATH = '/tmp/hohmann_full_mesh_checkpoint.json'
ARCHIVE_PATH = '/tmp/hohmann_full_mesh_archive.json'


def Ipt(x): return mp.iv.mpf(str(x))
def Ibox(lo, hi): return mp.iv.mpf([str(lo), str(hi)])


def s(x):
    """Serialize for the JSON archive, with an EXPLICIT SAFETY MARGIN
    toward the outside. Note: mp.nstr() performs ORDINARY decimal
    (nearest) rounding, which does NOT guarantee that the printed lower
    bound stays <= the true lower bound (resp. >= for the upper bound) --
    the archive could therefore lose the outward-rounding property of the
    intervals.
    Fix: many more digits are printed than the actual precision (dps+20),
    so that the string's rounding error is astronomically smaller than
    any interval width handled here (~1e-45), AND it is explicitly
    documented that any audit script reading this archive MUST widen
    every bound read by security_eps = 10^-(dps+10) toward the outside
    before any arithmetic use -- which restores an outward-inclusion
    guarantee independently of nstr()'s rounding direction."""
    extra_digits = mp.mp.dps + 20
    if hasattr(x, '_mpi_'):
        lo = mp.mp.mpf(x._mpi_[0]); hi = mp.mp.mpf(x._mpi_[1])
        if lo == hi:
            return mp.nstr(lo, extra_digits)
        return [mp.nstr(lo, extra_digits), mp.nstr(hi, extra_digits)]
    return mp.nstr(x, extra_digits)


def seed_a(t0_val, hint):
    return mp.findroot(lambda z: Ghat_mp(t0_val, z), hint)


def Ghat_edge_centered(Tbox, a_edge, t_m, ht):
    """Ghat(Tbox,a_edge) via a form centered in t ONLY (a is a fixed
    point here). The central value is evaluated as a DEGENERATE interval
    -- not an ordinary float -- so the result remains a certified bound
    end to end."""
    Gmid = Ghat_iv(Ipt(t_m), Ipt(a_edge))
    Gt_box = Ghat_t_iv(Tbox, Ipt(a_edge))
    return Gmid + Gt_box * Ibox(-ht, ht)


def certify_tile(t_lo, t_hi, hint_a):
    """Attempts to certify the tile [t_lo,t_hi] by trying the
    R_CANDIDATES (from widest to narrowest). Returns (record_or_None,
    reason, a_m) -- a_m is returned EVEN on failure, to seed the next
    tile's root without losing branch tracking."""
    t_m = (t_lo + t_hi) / 2
    ht = (t_hi - t_lo) / 2
    Tbox = Ibox(t_lo, t_hi)
    a_m = seed_a(t_m, hint_a)

    # ALL candidates are tried (no stopping at the first success), and the
    # one giving the LARGEST comfort margin on t*Hhat is kept -- avoiding a
    # razor-thin certificate (e.g. 0.0008) when a smaller r, tried later in
    # the loop, would give a comfortable margin (0.19) at negligible extra
    # cost. Rigor is unchanged: every retained tile STRICTLY satisfies the
    # 4 criteria; only the best among the valid candidates is chosen.
    last_reason = "NO_CANDIDATE_TRIED"
    best_record = None
    best_margin_tH = None
    for r in R_CANDIDATES:
        Abox = Ibox(a_m - r, a_m + r)
        try:
            S1box = S1_iv(Tbox, Abox); S2box = S2_iv(Tbox, Abox)
        except mp.libmp.libmpf.ComplexResult:
            last_reason = "COMPLEX_RESULT_S1S2"; continue
        if not (S1box.a > 0 and S2box.a > 0):
            last_reason = "RADICAND_UNCERTAIN"; continue

        try:
            GaBox = Ghat_a_iv(Tbox, Abox)
        except mp.libmp.libmpf.ComplexResult:
            last_reason = "COMPLEX_RESULT_Ga"; continue
        if not (GaBox.a > 0):
            last_reason = "GA_CONTAINS_ZERO"; continue

        try:
            Glo = Ghat_edge_centered(Tbox, a_m - r, t_m, ht)
            Ghi = Ghat_edge_centered(Tbox, a_m + r, t_m, ht)
        except mp.libmp.libmpf.ComplexResult:
            last_reason = "COMPLEX_RESULT_Ghat_edges"; continue
        if not (Glo.b < 0 and Ghi.a > 0):
            last_reason = "BRANCH_NOT_BRACKETED"; continue

        try:
            tHbox = tHhat_iv(Tbox, Abox)
        except mp.libmp.libmpf.ComplexResult:
            last_reason = "COMPLEX_RESULT_tH"; continue
        if not (tHbox.a > 0):
            last_reason = "H_LOWER_NONPOSITIVE"; continue

        # this candidate r strictly satisfies the 4 criteria: keep it as
        # the best if its tH margin exceeds the best margin seen so far.
        m_G = min(abs(float(Glo.b)), abs(float(Ghi.a)))
        cand_margin_tH = float(tHbox.a)
        if best_margin_tH is None or cand_margin_tH > best_margin_tH:
            best_margin_tH = cand_margin_tH
            best_record = {
                "t_lo": s(t_lo), "t_hi": s(t_hi), "t_m": s(t_m), "ht": s(ht),
                "a_m": s(a_m), "r": s(r),
                "A": [s(a_m - r), s(a_m + r)],
                "Ghat_a": [s(GaBox.a), s(GaBox.b)],
                "Ghat_lo": [s(Glo.a), s(Glo.b)],
                "Ghat_hi": [s(Ghi.a), s(Ghi.b)],
                "tHhat": [s(tHbox.a), s(tHbox.b)],
                "S1": [s(S1box.a), s(S1box.b)],
                "S2": [s(S2box.a), s(S2box.b)],
                "margin_G": m_G,
                "margin_Ga": float(GaBox.a),
                "margin_tH": cand_margin_tH,
            }
        # once a comfortable margin is reached (> COMFORT_TH), no need to
        # try even smaller r -- marginal gain, increasing cost
        if best_margin_tH is not None and best_margin_tH > 0.05:
            break

    if best_record is not None:
        return best_record, None, a_m
    return None, last_reason, a_m


def process_interval(t_lo, t_hi, hint_a, depth, accepted, failures):
    if depth > MAX_SUBDIV_DEPTH:
        failures.append({"t_lo": s(t_lo), "t_hi": s(t_hi), "reason": "MAX_DEPTH_EXCEEDED"})
        raise RuntimeError(f"max subdivision depth exceeded on [{t_lo},{t_hi}]")
    record, reason, a_m = certify_tile(t_lo, t_hi, hint_a)
    if record is not None:
        accepted.append(record)
        return a_m
    failures.append({"t_lo": s(t_lo), "t_hi": s(t_hi), "reason": reason})
    mid = (t_lo + t_hi) / 2
    a_mid = process_interval(t_lo, mid, hint_a, depth + 1, accepted, failures)
    a_right = process_interval(mid, t_hi, a_mid, depth + 1, accepted, failures)
    return a_right


def check_junction(prev_r, cur_r):
    """Junction argument: A_- and A_+ do NOT need to overlap beforehand.
    The proof is:
      - the previous tile guarantees (opposite signs at its edges) a root
        of Ghat(t_s, .) in A_- ;
      - the next tile guarantees, likewise, a root in A_+ ;
      - Ghat_a(t_s, hull(A_-,A_+)) > 0 (STRICT positivity, certified by
        interval arithmetic, with positive radicands included via the
        evaluation itself) gives strict monotonicity of Ghat in a on the
        hull hull(A_-,A_+) ;
      - a strictly monotone function has AT MOST one root on that domain;
      - the two roots, each found in a subset of hull(A_-,A_+), are
        therefore IDENTICAL.
    It is this identity (not mere "proximity") that justifies joining the
    two tiles into a single continuous branch."""
    if mp.mpf(prev_r["t_hi"]) != mp.mpf(cur_r["t_lo"]):
        return {"status": "GAP_OR_OVERLAP_IN_T"}
    t_s = mp.mpf(prev_r["t_hi"])
    a_prev = (mp.mpf(prev_r["A"][0]), mp.mpf(prev_r["A"][1]))
    a_cur = (mp.mpf(cur_r["A"][0]), mp.mpf(cur_r["A"][1]))
    a_lo_hull = min(a_prev[0], a_cur[0]); a_hi_hull = max(a_prev[1], a_cur[1])
    Tpt = Ipt(t_s); Ahull = Ibox(a_lo_hull, a_hi_hull)
    try:
        GaBox = Ghat_a_iv(Tpt, Ahull)
        ok = GaBox.a > 0
    except mp.libmp.libmpf.ComplexResult:
        ok = False
    return {"status": "OK" if ok else "FAIL_Ga_on_hull",
            "t_s": s(t_s), "A_hull": [s(a_lo_hull), s(a_hi_hull)]}


def write_checkpoint(accepted, failures, junctions, next_idx, hint, elapsed):
    ckpt = {
        "script_version": SCRIPT_VERSION, "frozen_sha256": FROZEN_SHA256,
        "next_boundary_idx": next_idx, "hint_a": s(hint),
        "accepted": accepted, "failures": failures, "junctions": junctions,
        "elapsed_seconds_so_far": elapsed,
    }
    tmp_path = CHECKPOINT_PATH + '.tmp'
    with open(tmp_path, 'w') as f:
        json.dump(ckpt, f, indent=1)
    os.replace(tmp_path, CHECKPOINT_PATH)


def load_checkpoint():
    if not os.path.exists(CHECKPOINT_PATH):
        return None
    with open(CHECKPOINT_PATH, 'r') as f:
        return json.load(f)


if __name__ == "__main__":
    t0 = time.time()
    print(f"Frozen expressions, SHA-256 = {FROZEN_SHA256}", flush=True)
    print(f"Bridge domain: T in [{T_START}, {T_TARGET}]  (exact t*)", flush=True)

    boundaries = [T_START]
    tc = T_START
    while tc < T_TARGET:
        tn = min(tc + 2 * HT_INIT, T_TARGET)
        boundaries.append(tn)
        tc = tn
    print(f"Number of top-level segments: {len(boundaries)-1}", flush=True)

    ckpt = load_checkpoint()
    a1 = mp.mpf('-0.0751885346646554796133223180548')
    if ckpt is not None and ckpt.get("frozen_sha256") == FROZEN_SHA256:
        accepted = ckpt["accepted"]; failures = ckpt["failures"]; junctions = ckpt["junctions"]
        start_idx = ckpt["next_boundary_idx"]; hint = mp.mpf(ckpt["hint_a"])
        elapsed_prior = ckpt["elapsed_seconds_so_far"]
        print(f"RESUMING from checkpoint: {len(accepted)} tiles, segment {start_idx}/{len(boundaries)-1}, "
              f"{elapsed_prior:.1f}s already elapsed", flush=True)
    else:
        accepted = []; failures = []; junctions = []
        start_idx = 0; hint = a1; elapsed_prior = 0.0
        if ckpt is not None:
            print("Checkpoint found but fingerprint differs -- ignored, restarting from zero.", flush=True)

    for i in range(start_idx, len(boundaries) - 1):
        t_lo, t_hi = boundaries[i], boundaries[i + 1]
        n_before = len(accepted)
        try:
            hint = process_interval(t_lo, t_hi, hint, 0, accepted, failures)
        except RuntimeError as e:
            print(f"UNRECOVERABLE FAILURE on [{t_lo},{t_hi}]: {e}", flush=True)
        # junctions for any newly accepted tile in this segment
        accepted.sort(key=lambda r: mp.mpf(r["t_lo"]))
        for j in range(max(1, n_before), len(accepted)):
            if j - 1 >= 0 and not any(
                jc.get("t_s") == accepted[j-1]["t_hi"] for jc in junctions
            ):
                jr = check_junction(accepted[j-1], accepted[j])
                jr["between_idx"] = [j-1, j]
                junctions.append(jr)
                if jr["status"] != "OK":
                    print(f"  [JUNCTION FAILURE] between tiles {j-1},{j}: {jr}", flush=True)
        write_checkpoint(accepted, failures, junctions, i + 1, hint, elapsed_prior + (time.time() - t0))
        n_fail_junc = sum(1 for jc in junctions if jc["status"] != "OK")
        print(f"[checkpoint segment {i+1}/{len(boundaries)-1}] "
              f"tiles={len(accepted)} junction_failures={n_fail_junc} "
              f"t_hi_reached={t_hi}", flush=True)

    write_checkpoint(accepted, failures, junctions, len(boundaries) - 1, hint,
                      elapsed_prior + (time.time() - t0))

    accepted.sort(key=lambda r: mp.mpf(r["t_lo"]))

    # --- BLOCKING assertions: an incomplete result must never be
    # presented as a success. ---
    assert accepted, "No tile accepted -- total failure"
    assert mp.mpf(accepted[0]["t_lo"]) == T_START, \
        f"First tile does not start at T_START: {accepted[0]['t_lo']} != {T_START}"
    assert mp.mpf(accepted[-1]["t_hi"]) == T_TARGET, \
        f"Last tile does not reach T_TARGET: {accepted[-1]['t_hi']} != {T_TARGET}"
    for j in range(1, len(accepted)):
        assert mp.mpf(accepted[j-1]["t_hi"]) == mp.mpf(accepted[j]["t_lo"]), \
            f"Gap or overlap between tiles {j-1} and {j}"
    assert len(junctions) == len(accepted) - 1, \
        f"Number of junctions ({len(junctions)}) != number of tiles - 1 ({len(accepted)-1})"
    assert all(jc["status"] == "OK" for jc in junctions), \
        "At least one junction is not certified OK"
    print("\n=== ALL BLOCKING COVERAGE/JUNCTION ASSERTIONS PASS ===", flush=True)

    print(f"\nTiles accepted: {len(accepted)}", flush=True)
    print(f"Failures (including recovered subdivisions): {len(failures)}", flush=True)
    n_junc_fail = sum(1 for jc in junctions if jc["status"] != "OK")
    print(f"Junctions tested: {len(junctions)}, failures: {n_junc_fail}", flush=True)
    if accepted:
        worst_tH = min(mp.mpf(r["tHhat"][0]) for r in accepted)
        worst_Ga = min(mp.mpf(r["Ghat_a"][0]) for r in accepted)
        print(f"Worst lower bound of t*Hhat over the whole bridge: {worst_tH}", flush=True)
        print(f"Worst lower bound of Ghat_a over the whole bridge : {worst_Ga}", flush=True)
        # coverage: check that there is no gap between consecutive tiles
        gaps = []
        for j in range(1, len(accepted)):
            if mp.mpf(accepted[j-1]["t_hi"]) != mp.mpf(accepted[j]["t_lo"]):
                gaps.append((accepted[j-1]["t_hi"], accepted[j]["t_lo"]))
        print(f"Coverage gaps detected: {len(gaps)}", flush=True)
        print(f"t_lo of the first tile: {accepted[0]['t_lo']}", flush=True)
        print(f"t_hi of the last tile: {accepted[-1]['t_hi']}", flush=True)

    archive = {
        "script_version": SCRIPT_VERSION, "frozen_sha256": FROZEN_SHA256,
        "mp_dps": mp.mp.dps, "iv_dps": IV_DPS,
        "domain_T": [s(T_START), s(T_TARGET)],
        "n_tiles": len(accepted), "n_failures_total_incl_recovered": len(failures),
        "n_junctions": len(junctions), "n_junction_failures": n_junc_fail,
        "tiles": accepted, "failures": failures, "junctions": junctions,
        "elapsed_seconds": elapsed_prior + (time.time() - t0),
    }
    with open(ARCHIVE_PATH, 'w') as f:
        json.dump(archive, f, indent=1)
    print(f"\nFinal archive written: {ARCHIVE_PATH}", flush=True)
