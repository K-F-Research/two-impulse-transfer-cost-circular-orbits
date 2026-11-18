# Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina
#
# Junction certificate at t = 1e-4 between the local Hohmann-endpoint
# certificate (t in [0, 1e-4], fixed box a in [-0.25, 0.25]) and the first
# tile of the adaptive Hohmann bridge mesh (t in [1e-4, t*]).
#
# The two constructions share the symbolic model frozen in
# hohmann_FROZEN_v1.pkl. The junction is established by set inclusion
# rather than by a two-sided overlap: the local certificate proves
# existence and uniqueness of the root of Ghat(t, .) = 0 on the full box
# a in [-0.25, 0.25] for every t in [0, 1e-4], including the boundary
# t = 1e-4; the mesh's first tile proves existence of a root in an
# a-box strictly contained in that same box, at the same t = 1e-4. By
# uniqueness on the larger box, the two roots coincide.
#
# Four checks are performed:
#   1. exact equality of the t = 1e-4 boundary on both sides;
#   2. strict inclusion of the mesh tile's a-box in the local a-box;
#   3. identity of the frozen symbolic expressions shared by both
#      constructions, verified by SHA-256;
#   4. independent re-verification of the existence/uniqueness
#      certificate on both sides.

import sympy as sp
import mpmath as mp
import pickle, json, hashlib, sys

mp.iv.dps = 50
mp.mp.dps = 60

FROZEN_PATH = '/tmp/hohmann_FROZEN_v1.pkl'
ARCHIVE_PATH = '/tmp/hohmann_full_mesh_archive.json'
T_JUNCTION = mp.mpf('1e-4')
SECURITY_EPS = mp.mpf(10) ** (-40)

passed = []
failed = []


def check(name, cond, detail=""):
    if cond:
        passed.append(name)
        print(f"[PASS] {name}" + (f"  ({detail})" if detail else ""))
    else:
        failed.append(name)
        print(f"[FAIL] {name}" + (f"  ({detail})" if detail else ""))


# ---------------------------------------------------------------------
# 1. Exact equality of the t = 1e-4 boundary.
# ---------------------------------------------------------------------
print("=== 1. Equality of the t = 1e-4 boundary ===")

t0_local = mp.mpf('1e-4')

with open(ARCHIVE_PATH, 'r') as f:
    archive = json.load(f)

tiles = sorted(archive['tiles'], key=lambda tl: mp.mpf(tl['t_lo']))
first_tile = tiles[0]
t_lo_mesh = mp.mpf(first_tile['t_lo'])

check(
    "t0_local == t_lo_mesh within SECURITY_EPS = 1e-40",
    abs(t0_local - t_lo_mesh) < SECURITY_EPS,
    f"t0_local={t0_local}  t_lo_mesh={t_lo_mesh}  diff={abs(t0_local - t_lo_mesh)}"
)

# ---------------------------------------------------------------------
# 2. Strict inclusion of the a-box.
# ---------------------------------------------------------------------
print("\n=== 2. Strict inclusion of the a-box ===")

a_lo_local, a_hi_local = mp.mpf('-0.25'), mp.mpf('0.25')
a_m_mesh = mp.mpf(first_tile["a_m"])
r_mesh = mp.mpf(first_tile["r"])
a_lo_mesh = a_m_mesh - r_mesh
a_hi_mesh = a_m_mesh + r_mesh

check(
    "[a_lo_mesh, a_hi_mesh] strictly contained in [a_lo_local, a_hi_local]",
    (a_lo_mesh > a_lo_local) and (a_hi_mesh < a_hi_local),
    f"mesh=[{a_lo_mesh}, {a_hi_mesh}]  local=[{a_lo_local}, {a_hi_local}]"
)

# ---------------------------------------------------------------------
# 3. Identity of the frozen symbolic expressions (SHA-256).
# ---------------------------------------------------------------------
print("\n=== 3. Identity of the frozen expressions (SHA-256) ===")

with open(FROZEN_PATH, 'rb') as f:
    EXPRS_FROZEN = pickle.load(f)
with open('/tmp/hohmann_regular_exprs.pkl', 'rb') as f:
    EXPRS_LOCAL = pickle.load(f)

# Both pickles encode the same object under different key names; alias
# before computing the common-key set.
EXPRS_LOCAL['Ghat_a'] = EXPRS_LOCAL['tCpphat']

common_keys = sorted(set(EXPRS_LOCAL.keys()) & set(EXPRS_FROZEN.keys()))
check(
    "common keys between hohmann_regular_exprs.pkl and hohmann_FROZEN_v1.pkl non-empty",
    len(common_keys) > 0,
    f"common keys: {common_keys}"
)

all_identical = True
for k in common_keys:
    sha_local = hashlib.sha256(str(EXPRS_LOCAL[k]).encode()).hexdigest()
    sha_frozen = hashlib.sha256(str(EXPRS_FROZEN[k]).encode()).hexdigest()
    sha_ok = sha_local == sha_frozen
    sym_ok = sp.simplify(EXPRS_LOCAL[k] - EXPRS_FROZEN[k]) == 0 if not sha_ok else True
    all_identical = all_identical and (sha_ok or sym_ok)
    print(f"  [{k}] sha_local={sha_local[:16]}...  sha_frozen={sha_frozen[:16]}...  "
          f"identical(SHA)={sha_ok}  identical(symbolic)={sym_ok}")

check("all common expressions are identical (SHA-256 or symbolically)", all_identical)

_h = hashlib.sha256()
for _k in sorted(EXPRS_FROZEN.keys()):
    _h.update(_k.encode())
    _h.update(str(EXPRS_FROZEN[_k]).encode())
frozen_sha256 = _h.hexdigest()
check(
    "global FROZEN fingerprint matches the one archived in hohmann_full_mesh_archive.json",
    archive["frozen_sha256"] == frozen_sha256,
    f"archive={archive['frozen_sha256'][:16]}...  recomputed={frozen_sha256[:16]}..."
)

# ---------------------------------------------------------------------
# 4. Independent re-verification of the existence/uniqueness certificate
#    on both sides.
# ---------------------------------------------------------------------
print("\n=== 4. Independent re-verification of both certificates ===")

t, a = sp.symbols('t a')
Ghat = EXPRS_FROZEN['Ghat']
tCpphat = EXPRS_FROZEN['Ghat_a']
tHhat = EXPRS_FROZEN['tHhat']
S1_hat = EXPRS_FROZEN['S1_hat']
S2_hat = EXPRS_FROZEN['S2_hat']

modules_iv = [{'sqrt': mp.iv.sqrt, 'mpf': mp.iv.mpf}, 'mpmath']
f_Ghat = sp.lambdify((t, a), Ghat, modules=modules_iv, cse=True)
f_tCpphat = sp.lambdify((t, a), tCpphat, modules=modules_iv, cse=True)
f_tHhat = sp.lambdify((t, a), tHhat, modules=modules_iv, cse=True)
f_S1 = sp.lambdify((t, a), S1_hat, modules=modules_iv, cse=True)
f_S2 = sp.lambdify((t, a), S2_hat, modules=modules_iv, cse=True)


def Ipt(x):
    return mp.iv.mpf(str(x))


def Ibox(lo, hi):
    return mp.iv.mpf([str(lo), str(hi)])


# -- local certificate: box [0, 1e-4] x [-0.25, 0.25] --
Bt_local = Ibox(0, t0_local)
Ba_local = Ibox(a_lo_local, a_hi_local)

S1_box = f_S1(Bt_local, Ba_local)
S2_box = f_S2(Bt_local, Ba_local)
check("local: S1, S2 > 0 on the full box", mp.mpf(S1_box.a) > 0 and mp.mpf(S2_box.a) > 0)

tCpp_box = f_tCpphat(Bt_local, Ba_local)
check("local: t*Cpp > 0 on the full box (monotonicity => uniqueness at fixed t)",
      mp.mpf(tCpp_box.a) > 0, f"tCpp in {tCpp_box}")

Gm = f_Ghat(Bt_local, Ipt(a_lo_local))
Gp = f_Ghat(Bt_local, Ipt(a_hi_local))
check("local: Ghat(., a_lo) < 0 for all t in [0, 1e-4]", mp.mpf(Gm.b) < 0, f"{Gm}")
check("local: Ghat(., a_hi) > 0 for all t in [0, 1e-4]", mp.mpf(Gp.a) > 0, f"{Gp}")

tH_box_local = f_tHhat(Bt_local, Ba_local)
tH_full_box_positive = mp.mpf(tH_box_local.a) > 0
print(f"  [info] t*Hhat on the full box [0,1e-4]x[-0.25,0.25] (no subdivision): {tH_box_local}")
if not tH_full_box_positive:
    print("  (expected: the full box loses the t-a correlation and does not yield a")
    print("   positive lower bound; the retained sign certificate subdivides the")
    print("   a-interval, see below. Existence and uniqueness of the root do not")
    print("   depend on this subdivision.)")

n_a = 100
a_edges = [a_lo_local + (a_hi_local - a_lo_local) * mp.mpf(i) / n_a for i in range(n_a + 1)]
tH_sub_lo = None
all_sub_positive = True
for i in range(n_a):
    Ba_sub = Ibox(a_edges[i], a_edges[i + 1])
    tH_sub = f_tHhat(Bt_local, Ba_sub)
    lo_i = mp.mpf(tH_sub.a)
    if lo_i <= 0:
        all_sub_positive = False
    tH_sub_lo = lo_i if tH_sub_lo is None else min(tH_sub_lo, lo_i)
check(
    f"local: t*Hhat > 0 on [0,1e-4]x[-0.25,0.25] after subdivision in a (n_a={n_a})",
    all_sub_positive,
    f"min over the {n_a} sub-boxes = {mp.nstr(tH_sub_lo, 12)}"
)

# -- mesh certificate: first tile --
t_lo_t, t_hi_t = mp.mpf(first_tile["t_lo"]), mp.mpf(first_tile["t_hi"])
Bt_mesh = Ibox(t_lo_t, t_hi_t)
Ba_mesh = Ibox(a_lo_mesh, a_hi_mesh)

S1_box_m = f_S1(Bt_mesh, Ba_mesh)
S2_box_m = f_S2(Bt_mesh, Ba_mesh)
check("mesh (first tile): S1, S2 > 0", mp.mpf(S1_box_m.a) > 0 and mp.mpf(S2_box_m.a) > 0)

tCpp_box_m = f_tCpphat(Bt_mesh, Ba_mesh)
check("mesh (first tile): t*Cpp > 0 (uniqueness at fixed t)",
      mp.mpf(tCpp_box_m.a) > 0, f"{tCpp_box_m}")

Gm_mesh = f_Ghat(Bt_mesh, Ipt(a_lo_mesh))
Gp_mesh = f_Ghat(Bt_mesh, Ipt(a_hi_mesh))
check("mesh (first tile): Ghat(., a_lo) < 0", mp.mpf(Gm_mesh.b) < 0, f"{Gm_mesh}")
check("mesh (first tile): Ghat(., a_hi) > 0", mp.mpf(Gp_mesh.a) > 0, f"{Gp_mesh}")

tH_box_mesh = f_tHhat(Bt_mesh, Ba_mesh)
check("mesh (first tile): t*Hhat > 0", mp.mpf(tH_box_mesh.a) > 0, f"{tH_box_mesh}")

# ---------------------------------------------------------------------
# Summary.
# ---------------------------------------------------------------------
print("\n" + "=" * 60)
print("SUMMARY: junction at t = 1e-4 (local Hohmann certificate <-> mesh)")
print("=" * 60)
print(f"PASS: {len(passed)}")
print(f"FAIL: {len(failed)}")
if failed:
    print("\nUnverified checks:")
    for name in failed:
        print(" -", name)
    sys.exit(1)

print("\nALL REQUIRED CHECKS PASSED.")
print("\nConclusion: the local certificate (t in [0,1e-4], a in [-0.25,0.25])")
print("proves existence and uniqueness of a root of Ghat(t,.) = 0 on the full")
print("box, for every t in the closed domain, including t = 1e-4. The first")
print("mesh tile (t in [1e-4, "
      f"{mp.nstr(t_hi_t, 8)}], a in [{mp.nstr(a_lo_mesh, 8)},{mp.nstr(a_hi_mesh, 8)}])")
print("proves existence of a root in an a-box strictly contained in the local")
print("certificate's box, at the same t = 1e-4. By inclusion, this root is")
print("necessarily the local certificate's unique root: the two constructions")
print("describe the same continuous branch at the boundary t = 1e-4.")
