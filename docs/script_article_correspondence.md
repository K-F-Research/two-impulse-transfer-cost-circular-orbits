# Script-to-Article Correspondence Tables

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

One table per proof folder, mapping every script-internal variable name to
its article-notation symbol and meaning. See `mathematical_conventions.md`
for full definitions; this document only maps names, it does not redefine
them.

## `circular/existence_uniqueness_minimizer/` (supplementary result, not a numbered lemma of the published manuscript -- the article's own Theorem 1 does not require uniqueness; see this folder's own README.md)

| Script name | Article symbol | Meaning |
|---|---|---|
| `p` | `p` | semi-latus rectum |
| `s` | `s := sqrt(p)` | proof substitution local to this lemma (not the article's `z`) |
| `w` | `w` (local use, distinct role) | auxiliary quantity in the stationary-point closure, scoped to this script |
| `a`, `b`, `q` | (script-local) | intermediate quadratic-form coefficients, no article-level symbol |
| `C` | `C(s)` | two-impulse cost as a function of `s` |
| `Cz` | `dC/ds` | derivative of `C` with respect to `s` (variable named `Cz` for historical reasons; it is a derivative with respect to `s`, not `z`) |

## `circular/hohmann_local_minimum/` (article Theorem 2: local curvature at Hohmann; filed internally as "Theorem 1" -- see `mathematical_conventions.md` numbering note)

| Script/doc name | Article symbol | Meaning |
|---|---|---|
| `R` | `R` | `r2/r1` |
| `t` | `t := sqrt(R)` | |
| `phi` | `phi` | angular deviation from Hohmann |

(This folder contains the proof document only, no standalone verification
script.)

## `circular/global_monotonicity/` (article Theorem 1: global monotonicity; filed internally as "Theorem 2" -- see `mathematical_conventions.md` numbering note)

| Script name | Article symbol | Meaning |
|---|---|---|
| `R`, `phi`, `p` | `R`, `phi`, `p` | as in `mathematical_conventions.md` |
| `c` | `cos(phi)` | |
| `e2` | `e^2` | squared eccentricity of the transfer orbit |
| `D` | `D(p,phi)` | quantity whose derivative's sign localizes the minimizer (Lemma 6.1 of the proof document) |
| `g(p)` | `g(p)` | sign-determining factorized quantity (Lemma 8.3) |
| `p1star`, `p2star` | `p1*`, `p2*` | roots of the localization quadratic |
| `S1`, `S2` | `S1`, `S2` | the two impulse-cost radicands |

## `circular/universal_convexity_z/` (article Theorem 3: universal convexity in z)

| Script/doc name | Article symbol | Meaning |
|---|---|---|
| `R`, `t`, `w`, `phi`, `z`, `rho` | as in `mathematical_conventions.md` | |
| `F` | `F(t,w,rho)` | governing quartic whose physical root is `rho*` |
| `N_rho` | `N_rho(t,w,rho)` | sign-selector quadratic identifying the physical branch |
| `F4` | leading coefficient of `F` | |
| `rho*` / `rho_c` | `rho*(t,w)` | the physical branch (Lemma 1.1 of the proof document) |
| `P` | `P(t,w,rho)` | sign-determining quantity reducing to strict convexity of `V(z)` |
| `Q`, `H` (in `derivation_scripts/` and `certificate_data/`) | **script-internal only, no article symbol** | `Q` and `H` ("auxiliary_H_poly") are this derivation stage's own working names for two source polynomials in the resultant chain that closes the proof's last lemma; they are not renotations of the proof document's own `F`, `N_rho`, `P` -- see `certificate_data/README.md`. **`H` here is unrelated to the curvature `H(p,phi)` of `circular/theorem_4/`** (see `mathematical_conventions.md`, "`H`: two unrelated objects sharing a letter"). |
| `RH`, `RJ` | script-internal | resultants `Res_rho(Q,H)` and `Res_rho(Q,J)` respectively (`J` itself script-internal, built from `H,Q` derivatives) |
| `G` | script-internal | `gcd(RH, RJ)` |
| `A`, `B` | script-internal | cofactors `RH/G`, `RJ/G`; their coprimality (`gcd(A,B)=1`) is the fact underlying the proof document's Lemma 4.2 |
| discriminant (degree 1807) | `Res_w(Big, d(Big)/dw)` in the proof document's own prose | certified via `python-flint`; no separate script name |

**Note**: because `Q, H, RH, RJ, G, A, B` are expensive to recompute (the
`RJ` resultant is the single most costly step of the whole chain -- see
`certificate_data/README.md`), these script-internal variables are kept
as-is rather than renamed. This correspondence table is the authoritative
bridge between this derivation stage's own working names and the proof
document's `F, N_rho, P` notation.

## `circular/theorem_4/` (article Theorem 4: certified curvature sign change, R=2)

| Script/doc name | Article symbol | Meaning |
|---|---|---|
| `R` (fixed to 2), `phi`, `p`, `w := sqrt(p)`, `x := cos(phi)` | as in `mathematical_conventions.md` | |
| `H(p,phi)` | `H(p,phi)` | curvature (Hessian-determinant) quantity, `H := C_phiphi*C_pp - C_pphi^2`; the single, consistent definition used across all six regions of this folder (see `circular/theorem_4/README.md` and `06_global_audit/global_coverage_summary.md`) |
| `P`, `P2`, `Ppp` (`01_lower_endpoint/scripts/`) | components of `H` near the lower endpoint | intermediate quantities built by `build_P.py`, `build_P2.py`, `build_Ppp_chainrule.py`; certified and cross-checked by `certify_Ppp_v2.py` and `crosscheck_Ppp_methods.py` |
| mesh/bridge certificates (`02_lower_bridge/`, `03_central_mesh/`, `05_upper_bridge/`) | certified sign of `H` on each subinterval | produced by `full_mesh_bridge.py`, `pilot_mesh_audit.py`/`pilot_mesh_audit_level2.py`, `mesh_certificate_R2.py`, `hohmann_central_join.py`, and `hohmann_full_mesh.py` respectively; each script name is local to the region it certifies, not a repository-wide symbol |
| Hohmann-endpoint certificates (`04_hohmann_endpoint/`) | certified limit of `H` as `phi -> 0` | produced by `hohmann_certificate_v3.py`, `hohmann_raccord_t1e4.py`, and `hohmann_w_bound_exact.py` |

This folder's six subregions (`01_lower_endpoint/` through `06_global_audit/`)
are the region-by-region decomposition of the full physical domain; the
precise statement, the exact role of each script, and how the six regional
certificates compose into the global result are given in
`circular/theorem_4/README.md`, `circular/theorem_4/THEOREM4_PROOF_MAP.md`,
and `06_global_audit/global_coverage_summary.md`, which take precedence over
this correspondence table for any detail not repeated here.

## `circular/curvature_sign_change_general_R/` (repository's own formalization of a one-sentence remark in article Section 7.3 -- not an article-numbered conjecture; see this folder's own README.md)

| Symbol | Article symbol | Meaning |
|---|---|---|
| `R`, `phi`, `w`, `t`, `z` | as in `mathematical_conventions.md` | this folder explicitly reuses the full `t,w,z` reduction of `universal_convexity_z` plus the `w` of `theorem_4`, stated as such in its own README |

No script is present in this folder (document only).

## Software package (`2i-OTCA-C-v1.0-algorithm/`)

Not covered here: this is an implementation package (api/core/derived/
geometry/hohmann/lambert/rk4/units/surrogate modules) serving a different
purpose (a usable transfer-cost calculator) rather than a proof document,
and is out of scope for this document. A separate pass would be needed if
article-notation consistency with this package's public API is also
required.
