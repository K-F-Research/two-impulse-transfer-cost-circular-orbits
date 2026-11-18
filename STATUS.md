# Status

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Two independent fields are tracked for every result in this repository.
**Scientific status** describes the mathematical result itself.
**Reproducibility status** describes whether the computational part (if any)
has been independently re-executed by someone other than the original author.

## Nomenclature

Scientific status:
- `PROVEN -- ANALYTICALLY` -- established by handwritten mathematical proof.
- `PROVEN -- COMPUTATIONALLY CERTIFIED` -- established using exact symbolic
  algebra or certified (rigorous) interval arithmetic, without approximation.
- `OPEN CONJECTURE` -- the general statement has not been proven.
- `SUPPLEMENTARY (not part of the published manuscript)` -- a result proved
  in this repository that is correct and independently reproduced, but is
  neither one of the four numbered theorems of the published manuscript
  (Reference [8]) nor assigned a theorem or lemma number of its own there.

No result is published in this repository unless its reproducibility status
is `REPRODUCED` or `NOT APPLICABLE`.

## Numbering note (read this before matching folder names to the article)

This repository's internal folder names and historical file names for the
first two theorems were set before a later renumbering of the manuscript and
were never renamed to match it. Concretely:

- `circular/hohmann_local_minimum/` is filed internally as "Theorem 1", but
  its content is the article's **Theorem 2** (closed-form local curvature at
  the Hohmann configuration).
- `circular/global_monotonicity/` is filed internally as "Theorem 2", but its
  content is the article's **Theorem 1** (global monotonicity).

Theorem 3 (`circular/universal_convexity_z/`) and Theorem 4
(`circular/theorem_4/`) are numbered consistently with the article. See
`docs/mathematical_conventions.md` ("Flagged cross-document numbering
inconsistency") for the full detail. The table below uses the **article's
own numbering** throughout (Theorem 1 = global monotonicity, Theorem 2 =
local curvature, matching Reference [8] exactly), regardless of which folder
contains the supporting files for a given result.

## Circular case

| # | Result | Scientific status | Reproducibility status |
|---|---|---|---|
| 1 | Theorem 1 -- global monotonicity of V_R on (0,pi), for every R>1 (content lives in `circular/global_monotonicity/`, filed there as "Theorem 2" -- see numbering note above) | PROVEN -- ANALYTICALLY | REPRODUCED |
| 2 | Theorem 2 -- closed-form local curvature coefficient K(r1,r2) at the Hohmann configuration, with remainder O(phi^4) (content lives in `circular/hohmann_local_minimum/`, filed there as "Theorem 1" -- see numbering note above) | PROVEN -- ANALYTICALLY | NOT APPLICABLE (hand proof, no script) |
| 3 | Theorem 3 -- universal convexity in z | PROVEN -- COMPUTATIONALLY CERTIFIED | REPRODUCED (see `circular/universal_convexity_z/certificate_data/`) |
| 4 | Theorem 4, R=2 only -- existence and uniqueness of the curvature inflection on the full physical domain | PROVEN -- COMPUTATIONALLY CERTIFIED: a Krawczyk-certified isolated sign-change point, flanked on each side by a certified interval mesh, with the two domain endpoints handled by an exact symbolic identity (lower endpoint) and a certified asymptotic limit of the correct sign (Hohmann endpoint) (see `circular/theorem_4/README.md` and `circular/theorem_4/06_global_audit/global_coverage_summary.md` for the precise statement and region-by-region coverage) | REPRODUCED |
| 5 | Global existence and uniqueness of the minimizer p*(phi) for every R>1, phi in (0,pi) (`circular/existence_uniqueness_minimizer/`) -- **not a numbered result of the published manuscript**: the article's own Theorem 1 is proved by a stronger, more direct argument that explicitly does not require uniqueness (see the article's abstract and Section 3). This result is true, independently proved, and reproduced, and is kept here because the package's numerical solver (`otca_c.core`) relies on it to justify that root-finding converges to the unique global minimizer rather than to an arbitrary local one. | SUPPLEMENTARY (not part of the published manuscript); scientifically PROVEN -- ANALYTICALLY on its own terms | REPRODUCED (technical audit, 2026-10-25) |
| 6 | General-R extension of Theorem 4 -- single curvature sign change of V_R for every R>1, not only R=2 (`circular/curvature_sign_change_general_R/`) -- **not a numbered conjecture of the published manuscript**: Section 7.3 of the article mentions this only in passing, in a single sentence under "Limitations and future work", with no conjecture number. This repository states it as a precise, well-posed target because the R=2 case (#4 above) provides a fully worked model of what a complete proof would require. | OPEN CONJECTURE (repository's own formalization of an article remark; reduction unconditional and proven; R=2 case itself proven as #4; general R strongly corroborated numerically but not closed) | NOT APPLICABLE |

## Software

| Component | Scientific status | Reproducibility status |
|---|---|---|
| `2i-OTCA-C-v1.0-algorithm/` (2i-OTCA-C) | Implements results #1, #2, #3, and #4 above (circular case only; article numbering), plus a closed-form surrogate valid near the Hohmann transfer | REPRODUCED -- full automated test suite included (`tests/`); three independent validation layers (geometry, Lambert, RK4 propagation) are run on every solved transfer; the JavaScript port is kept in sync with the Python library by a parity test. Full suite (585 cases, including the JS/Python parity tests) was re-executed from the repository root on 2026-11-18: 585 passed, 0 skipped, 0 failed -- see `2i-OTCA-C-v1.0-algorithm/docs_atbd/ATBD_2i-OTCA-C_surrogate_v1.0.docx`, Section 10.6. |
