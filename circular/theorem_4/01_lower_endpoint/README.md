# Theorem 4 (R=2) — Region 01_lower_endpoint

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

**Scope:** the physical stationary branch segment `0 < s <= 0.00015`
(`s := sqrt(p)`), the lower endpoint of the domain covered by `theorem_4`
at R=2. Establishes `H < 0` (strict negativity of the curvature
determinant) throughout this segment.

## Contents

- `lower_endpoint_sign_conclusion.md` — the complete argument:
  exact endpoint identities `P(0)=0`, `P'(0)=0`, the interval-certified
  bound `P''<0` with branch inclusion, and the Taylor-Lagrange step
  concluding `H<0`.
- `scripts/` — all scripts cited by the proof document:
  - `build_P.py`, `build_P2.py` — construct `P := s^6*H` and verify
    `P(0,1/2)=0`.
  - `build_Ypp_Yppp_G.py` — derives `Y'(0)`, `Y''(0)`, `Y'''(0)` by
    implicit differentiation of `G`.
  - `verify_Pprime_symbolic.py` — method A for `P'(0)` and `P''(0)`.
  - `build_Ppp_chainrule.py` — method B (chain rule) for the same second
    derivatives.
  - `crosscheck_Ppp_methods.py` — confronts methods A and B term by
    term.
  - `certify_Ppp_v2.py` — interval certificate for `P''<0` and for
    stationary-branch inclusion in the certified rectangle.
  - `F0_full_pipeline.py` — builds `H` and `F0 := s^4*H` from the
    original cost formulas; numeric corroboration of `F0`'s limit along
    the branch.
- `lower_endpoint_certificate.json` — certification record:
  domain, claims, method distinctions, artifact hashes, dependencies.

## Notation

`s := sqrt(p)` throughout this region's documents. The historical
scripts use `w` for this same quantity (`w = s`), consistent with the
project-wide convention `w := sqrt(p)` in `docs/mathematical_conventions.md`
(shared with `inflection_uniqueness_R2`, `universal_convexity_z`,
`curvature_sign_change_general_R`). This is a **different** quantity
from the `w := R^(-1/2)` used in `theorem_4/05_upper_bridge`'s
certificate. See `lower_endpoint_sign_conclusion.md`, Section 1,
for the full disambiguation.

## How to reproduce

```
python3 scripts/build_P.py
python3 scripts/build_P2.py
python3 scripts/build_Ypp_Yppp_G.py
python3 scripts/verify_Pprime_symbolic.py
python3 scripts/build_Ppp_chainrule.py
python3 scripts/crosscheck_Ppp_methods.py
python3 scripts/certify_Ppp_v2.py
python3 scripts/F0_full_pipeline.py
```

`build_P.py` and `crosscheck_Ppp_methods.py` depend on pickled symbolic
intermediates (`P_built.pkl`, `Q_second_partials.pkl`) produced by
`build_P.py` and `build_Ppp_chainrule.py` respectively; run them first.
`verify_Pprime_symbolic.py` takes several minutes (symbolic
differentiation without intermediate simplification). Dependencies:
`sympy`, `mpmath`.
