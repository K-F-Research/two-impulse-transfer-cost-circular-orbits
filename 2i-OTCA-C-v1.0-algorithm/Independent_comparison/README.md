# Independent Comparison — Lambert p-Iteration Cross-Check

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina.
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See LICENSE-CODE.md at the repository root for the full terms.

This folder contains an independent numerical cross-check of the `otca_c.core` exact
two-impulse transfer solver (2i-OTCA-C v1.0), performed with a second, independently
coded implementation based on the classical semi-latus-rectum (p-iteration) Lambert
formulation of Shen, H. & Tsiotras, P. (AAS/AIAA Astrodynamics Specialist Conference,
1999). The corresponding write-up is reproduced in Section 9.7 of the
Algorithm Theoretical Basis Document (ATBD) for this component.

See `FINAL_REPORT.md` for the full methodology, the four verification controls (A-D),
and the numerical results.

## Quick start

```
pip install mpmath matplotlib
python3 lambert_independent_mp.py      # sanity check vs. closed-form Hohmann limit
python3 vector_reconstruction.py       # Control A
python3 global_min_check.py            # Control B
python3 scaling_check.py               # Control C
python3 k_convergence.py               # independent finite-difference check of K
python3 grid_compare.py                # Run A: 36-case comparison grid
python3 grid_compare_full.py           # Run B: 230-case comparison grid (Annex B.3)
python3 plot_dv_curves.py              # Figure 9.1 (png, svg, pdf)
python3 -m pytest ../tests/test_surrogate.py ../tests/test_surrogate_properties.py -v
```

`grid_compare.py`, `grid_compare_full.py` and `k_convergence.py` import `otca_c` from
the main package (`otca_c.api`, `otca_c.surrogate`); the other scripts are fully
self-contained and do not depend on `otca_c` at all.
