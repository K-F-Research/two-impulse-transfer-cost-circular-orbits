# Independent Numerical Cross-Check Report — Lambert p-Iteration Method

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina.
Licensed under the Research, Education and Non-Commercial Use License (NCSL-1.0).
See LICENSE-CODE.md at the repository root for the full terms.

## Run identification (to remove any ambiguity)

- **Run A ("36 cases")**: R in {1.2, 1.5, 2.0, 3.0, 5.0, 10.0} x deviation in {5,10,15,20,25,30} deg.
  Files: `grid_compare.py` / `grid_compare_results.json`. Max relative difference = 5.812e-14.
- **Run B ("230 cases", aligned with ATBD Annex B.3)**: R in {1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.6,1.75,1.9,2.0,2.25,2.5,2.75,3.0,3.5,4.0,4.5,5.0,7.0,10.0,20.0} x deviation in {1,2,3,5,7,10,15,20,25,30} deg.
  Files: `grid_compare_full.py` / `grid_compare_full_results.json`. Max relative difference = 5.747e-12, at case R=1.02, dev=3 deg.

The previously reported worst case (R=1.02) belongs exclusively to Run B, not Run A. Run A does not contain any R below 1.2.

## Additional verification controls performed (before ATBD integration)

### A. Independent vector reconstruction of the impulses
Script: `vector_reconstruction.py`. Computed via explicit radial/tangential component
subtraction, without reusing the scalar law-of-cosines formula (Shen & Tsiotras Eq. 4).
Result: EXACT agreement (0 difference at 50 significant digits) on 4 representative cases.
Confirms there is no algebraic transcription error in the scalar formula used.

### B. Global-minimum verification
Script: `global_min_check.py`. Dense scan (4000, then 40000 points near the limiting case)
of the full admissible interval for 5 representative cases, including a near-degenerate
case (R=1.02) and a large-radius-ratio case (R=20).
Result: a single local minimum in every case (unimodal function on the short-path branch),
no minimum near the admissible-interval boundaries. The R=1.02 case initially showed a
1.27e-3 gap between the golden-section minimum and the 4000-point scan minimum; refining
the scan to 40000 points reduced this gap to 1.2e-5, confirming it was a discrete-scan
resolution limitation in an ill-conditioned region, not a second minimum nor a defect of
the golden-section search.

### C. Scaling laws (tested on the independent script itself, not on otca_c)
Script: `scaling_check.py`.
- Length scaling (r1,r2 -> lambda*r1,lambda*r2): dV scales as 1/sqrt(lambda).
  Verified for lambda in {0.25,0.5,2,4,10}: zero difference or at floating-point noise level (1e-16).
- Gravitational-parameter scaling (mu -> q*mu): dV scales as sqrt(q).
  Verified for q in {0.25,0.5,2,4,10}: zero difference across the whole range.

### D. Code inspection (f1/f2 angles)
Verified directly in `lambert_independent_mp.py`: `speed_and_fpa(r1, f1)` and
`speed_and_fpa(r2, f2)` are called separately with distinct f1 and f2=f1+theta.
No erroneous reuse of f1 at both points in the CODE.

## Explicit scope limitations

- Only the short-path (single-revolution) branch was implemented and tested, consistent
  with the phi = pi - delta_theta convention already used throughout the ATBD, and must
  be stated explicitly there.
- This is an independent-method numerical cross-check (Lambert/p-iteration formulation vs.
  the internal otca_c.core implementation), NOT a physical validation against flight data.
  An agreement at 1e-12 to 1e-14 characterizes the numerical consistency of two ways of
  solving the same idealized two-body problem; it says nothing about how well that
  idealized problem represents a real mission's physics.
- Citation: Shen, H. & Tsiotras, P., "Optimal Two-Impulse Rendezvous Between Two Circular
  Orbits Using Multiple-Revolution Lambert's Solutions," AAS/AIAA Astrodynamics Specialist
  Conference, Girdwood, Alaska, August 16-18, 1999 (conference version; Eqs. (3)-(4) and
  the free-time auxiliary-transfer-problem Appendix are the ones implemented here). Not to
  be conflated with the later journal version (Journal of Guidance, Control, and Dynamics,
  Vol. 26, No. 1, 2003, pp. 50-61), which additionally treats a fixed-time rendezvous
  problem not used here.

## Files in this folder

| File | Role |
|---|---|
| `lambert_independent_mp.py` | Core independent method (p-iteration, mpmath 50-digit precision) |
| `vector_reconstruction.py` | Control A |
| `global_min_check.py` | Control B |
| `scaling_check.py` | Control C |
| `k_convergence.py` | Independent finite-difference check of curvature coefficient K |
| `grid_compare.py` / `grid_compare_results.json` | Run A (36 cases) |
| `grid_compare_full.py` / `grid_compare_full_results.json` | Run B (230 cases, matches ATBD Annex B.3) |
| `k_convergence_results.json` | K convergence sweep results (4 R values x 10 phi steps) |
| `plot_dv_curves.py` / `dv_curves_global_min_check.png` | Figure for Control B |
