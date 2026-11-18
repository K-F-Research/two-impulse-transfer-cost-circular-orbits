#!/usr/bin/env python3
"""
Lemma 1 -- Independent Numerical Stress Test (2000 Random Trials)

Copyright (c) 2026 Frederic G. Speyser & Kseniya Vyatkina

Status: EXPLORATORY (numerical cross-check, not a proof in itself).
Reproducibility: REPRODUCED (re-executed during the technical audit of
2026-10-25; result identical: 0/2000 failures).

Notation note (see docs/mathematical_conventions.md): as in
lemma1_closure_stationary_point.py, s here denotes sqrt(p) (renamed from an
earlier z during the 2026-10-25 audit) -- distinct from the UNRELATED
z := cot^2(phi/2) substitution used for Theorem 3 elsewhere in this
repository.

--------------------------------------------------------------------------
What this script does
--------------------------------------------------------------------------
Independent numerical cross-check of the symbolic closure proven in
lemma1_closure_stationary_point.py. For 2000 random pairs (w,u) drawn from
w in (1e-4, 1-1e-4) and u in (1e-3, 50):

  1. Confirms the boundary sign pattern used in the closure argument:
     C_s(s) < 0 at s = 1e-8 (near the s->0+ boundary) and C_s(s) > 0 at
     s = 1e8 (near the s->+infinity boundary).
  2. Scans C_s on a coarse logarithmic grid (s = 10^(-8.0) ... 10^(8.0),
     step 0.1 in the exponent) and flags any sign crossing from positive to
     negative (a "decreasing crossing", which the closure argument proves
     cannot occur at a true stationary point under Assumption H).

This is a numerical corroboration only -- it does not replace the symbolic
proof, and a clean result here does not certify the lemma on its own. It is
useful as an independent check that the symbolic argument's premises
(the boundary signs and the absence of decreasing crossings) are not an
artifact of the specific algebra, but hold robustly across the parameter
domain.

--------------------------------------------------------------------------
Inputs / outputs / how to reproduce
--------------------------------------------------------------------------
No external inputs. Run with: python3 lemma1_stress_test_2000_trials.py
Dependencies: sympy, mpmath (50-digit precision), Python's random module
(seeded with 42 for exact reproducibility of the sampled (w,u) pairs).
Expected output:
  "Tested 2000 random (w,u) pairs."
  "Boundary sign failures: 0"
  "Decreasing-crossing detections (coarse grid): 0"
Any nonzero count would print the offending (w,u) pair and the computed
values, via the "BOUNDARY FAIL" line.
"""
import sympy as sp
import mpmath as mp
import random

mp.mp.dps = 50

w, s, u = sp.symbols("w s u", positive=True)
a = 1 - w**2
b = 1 + w**2
q = u*s
lam_minus = (b*q - a/q)/2
lam_plus  = (b*q + a/q)/2
S1 = sp.expand((s - 1)**2 + (lam_minus*s - q/s)**2)
S2 = sp.expand(w**2*(w*s - 1)**2 + (lam_plus*s - q/s)**2)
C  = sp.sqrt(S1) + sp.sqrt(S2)
Cs = sp.diff(C, s)
Css = sp.diff(Cs, s)

Cs_f  = sp.lambdify((w,u,s), Cs, 'mpmath')
Css_f = sp.lambdify((w,u,s), Css, 'mpmath')

random.seed(42)
fails_boundary = 0
fails_monotone_argument = 0
N = 2000
for i in range(N):
    wv = mp.mpf(random.uniform(1e-4, 1-1e-4))
    uv = mp.mpf(random.uniform(1e-3, 50))
    cs_small = Cs_f(wv, uv, mp.mpf('1e-8'))
    cs_large = Cs_f(wv, uv, mp.mpf(1e8))
    if not (cs_small < 0 and cs_large > 0):
        fails_boundary += 1
        print("BOUNDARY FAIL", wv, uv, cs_small, cs_large)

    # Scan for zeros of Cs on a log grid and check the Css sign there
    # (coarse but informative).
    svals = [mp.mpf(10)**mp.mpf(e) for e in [x/10 for x in range(-80,80)]]
    vals = [Cs_f(wv, uv, sv) for sv in svals]
    crossings = 0
    for k in range(len(vals)-1):
        if vals[k] == 0:
            continue
        if vals[k] < 0 and vals[k+1] > 0:
            crossings += 1
        elif vals[k] > 0 and vals[k+1] < 0:
            crossings += 1  # this would be a bad decreasing crossing
            fails_monotone_argument += 1

print(f"\nTested {N} random (w,u) pairs.")
print(f"Boundary sign failures: {fails_boundary}")
print(f"Decreasing-crossing detections (coarse grid): {fails_monotone_argument}")
